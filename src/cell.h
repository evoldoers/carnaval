#ifndef CELL_INCLUDED
#define CELL_INCLUDED

#include <random>
#include "json.hpp"
#include "vguard.h"

using namespace std;
using json = nlohmann::json;

struct Vec {
  int xyz[3];
  Vec() { xyz[0] = xyz[1] = xyz[2] = 0; }
  Vec (int x, int y, int z) { xyz[0] = x; xyz[1] = y; xyz[2] = z; }
  inline int x() const { return xyz[0]; }
  inline int& x() { return xyz[0]; }
  inline int y() const { return xyz[1]; }
  inline int& y() { return xyz[1]; }
  inline int z() const { return xyz[2]; }
  inline int& z() { return xyz[2]; }
  inline Vec operator+ (const Vec& d) const { return Vec (x()+d.x(), y()+d.y(), z()+d.z()); }
  inline Vec operator- (const Vec& d) const { return Vec (x()-d.x(), y()-d.y(), z()-d.z()); }
  inline bool isZero() const { return x() == 0 && y() == 0 && z() == 0; }
  friend ostream& operator<< (ostream& out, const Vec& v) { return out << "(" << v.x() << "," << v.y() << "," << v.z() << ")"; }
};

struct Unit {
  int base;
  Vec pos;
  bool rev;
  int index, prev, next;
  static string alphabet;  // acgu
  Unit() { }
  Unit (char b, int x, int y, int z, bool r, int i, int p, int n)
    : base(b), pos(x,y,z), rev(r), index(i), prev(p), next(n)
  { }
  inline char baseChar() const { return base2char (base); }
  inline static bool isRNA (char c) {
    return alphabet.find (c) < alphabet.size();
  }
  inline static int char2base (char c) {
    if (!isRNA (c))
      throw runtime_error ("Not a base");
    return alphabet.find(c);
  }
  inline static char base2char (int b) {
    return alphabet[b];
  }
  inline bool isComplementOrWobble (const Unit& u) const {
    const int x = base, y = u.base;
    return (x == 3 - y) || (x * y == 6);
  }
};

struct Params {
  double splitProb;  // probability that a move is a split, given that the Unit is paired
  double stackEnergy, auEnergy, gcEnergy, guEnergy, temp;  // simplified basepair stacking model
  double bondProb;  // probability that two adjacent template-bound monomers will form a covalent bond
  Params() : splitProb(.5), stackEnergy(4), auEnergy(-2), gcEnergy(2), guEnergy(-3), temp(1), bondProb(.01) { }
  static Params fromJson (json&);
  json toJson() const;
};

struct Board {
  typedef pair<int,int> IndexPair;
  vguard<int> cellStorage;
  vguard<Vec> neighborhood;
  uniform_real_distribution<> dist;  // real distribution over [0,1)
  uniform_int_distribution<> baseDist;  // integer distribution over [0,4)
  static string leftFoldChar, rightFoldChar;
  inline static int boardCoord (int val, int size) {
    const int m = val % size;
    return m < 0 ? (m + size) : m;
  }
  inline bool boardCoordsEqual (const Vec& a, const Vec& b) const {
    return boardCoord (a.x(), xSize) == boardCoord (b.x(), xSize)
      && boardCoord (a.y(), ySize) == boardCoord (b.y(), ySize)
      && boardCoord (a.z(), zSize) == boardCoord (b.z(), zSize);
  }
  inline int cellIndex (int x, int y, int z, bool rev) const {
    return (rev ? 1 : 0) + 2 * (boardCoord(x,xSize) + xSize * (boardCoord(y,ySize) + ySize * boardCoord(z,zSize)));
  }
  inline static int nbrRange (int size) {
    return size > 1 ? 1 : 0;
  }
  int xSize, ySize, zSize;
  Params params;
  vguard<Unit> unit;
  
  Board (int, int, int);
  Board();

  static Board fromJson (json&);
  json toJson() const;

  void addSeq (const string&);  // adds sequence along x-axis starting at origin
  void addBases (double, mt19937&);  // adds random monomeric bases with given density
  
  const Vec& rndNbrVec (mt19937&) const;

  void assertValid() const;

  inline static int shortestDistance (int c1, int c2, int size) {
    const int d = boardCoord (c1 - c2, size);
    return min (d, size - d);
  }
  inline static bool coordAdjacent (int c1, int c2, int size) {
    return shortestDistance (c1, c2, size) <= 1;
  }
  inline bool adjacent (const Vec& a, const Vec& b) const {
    return coordAdjacent (a.x(), b.x(), xSize)
      && coordAdjacent (a.y(), b.y(), ySize)
      && coordAdjacent (a.z(), b.z(), zSize);
  }
  inline bool canMoveTo (const Unit& u, const Vec& newPos) const {
    return (u.next < 0 || adjacent (unit[u.next].pos, newPos))
      && (u.prev < 0 || adjacent (unit[u.prev].pos, newPos));
  }
  inline bool isPaired (const Unit& u) const {
    return pairedIndex(u) >= 0;
  }
  inline int pairedIndex (const Unit& u) const {
    return cell (u.pos, !u.rev);
  }
  inline bool indicesPaired (int i, int j) const {
    return i >= 0 && j >= 0 && boardCoordsEqual (unit[i].pos, unit[j].pos);
  }
  inline bool canMerge (const Unit& u, const Unit& v) const {
    if (!u.isComplementOrWobble(v))
      return false;
    if (u.next == v.index || v.next == u.index)  // disallow neighbors
      return false;
    const int u_next2 = u.next >= 0 ? unit[u.next].next : -1;
    const int u_prev2 = u.prev >= 0 ? unit[u.prev].prev : -1;
    if (u_next2 >= 0
	&& (u_next2 == v.index  // disallow next-but-one neighbors
	    || u_next2 == v.prev))  // disallow next-but-two neighbors
      return false;
    if (u_prev2 >= 0
	&& (u_prev2 == v.index
	    || u_prev2 == v.next))
      return false;
  return (!indicesPaired (u.prev, v.prev)  // disallow parallel stacking
	  && !indicesPaired (u.next, v.next));
  }
  inline double calcEnergy (const Unit& u, const Unit& v, double stackWeight) const {
    double e = 0;
    const int bprod = u.base * v.base;
    switch (bprod) {
    case 0: e += params.auEnergy; break;
    case 2: e += params.gcEnergy; break;
    case 6: e += params.guEnergy; break;
    default: throw runtime_error("Not a basepair"); break;
    }
    if (indicesPaired (u.prev, v.next))
      e += params.stackEnergy * stackWeight;
    if (indicesPaired (u.next, v.prev))
      e += params.stackEnergy * stackWeight;
    return e;
  }
  inline double pairingEnergy (const Unit& u, const Unit& v) const {
    return calcEnergy (u, v, 1);
  }
  inline bool acceptMove (double energyDelta, double fwdBackRatio, mt19937& mt) {
    const double p = exp (energyDelta / params.temp) / fwdBackRatio;
    //    cerr << "delta=" << energyDelta << " ratio=" << fwdBackRatio << " p=" << p << endl;
    return p >= 1 || dist(mt) < p;
  }
  inline void moveUnit (Unit& u, const Vec& pos, bool rev) {
    //    cerr << "before move..." << endl; dump(cerr);
    cell (u.pos, u.rev) = -1;
    u.pos = pos;
    u.rev = rev;
    cell (u.pos, u.rev) = u.index;
    //    cerr << u.pos << "." << u.rev << endl;
    //    cerr << "after move..." << endl; dump(cerr);
  }
  
  bool tryMove (mt19937&);
  void dump (ostream&) const;
  
  inline const int& cell (int x, int y, int z, bool rev) const {
    return cellStorage[cellIndex (x, y, z, rev)];
  }
  inline int& cell (int x, int y, int z, bool rev) {
    return cellStorage[cellIndex (x, y, z, rev)];
  }
  inline const int& cell (const Vec& v, bool rev) const {
    return cell (v.x(), v.y(), v.z(), rev);
  }
  inline int& cell (const Vec& v, bool rev) {
    return cell (v.x(), v.y(), v.z(), rev);
  }

  // single-chain logging
  void assertLinear() const;
  vguard<IndexPair> indexPairs() const;
  string sequence() const;
  vguard<Vec> unitPos() const;
  vguard<double> unitCentroid() const;
  double unitRadiusOfGyration() const;
  string foldString() const;
  double foldEnergy() const;
  string coloredFoldString() const;

  // multi-chain logging
  map<string,int> sequenceFreqs() const;
};

#endif /* CELL_INCLUDED */
