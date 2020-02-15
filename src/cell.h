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
  Params() : splitProb(.5), stackEnergy(1), auEnergy(-.5), gcEnergy(.5), guEnergy(-1), temp(.5) { }
  static Params fromJson (json&);
  json toJson() const;
};

struct Board {
  typedef pair<int,int> IndexPair;
  vguard<int> cellStorage;
  vguard<Vec> neighborhood;
  uniform_real_distribution<> dist;
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
    return !(u.next == v.index || v.next == u.index)
      && u.isComplementOrWobble(v)
      && u.next != v.index  // disallow neighbors
      && v.next != u.index
      && u.next != v.prev  // disallow next-but-one neighbors
      && u.prev != v.next
      && !indicesPaired (u.prev, v.prev)  // disallow parallel stacking
      && !indicesPaired (u.next, v.next);
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

  vguard<IndexPair> indexPairs() const;
  string foldString() const;
  double foldEnergy() const;
};

#endif /* CELL_INCLUDED */
