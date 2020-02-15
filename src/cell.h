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
  inline bool operator== (const Vec& v) const { return x() == v.x() && y() == v.y() && z() == v.z(); }
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
  Params() : splitProb(.1) { }
  static Params fromJson (json&);
  json toJson() const;
};

class Board {
private:
  vguard<int> cellStorage;
  vguard<Vec> neighborhood;
  uniform_real_distribution<> dist;
protected:
  inline static int boardCoord (int val, int size) {
    const int m = val % size;
    return m < 0 ? (m + size) : m;
  }
  inline int cellIndex (int x, int y, int z, bool rev) const {
    return (rev ? 1 : 0) + 2 * (boardCoord(x,xSize) + xSize * (boardCoord(y,ySize) + ySize * boardCoord(z,zSize)));
  }
  inline static int nbrRange (int size) {
    return size > 1 ? 1 : 0;
  }
public:
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
    return i >= 0 && j >= 0 && unit[i].pos == unit[j].pos;
  }
  inline bool canMerge (const Unit& u, const Unit& v, mt19937& mt) {
    return !(u.next == v.index || v.next == u.index)
      && u.isComplementOrWobble(v)
      && !indicesPaired (u.prev, v.prev)  // disallow parallel stacking
      && !indicesPaired (u.next, v.next)  // disallow parallel stacking
      && u.next != v.index  // disallow neighbors
      && v.next != u.index
      && u.next != v.prev  // disallow next-but-one neighbors
      && u.prev != v.next;
  }
  inline void moveUnit (Unit& u, const Vec& pos, bool rev) {
    //    cerr << "before move..." << endl; dump(cerr);
    cell (u.pos, u.rev) = -1;
    u.pos.x() = boardCoord (pos.x(), xSize);
    u.pos.y() = boardCoord (pos.y(), ySize);
    u.pos.z() = boardCoord (pos.z(), zSize);
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

  string foldString() const;
};

#endif /* CELL_INCLUDED */
