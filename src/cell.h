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
};

struct Unit {
  char base;
  Vec pos;
  bool isRev;
  int next, prev;
  Unit() { }
  Unit (char b, int x, int y, int z, bool r, int n, int p)
    : base(b), pos(x,y,z), isRev(r), next(n), prev(p)
  { }
};

struct Params {
  double splitProb;  // probability that a move is a split
  Params() : splitProb(.1) { }
  static Params fromJson (json&);
  json toJson() const;
};

class Board {
private:
  vguard<int> cellStorage;
  vguard<Vec> neighborhood;
protected:
  inline int mod (int val, int size) const {
    const int m = val % size;
    return m < 0 ? (m + size) : m;
  }
  inline int cellIndex (int x, int y, int z, bool rev) const {
    return (rev ? 1 : 0) + 2 * (mod(x,xSize) + xSize * (mod(y,ySize) + ySize * mod(z,zSize)));
  }
public:
  int xSize, ySize, zSize;
  Params params;
  vguard<Unit> unit;
  
  Board (int, int, int);
  Board();

  static Board fromJson (json&);
  json toJson() const;

  const Vec& rndNbrVec (mt19937&) const;

  bool cellIsValid (int, int, int) const;

  bool tryMove (const Vec&, const Vec&);
  bool trySplit (const Vec&, const Vec&);
  
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
};

#endif /* CELL_INCLUDED */
