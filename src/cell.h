#ifndef CELL_INCLUDED
#define CELL_INCLUDED

#include <cstdint>
#include <random>
#include <string>
#include <map>
#include "json.hpp"
#include "vguard.h"
#include "rng.h"

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
  uint32_t mask;   // bitmask over neighbourhood directions that keep both chain links adjacent
  bool pinned;     // never diffuses
  bool inert;      // never forms a basepair; used for surface sites and tether anchors
  static string alphabet;  // acgu
  static bool dnaMode;     // accept and emit 't' in place of 'u'
  Unit()
    : base(0), rev(false), index(-1), prev(-1), next(-1), mask(0), pinned(false), inert(false)
  { }
  Unit (int b, int x, int y, int z, bool r, int i, int p, int n)
    : base(b), pos(x,y,z), rev(r), index(i), prev(p), next(n), mask(0), pinned(false), inert(false)
  { }
  inline char baseChar() const { return inert ? '#' : base2char (base); }
  inline static bool isRNA (char c) {
    return alphabet.find (c) < alphabet.size() || c == 't';
  }
  inline static int char2base (char c) {
    if (c == 't') return 3;  // DNA thymine occupies the uracil slot
    if (!isRNA (c))
      throw runtime_error ("Not a base");
    return alphabet.find(c);
  }
  inline static char base2char (int b) {
    return (dnaMode && b == 3) ? 't' : alphabet[b];
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
  double secondsPerSweep;  // experimental time per Monte Carlo sweep; 0 means uncalibrated
  bool allowWobble;  // G-U (G-T) wobble pairs permitted; false for DNA
  Params()
    : splitProb(.5), stackEnergy(4), auEnergy(-2), gcEnergy(2), guEnergy(-3), temp(1),
      bondProb(0), secondsPerSweep(0), allowWobble(true)
  { }
  static Params fromJson (json&);
  json toJson() const;
};

struct Board {
  typedef pair<int,int> IndexPair;
  vguard<int> cellStorage;
  vguard<Vec> neighborhood;
  vguard<int> mobile;        // indices of units that are free to diffuse
  vguard<uint32_t> adjMask;  // indexed by offsetCode: directions adjacent to a given chain-link offset
  uint32_t allMask;          // all neighbourhood directions set
  uniform_real_distribution<> dist;  // real distribution over [0,1)
  uniform_int_distribution<> baseDist;  // integer distribution over [0,4)
  static string leftFoldChar, rightFoldChar;

  int xSize, ySize, zSize;
  int xWrap, yWrap, zWrap;   // size-1 when the size is a power of two, else -1 meaning "use modulo"
  bool legalMovesOnly;       // sample only from geometrically legal directions (alters kinetics)
  double attempts;           // virtual attempt count, on the uniform-proposal clock
  Params params;
  vguard<Unit> unit;

  Board (int, int, int);
  Board();

  static Board fromJson (json&);
  json toJson() const;

  void addSeq (const string&);               // adds sequence along x-axis starting at origin
  void addBases (double, Rng&);              // adds random monomeric bases with given density
  void addSurface (int axis);                // fills the plane at coord 0 normal to axis with inert sites
  int addTetheredArray (const string&, int spacing, int axis);  // returns number of strands placed
  void sortUnitsByPosition();                // renumbers units for cache locality; relabelling only

  void rebuildMobile();
  void rebuildMasks();
  void buildAdjMask();

  void assertValid() const;

  // ---- coordinate helpers -------------------------------------------------
  // A power-of-two size lets the wrap be a single mask instead of an integer
  // division; boardCoord sits on the hottest path in the simulator.
  inline static int wrapOf (int size) { return ((size & (size - 1)) == 0) ? (size - 1) : -1; }
  inline static int boardCoord (int val, int size, int wrap) {
    if (wrap >= 0) return val & wrap;
    const int m = val % size;
    return m < 0 ? (m + size) : m;
  }
  inline static int boardCoord (int val, int size) {
    return boardCoord (val, size, wrapOf (size));
  }
  inline int bcx (int v) const { return boardCoord (v, xSize, xWrap); }
  inline int bcy (int v) const { return boardCoord (v, ySize, yWrap); }
  inline int bcz (int v) const { return boardCoord (v, zSize, zWrap); }

  inline bool boardCoordsEqual (const Vec& a, const Vec& b) const {
    return bcx (a.x()) == bcx (b.x()) && bcy (a.y()) == bcy (b.y()) && bcz (a.z()) == bcz (b.z());
  }
  inline int cellIndex (int x, int y, int z, bool rev) const {
    return (rev ? 1 : 0) + 2 * (bcx(x) + xSize * (bcy(y) + ySize * bcz(z)));
  }
  inline int cellIndex (const Vec& v) const { return cellIndex (v.x(), v.y(), v.z(), false); }
  inline static int nbrRange (int size) { return size > 1 ? 1 : 0; }

  inline static int shortestDistance (int c1, int c2, int size, int wrap) {
    const int d = boardCoord (c1 - c2, size, wrap);
    return min (d, size - d);
  }
  inline static int shortestDistance (int c1, int c2, int size) {
    return shortestDistance (c1, c2, size, wrapOf (size));
  }
  inline bool adjacent (const Vec& a, const Vec& b) const {
    return shortestDistance (a.x(), b.x(), xSize, xWrap) <= 1
      && shortestDistance (a.y(), b.y(), ySize, yWrap) <= 1
      && shortestDistance (a.z(), b.z(), zSize, zWrap) <= 1;
  }

  // signed shortest step from `from` to `to`, in {-1,0,1} for adjacent coords
  inline static int signedStep (int from, int to, int size, int wrap) {
    if (size <= 1) return 0;
    const int d = boardCoord (to - from, size, wrap);
    return (d * 2 > size) ? (d - size) : d;
  }
  // offsets are per-axis values in {-1,0,1}; code them into [0,27)
  inline static int offsetCode (int dx, int dy, int dz) {
    return (dx + 1) + 3 * ((dy + 1) + 3 * (dz + 1));
  }
  inline int offsetCodeTo (const Vec& from, const Vec& to) const {
    int dx = signedStep (from.x(), to.x(), xSize, xWrap);
    int dy = signedStep (from.y(), to.y(), ySize, yWrap);
    int dz = signedStep (from.z(), to.z(), zSize, zWrap);
    if (dx < -1) dx = -1; if (dx > 1) dx = 1;
    if (dy < -1) dy = -1; if (dy > 1) dy = 1;
    if (dz < -1) dz = -1; if (dz > 1) dz = 1;
    return offsetCode (dx, dy, dz);
  }

  // The legal-direction mask for a unit depends only on where its chain
  // neighbours are, so it is two table lookups rather than 8 (or 26) geometry
  // tests.  It is recomputed only for the handful of units a move can affect.
  inline uint32_t computeMask (const Unit& u) const {
    uint32_t m = allMask;
    if (u.prev >= 0) m &= adjMask[offsetCodeTo (u.pos, unit[u.prev].pos)];
    if (u.next >= 0) m &= adjMask[offsetCodeTo (u.pos, unit[u.next].pos)];
    return m;
  }
  inline uint32_t maskAt (const Unit& u, const Vec& at) const {
    uint32_t m = allMask;
    if (u.prev >= 0) m &= adjMask[offsetCodeTo (at, unit[u.prev].pos)];
    if (u.next >= 0) m &= adjMask[offsetCodeTo (at, unit[u.next].pos)];
    return m;
  }
  inline void refreshMask (int i) {
    if (i >= 0) unit[i].mask = computeMask (unit[i]);
  }
  inline void refreshAround (const Unit& u) {
    refreshMask (u.index);
    refreshMask (u.prev);
    refreshMask (u.next);
  }

  inline bool canMoveTo (const Unit& u, const Vec& newPos) const {
    return (u.next < 0 || adjacent (unit[u.next].pos, newPos))
      && (u.prev < 0 || adjacent (unit[u.prev].pos, newPos));
  }
  inline bool isPaired (const Unit& u) const { return pairedIndex(u) >= 0; }
  inline int pairedIndex (const Unit& u) const { return cell (u.pos, !u.rev); }
  inline bool indicesPaired (int i, int j) const {
    return i >= 0 && j >= 0 && boardCoordsEqual (unit[i].pos, unit[j].pos);
  }
  inline bool canMerge (const Unit& u, const Unit& v) const {
    if (u.inert || v.inert)
      return false;
    if (!u.isComplementOrWobble(v))
      return false;
    if (!params.allowWobble && u.base * v.base == 6)
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
  inline bool acceptMove (double energyDelta, double fwdBackRatio, Rng& rng) {
    const double p = exp (energyDelta / params.temp) / fwdBackRatio;
    return p >= 1 || rng.real() < p;
  }
  inline void moveUnit (Unit& u, const Vec& pos, bool rev) {
    cell (u.pos, u.rev) = -1;
    u.pos = pos;
    u.rev = rev;
    cell (u.pos, u.rev) = u.index;
  }

  bool tryMove (Rng&);
  void dump (ostream&) const;

  // Monte Carlo sweeps: attempts per mobile unit.  This is the simulator's
  // time axis, and the only quantity an experiment can be mapped onto.
  inline double sweeps() const {
    return mobile.size() ? (attempts / (double) mobile.size()) : 0;
  }
  inline double seconds() const { return sweeps() * params.secondsPerSweep; }

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
  bool isLinear() const;
  void assertLinear() const;
  vguard<IndexPair> indexPairs() const;
  string sequence() const;
  vguard<Vec> unitPos() const;
  vguard<double> unitCentroid() const;
  double unitRadiusOfGyration() const;
  string foldString() const;
  string foldString (const vguard<IndexPair>&) const;
  double foldEnergy() const;
  double foldEnergy (const vguard<IndexPair>&) const;
  string coloredFoldString() const;
  string coloredFoldString (const vguard<IndexPair>&) const;

  // multi-chain logging
  map<string,int> sequenceFreqs() const;
};

#endif /* CELL_INCLUDED */
