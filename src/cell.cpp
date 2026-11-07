#include <set>
#include <algorithm>
#include "cell.h"

Params Params::fromJson (json& j) {
  Params p;
  p.splitProb = j["split"];
  p.stackEnergy = j["stack"];
  p.auEnergy = j["au"];
  p.gcEnergy = j["gc"];
  p.guEnergy = j["gu"];
  p.temp = j["temp"];
  if (j.count("bond")) p.bondProb = j["bond"];
  if (j.count("secondsPerSweep")) p.secondsPerSweep = j["secondsPerSweep"];
  if (j.count("wobble")) p.allowWobble = j["wobble"];
  return p;
}

json Params::toJson() const {
  json j;
  j["split"] = splitProb;
  j["stack"] = stackEnergy;
  j["au"] = auEnergy;
  j["gc"] = gcEnergy;
  j["gu"] = guEnergy;
  j["temp"] = temp;
  j["bond"] = bondProb;
  j["secondsPerSweep"] = secondsPerSweep;
  j["wobble"] = allowWobble;
  return j;
}

string Unit::alphabet ("acgu");
bool Unit::dnaMode = false;

static inline int popcount (uint32_t m) { return __builtin_popcount (m); }
static inline int nthSetBit (uint32_t m, int n) {
  while (n--) m &= m - 1;
  return __builtin_ctz (m);
}

Board::Board()
  : allMask(0), dist(0,1), baseDist(0,3),
    xSize(1), ySize(1), zSize(1), xWrap(0), yWrap(0), zWrap(0),
    legalMovesOnly(false), attempts(0)
{ }

Board::Board (int xs, int ys, int zs)
  : cellStorage (2*xs*ys*zs, -1), allMask(0), dist(0,1), baseDist(0,3),
    xSize(xs), ySize(ys), zSize(zs),
    xWrap (wrapOf (xs)), yWrap (wrapOf (ys)), zWrap (wrapOf (zs)),
    legalMovesOnly(false), attempts(0)
{
  for (int x = -nbrRange(xs); x <= nbrRange(xs); ++x)
    for (int y = -nbrRange(ys); y <= nbrRange(ys); ++y)
      for (int z = -nbrRange(zs); z <= nbrRange(zs); ++z)
	if (x != 0 || y != 0 || z != 0)
	  neighborhood.push_back (Vec (x, y, z));
  if (neighborhood.size() > 32)
    throw runtime_error ("Neighbourhood too large for a 32-bit move mask");
  allMask = (neighborhood.size() >= 32) ? 0xffffffffu : ((1u << neighborhood.size()) - 1);
  buildAdjMask();
}

// adjMask[c] is the set of neighbourhood directions d such that a unit stepping
// along d stays adjacent to a chain neighbour sitting at offset c.  Adjacency is
// translation invariant under the periodic wrap, so this can be tabulated once
// from relative coordinates and then reused for every unit on every move.
void Board::buildAdjMask() {
  adjMask.assign (27, 0);
  for (int dx = -1; dx <= 1; ++dx)
    for (int dy = -1; dy <= 1; ++dy)
      for (int dz = -1; dz <= 1; ++dz) {
	const Vec offset (dx, dy, dz);
	uint32_t m = 0;
	for (size_t d = 0; d < neighborhood.size(); ++d)
	  if (adjacent (offset, neighborhood[d]))
	    m |= (1u << d);
	adjMask[offsetCode (dx, dy, dz)] = m;
      }
}

void Board::rebuildMobile() {
  mobile.clear();
  for (const auto& u: unit)
    if (!u.pinned)
      mobile.push_back (u.index);
}

void Board::rebuildMasks() {
  for (auto& u: unit)
    u.mask = computeMask (u);
}

Board Board::fromJson (json& j) {
  json& js = j["size"];
  Board board (js[0], js[1], js[2]);
  board.params = Params::fromJson (j["params"]);
  for (auto& ju : j["unit"]) {
    const int index = board.unit.size();
    auto& jp = ju["pos"];
    Unit unit (Unit::char2base (ju["base"].get<string>()[0]),
	       jp[0].get<int>(),
	       jp[1].get<int>(),
	       jp[2].get<int>(),
	       ju.count("rev") && ju["rev"].get<bool>(),
	       index,
	       ju.count("prev") ? ju["prev"].get<int>() : -1,
	       -1);
    unit.pinned = ju.count("pinned") && ju["pinned"].get<bool>();
    unit.inert = ju.count("inert") && ju["inert"].get<bool>();
    board.unit.push_back (unit);
    board.cell (unit.pos, unit.rev) = index;
  }
  for (Unit& u: board.unit)
    if (u.prev >= 0)
      board.unit[u.prev].next = u.index;
  board.rebuildMobile();
  board.rebuildMasks();
  board.assertValid();
  return board;
}

json Board::toJson() const {
  assertValid();
  json j;
  j["size"] = { xSize, ySize, zSize };
  j["params"] = params.toJson();
  j["sweeps"] = sweeps();
  if (unit.size()) {
    json units;
    for (auto& u: unit) {
      json ju = {{ "base", string (1, Unit::base2char (u.base)) },
		 { "pos", { u.pos.x(), u.pos.y(), u.pos.z() } }};
      if (u.rev) ju["rev"] = true;
      if (u.prev >= 0) ju["prev"] = u.prev;
      if (u.next >= 0) ju["next"] = u.next;
      if (u.pinned) ju["pinned"] = true;
      if (u.inert) ju["inert"] = true;
      units.push_back (ju);
    }
    j["unit"] = units;
  }
  return j;
}

void Board::addSeq (const string& seq) {
  if (xSize < (int) seq.length())
    throw runtime_error ("Board is too small for sequence");
  for (size_t pos = 0; pos < seq.length(); ++pos) {
    if (cell (pos, 0, 0, false) != -1)
      throw runtime_error ("Cell occupied");
    const char c = tolower (seq.at(pos));
    if (!Unit::isRNA (c))
      throw runtime_error ("Sequence is not a nucleic acid");
    const int index = unit.size();
    Unit u (Unit::char2base (c), pos, 0, 0, false, index, pos ? (index - 1) : -1, -1);
    if (pos)
      unit.back().next = index;
    unit.push_back (u);
    cell (u.pos, false) = index;
  }
  rebuildMobile();
  rebuildMasks();
}

void Board::addBases (double density, Rng& rng) {
  for (int x = 0; x < xSize; ++x)
    for (int y = 0; y < ySize; ++y)
      for (int z = 0; z < zSize; ++z)
	if (cell(x,y,z,false) < 0
	    && cell(x,y,z,true) < 0
	    && rng.real() < density) {
	  const int index = unit.size();
	  Unit u (baseDist(rng), x, y, z, false, index, -1, -1);
	  unit.push_back (u);
	  cell (u.pos, false) = index;
	}
  rebuildMobile();
  rebuildMasks();
}

// An inert surface: every site in the plane at coordinate 0 normal to `axis` is
// filled with a unit that never diffuses and never base-pairs.  It is simply an
// excluded-volume wall that strands can also be tethered to.
void Board::addSurface (int axis) {
  if (axis < 0 || axis > 2)
    throw runtime_error ("Surface axis must be x, y or z");
  const int size[3] = { xSize, ySize, zSize };
  const int a1 = (axis + 1) % 3, a2 = (axis + 2) % 3;
  for (int i = 0; i < size[a1]; ++i)
    for (int k = 0; k < size[a2]; ++k) {
      int c[3] = { 0, 0, 0 };
      c[axis] = 0; c[a1] = i; c[a2] = k;
      if (cell (c[0], c[1], c[2], false) >= 0)
	continue;
      const int index = unit.size();
      Unit u (0, c[0], c[1], c[2], false, index, -1, -1);
      u.pinned = true;
      u.inert = true;
      unit.push_back (u);
      cell (u.pos, false) = index;
    }
  rebuildMobile();
  rebuildMasks();
}

// A 2D array of tethered strands: copies of `seq` on a grid of spacing `spacing`
// in the plane normal to `axis`, each extending away from the surface and
// covalently linked at its 5' end to an inert anchor site in the plane.
int Board::addTetheredArray (const string& seq, int spacing, int axis) {
  if (axis < 0 || axis > 2)
    throw runtime_error ("Tether axis must be x, y or z");
  if (spacing < 1)
    throw runtime_error ("Tether spacing must be at least 1");
  const int size[3] = { xSize, ySize, zSize };
  if (size[axis] < (int) seq.length() + 1)
    throw runtime_error ("Board is too small along the tether axis for this sequence");
  const int a1 = (axis + 1) % 3, a2 = (axis + 2) % 3;
  string lower;
  for (char c: seq) {
    const char l = tolower (c);
    if (!Unit::isRNA (l))
      throw runtime_error ("Sequence is not a nucleic acid");
    lower.push_back (l);
  }
  int strands = 0;
  for (int i = 0; i < size[a1]; i += spacing)
    for (int k = 0; k < size[a2]; k += spacing) {
      int c[3] = { 0, 0, 0 };
      c[a1] = i; c[a2] = k;
      // anchor at coordinate 0 along the axis; reuse a surface site if present
      c[axis] = 0;
      int anchor = cell (c[0], c[1], c[2], false);
      if (anchor < 0) {
	anchor = unit.size();
	Unit a (0, c[0], c[1], c[2], false, anchor, -1, -1);
	a.pinned = true;
	a.inert = true;
	unit.push_back (a);
	cell (a.pos, false) = anchor;
      } else if (!unit[anchor].inert || unit[anchor].next >= 0)
	continue;  // occupied by something that cannot serve as an anchor
      bool blocked = false;
      for (size_t n = 0; n < lower.size() && !blocked; ++n) {
	c[axis] = n + 1;
	if (cell (c[0], c[1], c[2], false) >= 0 || cell (c[0], c[1], c[2], true) >= 0)
	  blocked = true;
      }
      if (blocked)
	continue;
      int prev = anchor;
      for (size_t n = 0; n < lower.size(); ++n) {
	c[axis] = n + 1;
	const int index = unit.size();
	Unit u (Unit::char2base (lower[n]), c[0], c[1], c[2], false, index, prev, -1);
	unit.push_back (u);
	unit[prev].next = index;
	cell (u.pos, false) = index;
	prev = index;
      }
      ++strands;
    }
  rebuildMobile();
  rebuildMasks();
  return strands;
}

// Renumbering units in lattice order is a pure relabelling: with uniform unit
// selection the Markov chain is unchanged, but the random walk over the unit
// array touches far fewer cache lines once the system is large.
void Board::sortUnitsByPosition() {
  const size_t n = unit.size();
  vguard<int> order (n);
  for (size_t i = 0; i < n; ++i) order[i] = i;
  stable_sort (order.begin(), order.end(), [&](int a, int b) {
      return cellIndex (unit[a].pos) < cellIndex (unit[b].pos);
    });
  vguard<int> newIndex (n);
  for (size_t k = 0; k < n; ++k) newIndex[order[k]] = k;
  vguard<Unit> relabelled (n);
  for (size_t k = 0; k < n; ++k) {
    Unit u = unit[order[k]];
    u.index = k;
    if (u.prev >= 0) u.prev = newIndex[u.prev];
    if (u.next >= 0) u.next = newIndex[u.next];
    relabelled[k] = u;
  }
  unit = relabelled;
  fill (cellStorage.begin(), cellStorage.end(), -1);
  for (const auto& u: unit)
    cell (u.pos, u.rev) = u.index;
  rebuildMobile();
  rebuildMasks();
}

void Board::assertValid() const {
  set<int> seen;
  for (int x = 0; x < xSize; ++x)
    for (int y = 0; y < ySize; ++y)
      for (int z = 0; z < zSize; ++z)
	for (int rev = 0; rev <= 1; rev++) {
	  const int idx = cell (x, y, z, rev);
	  if (idx >= 0) {
	    if (seen.count(idx))
	      throw runtime_error ("Duplicate Unit index");
	    const Unit u = unit[idx];
	    if (u.index != idx)
	      throw runtime_error ("Incorrect Unit index");
	    if (!boardCoordsEqual (Vec(x,y,z), u.pos) || (rev ? !u.rev : u.rev))
	      throw runtime_error ("Mislocated Unit");
	    if (u.prev >= 0 && unit[u.prev].next != idx)
	      throw runtime_error ("Broken Unit.prev");
	    if (u.next >= 0 && unit[u.next].prev != idx)
	      throw runtime_error ("Broken Unit.next");
	    if (rev && cell (x, y, z, false) < 0)
	      throw runtime_error ("Unit.rev is true but no paired Unit exists");
	    seen.insert (idx);
	  }
	}
  if (seen.size() != unit.size())
    throw runtime_error ("Missing Unit");
}

bool Board::tryMove (Rng& rng) {
  if (mobile.empty()) {
    attempts += 1;
    return false;
  }
  bool moved = false;
  const int index = mobile [rng.index (mobile.size())];
  Unit& u = unit[index];
  const int nDir = neighborhood.size();
  const uint32_t uMask = u.mask;

  int dir;
  double propRatio = 1;   // q(forward)/q(reverse) contributed by the proposal restriction
  int mFwd = nDir;
  if (legalMovesOnly) {
    mFwd = popcount (uMask);
    if (!mFwd) {
      attempts += 1;
      return false;
    }
    dir = nthSetBit (uMask, rng.index (mFwd));
    attempts += (double) nDir / (double) mFwd;   // this proposal stands in for D/m uniform attempts
  } else {
    dir = rng.index (nDir);
    attempts += 1;
    if (!((uMask >> dir) & 1u))   // one bit test in place of the old canMoveTo call
      return false;
  }

  const Vec& delta = neighborhood[dir];
  const Vec newPos = u.pos + delta;
  const int ci = cellIndex (newPos);
  const int nbrIndex = cellStorage[ci];
  const int nbrPairIndex = cellStorage[ci + 1];

  // ratio of reverse to forward proposal weight for a move that relocates u
  if (legalMovesOnly)
    propRatio = (double) popcount (maskAt (u, newPos)) / (double) mFwd;

  if (isPaired (u)) {
    Unit& p = unit[pairedIndex(u)];
    const double oldEnergy = pairingEnergy (u, p);
    if (rng.real() < params.splitProb) {
      // attempt split
      if (nbrIndex < 0) {
	if (acceptMove (-oldEnergy, params.splitProb * propRatio, rng)) {
	  moveUnit (u, newPos, false);
	  moveUnit (p, p.pos, false);
	  refreshAround (u);
	  moved = true;
	}
      } else {
	Unit& nbr = unit[nbrIndex];
	if (nbrPairIndex < 0 && canMerge (u, nbr)) {
	  if (acceptMove (pairingEnergy(u,nbr) - oldEnergy, propRatio, rng)) {
	    moveUnit (u, newPos, true);
	    moveUnit (p, p.pos, false);
	    refreshAround (u);
	    moved = true;
	  }
	}
      }
    } else {  // paired and not attempting split
      if (nbrIndex < 0 && nbrPairIndex < 0 && !p.pinned && ((p.mask >> dir) & 1u)) {
	// translate the whole pair.  Either member could have been picked to
	// propose this, so the proposal weight is the sum of both.
	double pairRatio = 1;
	if (legalMovesOnly) {
	  const double fwd = 1.0 / mFwd + 1.0 / popcount (p.mask);
	  const double rev = 1.0 / popcount (maskAt (u, newPos)) + 1.0 / popcount (maskAt (p, newPos));
	  pairRatio = fwd / rev;
	}
	if (!legalMovesOnly || acceptMove (0, pairRatio, rng)) {
	  moveUnit (u, newPos, u.rev);
	  moveUnit (p, newPos, p.rev);
	  refreshAround (u);
	  refreshAround (p);
	  moved = true;
	}
      } else if (nbrIndex >= 0 && nbrPairIndex >= 0 && u.next < 0 && params.bondProb > 0) {
	// Template-directed ligation.  This is a reaction, not a diffusion move,
	// so it is gated by its own rate rather than by the energy model.  Under
	// restricted proposals the attempt rate is inflated by D/m, so the bond
	// probability is scaled back by m/D to leave the physical rate alone.
	const double pBond = legalMovesOnly
	  ? params.bondProb * (double) mFwd / (double) nDir
	  : params.bondProb;
	if (rng.real() < pBond) {
	  Unit& nbr = unit[nbrIndex];
	  Unit& nbrp = unit[nbrPairIndex];
	  if (p.prev == nbrIndex && nbrp.prev < 0 && !nbrp.inert && !u.inert) {
	    nbrp.prev = index;
	    u.next = nbrPairIndex;
	    refreshMask (index);
	    refreshMask (nbrPairIndex);
	    moved = true;
	  } else if (p.prev == nbrPairIndex && nbr.prev < 0 && !nbr.inert && !u.inert) {
	    nbr.prev = index;
	    u.next = nbrIndex;
	    refreshMask (index);
	    refreshMask (nbrIndex);
	    moved = true;
	  }
	}
      }
    }
  } else {  // not paired
    if (nbrIndex < 0) {
      if (!legalMovesOnly || acceptMove (0, propRatio, rng)) {
	moveUnit (u, newPos, false);
	refreshAround (u);
	moved = true;
      }
    } else {
      Unit& nbr = unit[nbrIndex];
      if (nbrPairIndex < 0 && canMerge (u, nbr)) {
	if (acceptMove (pairingEnergy(u,nbr), propRatio / params.splitProb, rng)) {
	  moveUnit (u, newPos, true);
	  refreshAround (u);
	  moved = true;
	}
      }
    }
  }
  return moved;
}

void Board::dump (ostream& out) const {
  for (int x = 0; x < xSize; ++x)
    for (int y = 0; y < ySize; ++y)
      for (int z = 0; z < zSize; ++z)
	for (int rev = 0; rev <= 1; rev++) {
	  const int idx = cell (x, y, z, rev);
	  out << "(" << x << "," << y << "," << z << ")." << (rev ? "1" : "0") << " #" << idx;
	  if (idx >= 0) {
	    const Unit& u = unit[idx];
	    out << ": " << u.base << " " << u.pos << "." << u.rev << " #" << u.index << " prev=#" << u.prev << " next=#" << u.next;
	  }
	  out << endl;
	}
}

vguard<Board::IndexPair> Board::indexPairs() const {
  vguard<IndexPair> p;
  for (size_t i = 0; i < unit.size(); ++i) {
    const Unit& u = unit[i];
    const int j = pairedIndex(u);
    if (j > (int) i)
      p.push_back (IndexPair (i, j));
  }
  return p;
}

bool Board::isLinear() const {
  for (size_t i = 0; i < unit.size(); ++i) {
    const Unit& u = unit[i];
    if (u.index != (int) i || (i > 0 && u.prev != (int) i-1) || (i < unit.size()-1 && u.next != (int) i+1))
      return false;
  }
  return true;
}

void Board::assertLinear() const {
  if (!isLinear())
    throw runtime_error ("Board does not contain a single linear chain");
}

string Board::sequence() const {
  string s;
  s.reserve (unit.size());
  for (auto& u: unit)
    s.push_back (u.baseChar());
  return s;
}

string Board::leftFoldChar ("<[{(abcdefghijklmnopqrstuvwxyz");
string Board::rightFoldChar (">]})ABCDEFGHIJKLMNOPQRSTUVWXYZ");

// Pairs arrive sorted by their left index, so the intervals already placed on a
// given level are nested around the current position once the closed ones are
// popped.  Only the innermost needs testing, which makes this linear in the
// number of pairs instead of quadratic.
string Board::foldString (const vguard<IndexPair>& pairs) const {
  string fs (unit.size(), '.');
  const size_t nLevels = leftFoldChar.size();
  vguard<vguard<IndexPair> > open (nLevels);
  bool warned = false;
  for (const auto& ij: pairs) {
    const int i = ij.first, j = ij.second;
    size_t lvl;
    for (lvl = 0; lvl < nLevels; ++lvl) {
      vguard<IndexPair>& st = open[lvl];
      while (!st.empty() && st.back().second < i)
	st.pop_back();
      if (st.empty() || st.back().second > j)
	break;
    }
    if (lvl < nLevels) {
      open[lvl].push_back (ij);
      fs[i] = leftFoldChar[lvl];
      fs[j] = rightFoldChar[lvl];
    } else if (!warned) {
      cerr << "Not enough fold characters!" << endl;
      warned = true;
    }
  }
  return fs;
}
string Board::foldString() const { return foldString (indexPairs()); }

string Board::coloredFoldString (const vguard<IndexPair>& pairs) const {
  const string fs = foldString (pairs);
  vguard<int> col (fs.size(), 7);
  int c = 1, last_i = -1, last_j = -1;
  for (auto& ij: pairs) {
    if (last_i >= 0 && (ij.first != last_i + 1 || ij.second != last_j - 1))
      c = (c % 6) + 1;
    col[ij.first] = col[ij.second] = c;
    last_i = ij.first;
    last_j = ij.second;
  }
  string cfs;
  cfs.reserve ((fs.size() + 1) * 6);
  for (size_t pos = 0; pos < fs.size(); ++pos)
    cfs.append (string("\033[") + to_string(30+col[pos]) + "m" + fs[pos]);
  cfs.append ("\033[37m");
  return cfs;
}
string Board::coloredFoldString() const { return coloredFoldString (indexPairs()); }

double Board::foldEnergy (const vguard<IndexPair>& pairs) const {
  double e = 0;
  for (const auto& ij: pairs)
    e += calcEnergy (unit[ij.first], unit[ij.second], 0.5);
  return e;
}
double Board::foldEnergy() const { return foldEnergy (indexPairs()); }

vguard<Vec> Board::unitPos() const {
  vguard<Vec> up;
  up.reserve (unit.size());
  for (auto& u: unit)
    up.push_back (u.pos);
  return up;
}

vguard<double> Board::unitCentroid() const {
  vguard<double> c (3);
  for (auto& u: unit)
    for (size_t n = 0; n < 3; ++n)
      c[n] += u.pos.xyz[n];
  for (size_t n = 0; n < 3; ++n)
    c[n] /= unit.size();
  return c;
}

double Board::unitRadiusOfGyration() const {
  const vguard<double> c = unitCentroid();
  double d2 = 0;
  for (auto& u: unit)
    for (size_t n = 0; n < 3; ++n) {
      const double d = u.pos.xyz[n] - c[n];
      d2 += d*d;
    }
  return sqrt (d2 / unit.size());
}

map<string,int> Board::sequenceFreqs() const {
  map<string,int> seqFreq;
  vguard<bool> unitSeen (unit.size());
  size_t nSeen = 0;
  for (size_t i = 0; i < unit.size(); ++i)
    if (!unitSeen[i]) {
      int j = i;
      bool cyclic = false;
      while (unit[j].prev >= 0) {
	j = unit[j].prev;
	if (j == (int) i) {
	  cyclic = true;
	  break;
	}
      }
      string s;
      while (j >= 0 && !unitSeen[j]) {
	unitSeen[j] = true;
	s.push_back (unit[j].baseChar());
	j = unit[j].next;
	++nSeen;
      }
      if (cyclic) {
	string canonical = s;
	for (size_t k = 1; k < s.size(); ++k) {
	  const string perm = s.substr(k,s.size()-k) + s.substr(0,k);
	  if (perm < canonical)
	    canonical = perm;
	}
	s = canonical + "*";
      }
      ++seqFreq[s];
    }
  if (nSeen != unit.size())
    throw runtime_error ("Missed Units");
  return seqFreq;
}
