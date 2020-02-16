#include <set>
#include "cell.h"

Params Params::fromJson (json& j) {
  Params p;
  p.splitProb = j["split"];
  p.stackEnergy = j["stack"];
  p.auEnergy = j["au"];
  p.gcEnergy = j["gc"];
  p.guEnergy = j["gu"];
  p.temp = j["temp"];
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
  return j;
}

string Unit::alphabet ("acgu");

Board::Board() : dist(0,1)
{ }

Board::Board (int xs, int ys, int zs)
  : xSize(xs), ySize(ys), zSize(zs), cellStorage (2*xs*ys*zs, -1), dist(0,1)
{
  for (int x = -nbrRange(xs); x <= nbrRange(xs); ++x)
    for (int y = -nbrRange(ys); y <= nbrRange(ys); ++y)
      for (int z = -nbrRange(zs); z <= nbrRange(zs); ++z)
	if (x != 0 || y != 0 || z != 0)
	  neighborhood.push_back (Vec (x, y, z));
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
	       jp.count("rev") && jp["rev"].get<bool>(),
	       index,
	       ju.count("prev") ? ju["prev"].get<int>() : -1,
	       ju.count("next") ? ju["next"].get<int>() : -1);
    board.unit.push_back (unit);
    board.cell (unit.pos, unit.rev) = index;
  }
  return board;
}

json Board::toJson() const {
  assertValid();
  json j;
  j["size"] = { xSize, ySize, zSize };
  j["params"] = params.toJson();
  if (unit.size()) {
    json units;
    for (auto& u: unit) {
      json ju = {{ "base", string (1, u.baseChar()) },
		 { "pos", { u.pos.x(), u.pos.y(), u.pos.z() } }};
      if (u.rev) ju["rev"] = true;
      if (u.prev >= 0) ju["prev"] = u.prev;
      if (u.next >= 0) ju["next"] = u.next;
      units.push_back (ju);
    }
    j["unit"] = units;
  }
  j["fold"] = foldString();
  j["energy"] = foldEnergy();
  j["sequence"] = sequence();
  return j;
}

void Board::addSeq (const string& seq) {
  if (xSize < seq.length())
    throw runtime_error ("Board is too small for sequence");
  for (size_t pos = 0; pos < seq.length(); ++pos) {
    if (cell (pos, 0, 0, false) != -1)
      throw runtime_error ("Cell occupied");
    const char c = tolower (seq.at(pos));
    if (!Unit::isRNA (c))
      throw runtime_error ("Sequence is not RNA");
    const int index = unit.size();
    Unit u (Unit::char2base (c),
	    pos,
	    0,
	    0,
	    false,
	    index,
	    pos ? (index - 1) : -1,
	    -1);
    if (pos)
      unit.back().next = index;
    unit.push_back (u);
    cell (u.pos, false) = index;
  }
}

const Vec& Board::rndNbrVec (mt19937& mt) const {
  return neighborhood [mt() % neighborhood.size()];
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

bool Board::tryMove (mt19937& mt) {
  bool moved = false;
  if (unit.size()) {
    const int index = mt() % unit.size();
    Unit& u = unit[index];
    const Vec& delta = rndNbrVec (mt);
    const Vec newPos = u.pos + delta;
    //    cerr << "Attempting to move unit #" << index << " from " << u.pos << " to " << newPos << endl;
    if (canMoveTo (u, newPos)) {
      const int nbrIndex = cell (newPos, false);
      const int nbrPairIndex = cell (newPos, true);
      if (isPaired (u)) {
	Unit& p = unit[pairedIndex(u)];
	const double oldEnergy = pairingEnergy (u, p);
	//	cerr << "Paired unit is at " << p.pos << endl;
	if (dist(mt) < params.splitProb) {
	  // attempt split
	  //	  cerr << "Attempting split" << endl;
	  if (nbrIndex < 0) {
	    if (acceptMove (-oldEnergy, params.splitProb, mt)) {
	      // split and move to forward slot
	      moveUnit (u, newPos, false);
	      moveUnit (p, p.pos, false);
	      //	    cerr << "Paired unit is now at " << p.pos << "." << p.rev << endl;
	      moved = true;
	    }
	  } else {
	    Unit& nbr = unit[nbrIndex];
	    if (nbrPairIndex < 0 && canMerge (u, nbr)) {
	      if (acceptMove (pairingEnergy(u,nbr) - oldEnergy, 1, mt)) {
		// split and move to rev slot
		moveUnit (u, newPos, true);
		moveUnit (p, p.pos, false);
		//	      cerr << "Paired unit is now at " << p.pos << "." << p.rev << endl;
		moved = true;
	      }
	    }
	  }
	} else {  // paired and not attempting split
	  if (nbrIndex < 0 && nbrPairIndex < 0 && canMoveTo (p, newPos)) {
	    // move both u & p
	    moveUnit (u, newPos, u.rev);
	    moveUnit (p, newPos, p.rev);
	    //	    cerr << "Paired unit is now at " << p.pos << "." << p.rev << endl;
	    moved = true;
	  }
	}
      } else {  // not paired
	if (nbrIndex < 0) {
	  // move to forward slot
	  moveUnit (u, newPos, false);
	  moved = true;
	} else {
	  Unit& nbr = unit[nbrIndex];
	  if (nbrPairIndex < 0 && canMerge (u, nbr)) {
	    if (acceptMove (pairingEnergy(u,nbr), 1. / params.splitProb, mt)) {
	      // move to rev slot
	      moveUnit (u, newPos, true);
	      moved = true;
	    }
	  }
	}
      }
    }
  }
  //  assertValid();
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
  for (int i = 0; i < unit.size(); ++i) {
    const Unit& u = unit[i];
    const int j = pairedIndex(u);
    if (j > i)
      p.push_back (IndexPair (i, j));
  }
  return p;
}

void Board::assertLinear() const {
  for (size_t i = 0; i < unit.size(); ++i) {
    const Unit& u = unit[i];
    if (u.index != i || (i > 0 && u.prev != i-1) || (i < unit.size()-1 && u.next != i+1))
      throw runtime_error ("foldString requires single linear chain");
  }
}

string Board::sequence() const {
  string s;
  s.reserve (unit.size());
  for (auto& u: unit)
    s.push_back (u.baseChar());
  return s;
}

string Board::foldString() const {
  string fs (unit.size(), '.');
  const string leftChar  = "<[{(0123456789abcdefghijklmnopqrstuvwxyz";
  const string rightChar = ">]})0123456789abcdefghijklmnopqrstuvwxyz";
  map<size_t,set<IndexPair>> offsetPairs;
  for (const auto& ij: indexPairs()) {
    const int i = ij.first, j = ij.second;
    //      cerr << ij.first << " <--> " << ij.second << endl;
    size_t offset;
    for (offset = 0; offset < leftChar.size(); ++offset) {
      bool intersects = false;
      for (const auto& op: offsetPairs[offset]) {
	const int a = op.first, b = op.second;
	if ((i < a && a < j && j < b)
	    || (a < i && i < b && b < j)) {
	  intersects = true;
	  break;
	}
      }
      if (!intersects)
	break;
    }
    if (offset < leftChar.size()) {
      offsetPairs[offset].insert (ij);
      fs[ij.first] = leftChar[offset];
      fs[ij.second] = rightChar[offset];
    } else
      cerr << "Not enough fold characters!" << endl;
  }
  return fs;
}

double Board::foldEnergy() const {
  double e = 0;
  for (const auto& ij: indexPairs())
    e += calcEnergy (unit[ij.first], unit[ij.second], 0.5);
  return e;
}

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
