#include <set>
#include "cell.h"

Params Params::fromJson (json& j) {
  Params p;
  p.splitProb = j["split"];
  p.mismatchProb = j["mismatch"];
  return p;
}

json Params::toJson() const {
  json j;
  j["split"] = splitProb;
  j["mismatch"] = mismatchProb;
  return j;
}

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
    Unit unit (ju["base"].get<string>()[0],
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
      json ju = {{ "base", string (1, u.base) },
		 { "pos", { u.pos.x(), u.pos.y(), u.pos.z() } }};
      if (u.rev) ju["rev"] = true;
      if (u.prev >= 0) ju["prev"] = u.prev;
      if (u.next >= 0) ju["next"] = u.next;
      units.push_back (ju);
    }
    j["unit"] = units;
  }
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
    Unit u (c,
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
	    if (x != u.pos.x() || y != u.pos.y() || z != u.pos.z() || (rev ? !u.rev : u.rev))
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
	Unit& p = unit[cell(u.pos,!u.rev)];
	//	cerr << "Paired unit is at " << p.pos << endl;
	if (dist(mt) < params.splitProb) {
	  // attempt split
	  //	  cerr << "Attempting split" << endl;
	  if (nbrIndex < 0) {
	    // split and move to forward slot
	    moveUnit (u, newPos, false);
	    moveUnit (p, p.pos, false);
	    //	    cerr << "Paired unit is now at " << p.pos << "." << p.rev << endl;
	    moved = true;
	  } else {
	    Unit& nbr = unit[nbrIndex];
	    if (nbrPairIndex < 0 && acceptMerge (u, nbr, mt)) {
	      // split and move to rev slot
	      moveUnit (u, newPos, true);
	      moveUnit (p, p.pos, false);
	      //	      cerr << "Paired unit is now at " << p.pos << "." << p.rev << endl;
	      moved = true;
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
	  if (nbrPairIndex < 0 && acceptMerge (u, nbr, mt)) {
	    // move to rev slot
	    moveUnit (u, newPos, true);
	    moved = true;
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
