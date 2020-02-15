#include <set>
#include "cell.h"

Params fromJson (json& j) {
  Params p;
  p.splitProb = j["split"];
  return p;
}

json Params::toJson() const {
  json j;
  j["split"] = splitProb;
  return j;
}

Board::Board()
{ }

Board::Board (int xs, int ys, int zs)
  : xSize(xs), ySize(ys), zSize(zs), cellStorage (2*xs*ys*zs, -1)
{
  for (int x = -1; x <= 1; ++x)
    for (int y = -1; y <= 1; ++y)
      for (int z = -1; z <= 1; ++z)
	if (x != 0 || y != 0 || z != 0)
	  neighborhood.push_back (Vec (x, y, z));
}

Board Board::fromJson (json& j) {
  json& js = j["size"];
  Board board (js[0], js[1], js[2]);
  for (auto& ju : j["unit"]) {
    auto& jp = ju["pos"];
    Unit unit (ju["base"].get<string>()[0],
	       jp[0].get<int>(),
	       jp[1].get<int>(),
	       jp[2].get<int>(),
	       jp.count("rev") && jp["rev"].get<bool>(),
	       ju.count("prev") ? ju["prev"].get<int>() : -1,
	       ju.count("next") ? ju["next"].get<int>() : -1);
    const int index = board.unit.size();
    board.unit.push_back (unit);
    board.cell (unit.pos, unit.isRev) = index;
  }
  return board;
}

json Board::toJson() const {
  assertValid();
  json j;
  j["size"] = { xSize, ySize, zSize };
  if (unit.size()) {
    json units;
    for (auto& u: unit) {
      json ju = {{ "base", string (1, u.base) },
		 { "pos", { u.pos.x(), u.pos.y(), u.pos.z() } }};
      if (u.isRev) ju["rev"] = true;
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
	for (bool rev = false; !rev; rev = true) {
	  const int idx = cell (x, y, z, rev);
	  if (idx >= 0) {
	    if (seen.count(idx))
	      throw runtime_error ("Duplicate Unit index");
	    const Unit u = unit[idx];
	    if (x != u.pos.x() || y != u.pos.y() || z != u.pos.z())
	      throw runtime_error ("Mislocated Unit");
	    if (u.prev >= 0 && unit[u.prev].next != idx)
	      throw runtime_error ("Broken Unit.prev");
	    if (u.next >= 0 && unit[u.next].prev != idx)
	      throw runtime_error ("Broken Unit.next");
	    seen.insert (idx);
	  }
	}
  if (seen.size() != unit.size())
    throw runtime_error ("Missing Unit");
}

bool Board::tryMove (int index) {
  
  return false;
}
