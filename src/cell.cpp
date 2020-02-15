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
  json j;
  j["size"] = { xSize, ySize, zSize };
  json units;
  for (auto& u: unit) {
    json ju = {{ "base", string (1, u.base),
		 "pos", { u.pos.x(), u.pos.y(), u.pos.z() } }};
    if (u.isRev) ju["rev"] = true;
    if (u.prev >= 0) ju["prev"] = u.prev;
    if (u.next >= 0) ju["next"] = u.next;
    units.push_back (ju);
  }
  return j;
}

const Vec& Board::rndNbrVec (mt19937& mt) const {
  return neighborhood [mt() % neighborhood.size()];
}

bool Board::cellIsValid (int, int, int) const {
  return false;
}

bool Board::tryMove (const Vec&, const Vec&) {
  return false;
}

bool Board::trySplit (const Vec&, const Vec&) {
  return false;
}
