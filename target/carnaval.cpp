#include <time.h>
#include <cstdlib>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <unordered_map>
#include <boost/program_options.hpp>

#include "../src/cell.h"
#include "../src/rng.h"
#include "../src/util.h"
#include "../src/bitmap_image.hpp"

using namespace std;
namespace po = boost::program_options;

// Base-pair occupancy counters.  A std::map keyed on an index pair costs a
// pointer chase and a comparison per pair per sample; for the single-chain runs
// that dominate analysis a dense triangle is both smaller and far faster, and
// the hash map keeps large multi-strand boards from degrading.
struct PairCounter {
  size_t n;
  bool dense;
  vguard<uint32_t> grid;
  unordered_map<uint64_t,uint32_t> sparse;
  explicit PairCounter (size_t units) : n(units), dense (units && units <= 2048) {
    if (dense)
      grid.assign (n * n, 0);
  }
  inline void add (int i, int j) {
    if (dense)
      ++grid[(size_t) i * n + j];
    else
      ++sparse[((uint64_t) i << 32) | (uint32_t) j];
  }
  template<typename F> void forEach (F f) const {
    if (dense) {
      for (size_t i = 0; i < n; ++i)
	for (size_t j = 0; j < n; ++j)
	  if (grid[i*n+j])
	    f ((int) i, (int) j, grid[i*n+j]);
    } else
      for (const auto& kv: sparse)
	f ((int) (kv.first >> 32), (int) (kv.first & 0xffffffffu), kv.second);
  }
};

static int axisFromName (const string& s) {
  if (s == "x" || s == "X" || s == "0") return 0;
  if (s == "y" || s == "Y" || s == "1") return 1;
  if (s == "z" || s == "Z" || s == "2") return 2;
  throw runtime_error ("Axis must be x, y or z");
}

int main (int argc, char** argv) {

  try {

    po::options_description opts("Options");
    opts.add_options()
      ("help,h", "display this help message")
      ("xsize,x", po::value<int>()->default_value(64), "size of board in X dimension")
      ("ysize,y", po::value<int>()->default_value(64), "size of board in Y dimension")
      ("zsize,z", po::value<int>()->default_value(1), "size of board in Z dimension")
      ("au",  po::value<double>(), "specify energy of A-U bond pairs")
      ("gc",  po::value<double>(), "specify energy of G-C bond pairs")
      ("gu",  po::value<double>(), "specify energy of G-U bond pairs")
      ("stack",  po::value<double>(), "specify energy of base-pair stacking interactions")
      ("init,i",  po::value<string>(), "specify initial template sequence")
      ("density,d",  po::value<double>(), "specify initial density of monomers")
      ("bond,B",  po::value<double>(), "polymerization rate (bond-formation probability per attempt); 0 disables ligation")
      ("dna",  "nucleic acid is DNA: read and write T in place of U, and forbid wobble pairs")
      ("no-wobble",  "forbid G-U (G-T) wobble base pairs")
      ("surface",  po::value<string>(), "add an inert non-diffusing surface filling the plane at coordinate 0 normal to this axis (x, y or z)")
      ("tether",  po::value<int>(), "tether copies of the --init sequence to the surface on a grid of this spacing")
      ("tether-axis",  po::value<string>(), "axis along which tethered strands extend (default z, or y for a 2D board)")
      ("pin-first",  "pin the first unit of the --init sequence to its lattice site")
      ("sort-units",  "renumber units in lattice order for cache locality (relabelling only; changes unit indices)")
      ("legal-moves-only",  "propose only geometrically legal moves. Preserves the equilibrium ensemble but distorts relative kinetics, so use for structure prediction, not rates")
      ("rnd,r",  po::value<int>(), "seed random number generator")
      ("total-moves,t",  po::value<long>()->default_value(0), "total number of moves")
      ("unit-moves,u",  po::value<long>()->default_value(0), "number of moves per mobile unit")
      ("seconds-per-sweep",  po::value<double>(), "experimental seconds per Monte Carlo sweep, for calibrated output")
      ("folds,f",  "periodically log move count, fold string, energy, radius of gyration, and centroid (single-chain simulations only)")
      ("seqs,S",  "periodically log sequences (for replication simulations)")
      ("monochrome,m",  "no ANSI color codes in logging, please")
      ("period,p", po::value<long>()->default_value(1000), "logging period")
      ("temp,T",  po::value<double>(), "specify temperature")
      ("load,l", po::value<string>(), "load board state from file")
      ("save,s", po::value<string>(), "save board state to file")
      ("json,j", po::value<string>(), "save base-pairing posterior probabilities to JSON file")
      ("bitmap,b", po::value<string>(), "save base-pairing probabilities to bitmap image file")
      ("csv,c", po::value<string>(), "save base-pairing probabilities to CSV file")
      ;

    po::variables_map vm;
    po::parsed_options parsed = po::command_line_parser(argc,argv).options(opts).allow_unregistered().run();
    po::store (parsed, vm);
    po::notify(vm);

    if (vm.count("help")) {
      cout << opts << endl;
      return 1;
    }

    time_t timer;
    time (&timer);
    const int seed = vm.count("rnd") ? vm.at("rnd").as<int>() : (int) timer;
    Rng rng (seed);
    cerr << "Random seed is " << seed << endl;

    if (vm.count("dna"))
      Unit::dnaMode = true;

    // create Board
    Board board;
    if (vm.count("load")) {
      ifstream infile (vm.at("load").as<string>());
      if (!infile)
	throw runtime_error ("Can't load board file");
      json j;
      infile >> j;
      board = Board::fromJson (j);
    } else {
      board = Board (vm["xsize"].as<int>(),
		     vm["ysize"].as<int>(),
		     vm["zsize"].as<int>());
      if (board.xWrap < 0 || board.yWrap < 0 || board.zWrap < 0)
	cerr << "Note: board dimensions that are powers of two let the periodic wrap use a bit mask"
	     << " instead of integer division, which is appreciably faster." << endl;
    }

    // parameters
    if (vm.count("temp")) board.params.temp = vm.at("temp").as<double>();
    if (vm.count("bond")) board.params.bondProb = vm.at("bond").as<double>();
    if (vm.count("stack")) board.params.stackEnergy = vm.at("stack").as<double>();
    if (vm.count("au")) board.params.auEnergy = vm.at("au").as<double>();
    if (vm.count("gc")) board.params.gcEnergy = vm.at("gc").as<double>();
    if (vm.count("gu")) board.params.guEnergy = vm.at("gu").as<double>();
    if (vm.count("seconds-per-sweep")) board.params.secondsPerSweep = vm.at("seconds-per-sweep").as<double>();
    if (vm.count("dna") || vm.count("no-wobble")) board.params.allowWobble = false;
    board.legalMovesOnly = vm.count("legal-moves-only") > 0;

    // surface, before any tethering so that anchors reuse surface sites
    if (vm.count("surface")) {
      const int axis = axisFromName (vm.at("surface").as<string>());
      board.addSurface (axis);
      cerr << "Added inert surface normal to axis " << vm.at("surface").as<string>()
	   << " (" << board.unit.size() << " sites)" << endl;
    }

    // initialization
    if (vm.count("tether")) {
      if (!vm.count("init"))
	throw runtime_error ("--tether needs a sequence from --init");
      const int axis = vm.count("tether-axis")
	? axisFromName (vm.at("tether-axis").as<string>())
	: (board.zSize > 1 ? 2 : 1);
      const int strands = board.addTetheredArray (vm.at("init").as<string>(),
						  vm.at("tether").as<int>(), axis);
      cerr << "Tethered " << strands << " strands of "
	   << vm.at("init").as<string>().length() << " bases" << endl;
    } else if (vm.count("init")) {
      board.addSeq (vm.at("init").as<string>());
      if (vm.count("pin-first") && board.unit.size()) {
	board.unit[0].pinned = true;
	board.rebuildMobile();
	cerr << "Pinned unit 0 at " << board.unit[0].pos << endl;
      }
    }

    if (vm.count("density"))
      board.addBases (vm.at("density").as<double>(), rng);

    if (vm.count("sort-units"))
      board.sortUnitsByPosition();

    // logging
    const long logPeriod = vm.at("period").as<long>();
    const bool logColors = !vm.count("monochrome");
    const bool logSeqs = vm.count("seqs");
    const bool countPairs = vm.count("bitmap") || vm.count("csv") || vm.count("json");
    bool logFolds = vm.count("folds");
    if (logFolds && !board.isLinear()) {
      cerr << "Board is not a single linear chain, so fold strings are omitted from the log." << endl;
      logFolds = false;
    }
    const bool logTime = board.params.secondsPerSweep > 0;

    // do the simulation
    const long moves = vm.at("total-moves").as<long>()
      + (long) board.mobile.size() * vm.at("unit-moves").as<long>();
    long succeeded = 0, samples = 0;
    PairCounter pairCount (board.unit.size());
    for (long move = 0; move < moves; ++move) {
      if (board.tryMove (rng))
	++succeeded;
      if (move % logPeriod == 0) {
	const vguard<Board::IndexPair> pairs = (logFolds || countPairs)
	  ? board.indexPairs() : vguard<Board::IndexPair>();
	if (logFolds)
	  cout << succeeded
	       << " (" << fixed << setprecision(1) << (100. * move / moves) << "%) "
	       << (logColors ? board.coloredFoldString(pairs) : board.foldString(pairs))
	       << " " << setw(5) << board.foldEnergy(pairs)
	       << " " << setw(5) << board.unitRadiusOfGyration()
	       << " (" << to_string_join (board.unitCentroid()) << ")"
	       << " sweeps=" << setprecision(3) << board.sweeps()
	       << (logTime ? (" t=" + to_string (board.seconds()) + "s") : string())
	       << endl;
	if (logSeqs) {
	  const auto seqFreqs = board.sequenceFreqs();
	  cout << succeeded
	       << " (" << fixed << setprecision(1) << (100. * move / moves) << "%)"
	       << " sweeps=" << setprecision(3) << board.sweeps();
	  if (logTime) cout << " t=" << board.seconds() << "s";
	  for (auto& sf: seqFreqs)
	    cout << " " << sf.first << "(" << sf.second << ")";
	  cout << endl;
	}
	if (countPairs)
	  for (const auto& ij: pairs)
	    pairCount.add (ij.first, ij.second);
	++samples;
      }
    }

    // report results
    if (moves) {
      cerr << "Tried " << moves << " moves, " << succeeded << " succeeded" << endl;
      cerr << "Simulated " << board.sweeps() << " sweeps";
      if (logTime) cerr << " = " << board.seconds() << " s";
      cerr << " over " << board.mobile.size() << " mobile units" << endl;
    }

    if (vm.count("bitmap")) {
      if (board.unit.size() > 8192)
	cerr << "Too many units for a bitmap; skipping." << endl;
      else {
	bitmap_image image (board.unit.size(), board.unit.size());
	image.set_all_channels (0, 0, 0);
	pairCount.forEach ([&](int i, int j, uint32_t n) {
	    const int level = samples ? (int) ((255 * (uint64_t) n + 128) / samples) : 0;
	    image.set_pixel (i, j, level, level, level);
	  });
	image.save_image (vm.at("bitmap").as<string>().c_str());
      }
    }

    if (vm.count("csv")) {
      vguard<vguard<string> > pp (board.unit.size(), vguard<string> (board.unit.size()));
      pairCount.forEach ([&](int i, int j, uint32_t n) {
	  pp[i][j] = to_string (n / (double) samples);
	});
      ofstream outfile (vm.at("csv").as<string>());
      if (!outfile)
	throw runtime_error ("Can't save basepair probabilities to CSV file");
      outfile << "*," << join (board.sequence(), ",") << endl;
      for (size_t n = 0; n < pp.size(); ++n)
	outfile << board.unit[n].baseChar() << "," << to_string_join (pp[n], ",") << endl;
    }

    if (vm.count("json")) {
      json js;
      js["samples"] = samples;
      js["sequence"] = board.sequence();
      js["sweeps"] = board.sweeps();
      if (logTime) js["seconds"] = board.seconds();
      pairCount.forEach ([&](int i, int j, uint32_t n) {
	  js["prob"][to_string(i)][to_string(j)] = ((double) n) / samples;
	});
      ofstream outfile (vm.at("json").as<string>());
      if (!outfile)
	throw runtime_error ("Can't save basepair probabilities to JSON file");
      outfile << js << endl;
    }

    if (vm.count("save")) {
      json j = board.toJson();
      ofstream outfile (vm.at("save").as<string>());
      if (!outfile)
	throw runtime_error ("Can't save board file");
      outfile << j << endl;
    } else if (!logFolds && !logSeqs)
      cout << board.toJson() << endl;

  } catch (const exception& e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
