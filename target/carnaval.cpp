#include <time.h>
#include <cstdlib>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <random>
#include <boost/program_options.hpp>

#include "../src/cell.h"

using namespace std;
namespace po = boost::program_options;

int main (int argc, char** argv) {

  try {

    // Declare the supported options.
    po::options_description opts("Options");
    opts.add_options()
      ("help,h", "display this help message")
      ("xsize,x", po::value<int>()->default_value(64), "size of board in X dimension")
      ("ysize,y", po::value<int>()->default_value(64), "size of board in Y dimension")
      ("zsize,z", po::value<int>()->default_value(1), "size of board in Z dimension")
      ("init,i",  po::value<string>(), "specify initial sequence")
      ("rnd,r",  po::value<int>(), "seed random number generator")
      ("total-moves,t",  po::value<long>()->default_value(0), "total number of moves")
      ("unit-moves,u",  po::value<long>()->default_value(0), "number of moves per unit")
      ("folds,f", "report fold strings every 1000 moves")
      ("temp,T",  po::value<double>(), "specify temperature")
      ("load,l", po::value<string>(), "load board state from file")
      ("save,s", po::value<string>(), "save board state to file")
      ;

    po::variables_map vm;
    po::parsed_options parsed = po::command_line_parser(argc,argv).options(opts).allow_unregistered().run();
    po::store (parsed, vm);
    po::notify(vm);    
      
    // parse args
    if (vm.count("help")) {
      cout << opts << endl;
      return 1;
    }

    time_t timer;
    time (&timer);
    const int seed = vm.count("rnd") ? vm.at("rnd").as<int>() : timer;
    mt19937 mt (seed);
    cerr << "Random seed is " << seed << endl;

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
    }

    if (vm.count("init"))
      board.addSeq (vm.at("init").as<string>());

    if (vm.count("temp"))
      board.params.temp = vm.at("temp").as<double>();

    const long moves = vm.at("total-moves").as<long>() + board.unit.size() * vm.at("unit-moves").as<long>();
    const bool reportFolds = vm.count("folds");
    int succeeded = 0;
    for (long move = 0; move < moves; ++move) {
      if (board.tryMove (mt))
	++succeeded;
      if (reportFolds && move % 1000 == 0)
	cerr << "Move " << move << ": " << board.foldString() << endl;
    }

    if (moves)
      cerr << "Tried " << moves << " moves, " << succeeded << " succeeded" << endl;
    
    if (vm.count("save")) {
      json j = board.toJson();
      ofstream outfile (vm.at("save").as<string>());
      if (!outfile)
	throw runtime_error ("Can't save board file");
      outfile << j << endl;
    } else
      cout << board.toJson() << endl;

  } catch (const exception& e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  
  return EXIT_SUCCESS;
}
