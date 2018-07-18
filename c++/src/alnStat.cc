#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <vector>
#include <time.h>
#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <boost/program_options.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/progress.hpp>
#include <Sequence/Fasta.hpp>
#include <Sequence/Alignment.hpp>
#include <Sequence/SimpleSNP.hpp>
#include <Sequence/PolySites.hpp>
#include <Sequence/PolySNP.hpp>
#include <Sequence/PolyTableSlice.hpp>
#include <Sequence/SeqExceptions.hpp>
#include "read.h"
using namespace std;
using boost::format;
using namespace Sequence;
namespace fs = boost::filesystem;
namespace po = boost::program_options;

int main(int argc, char *argv[]) {
  string f_id, d_in, f_out;
  clock_t time1 = clock();
  po::options_description cmdOpts("Allowed options");
  cmdOpts.add_options()
    ("help,h", "produce help message")
    ("f_id,i", po::value<string>(&f_id), "input (gene id) file")
    ("d_in,d", po::value<string>(&d_in), "input (fasta) directory")
    ("f_out,o", po::value<string>(&f_out), "output (gene statistics) file")
  ;
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, cmdOpts), vm);
  po::notify(vm);

  if(vm.count("help") || !vm.count("f_id") || !vm.count("d_in") || !vm.count("f_out")) {
    cout << cmdOpts << endl;
    return 1;
  }

  set<string> ids = read_ids(f_id);
  ofstream fho(f_out.c_str());
  if (!fho.is_open()) { cout << format("cannot write output: %s\n") % f_out; return 1; }

  for(set<string>::iterator it = ids.begin(); it != ids.end(); it ++) {
    string id = *it;
    string f_in = d_in + "/cds/" + id + ".fa";
    vector<Fasta> aln;
    Alignment::GetData(aln, f_in.c_str());
    PolySites pt( aln );
    PolySNP ps( &pt );
    cout << format("%s\t%d\t%d\t%.03f\t%.03f\n") % id % pt.numsites() % ps.NumPoly() % ps.ThetaW() % ps.ThetaPi();
  }
  fho.close();
  cout << right << setw(60) << format("(running time: %.01f minutes)\n") % ( (double)(clock() - time1) / ((double)CLOCKS_PER_SEC * 60) );
  return EXIT_SUCCESS;
}



