#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include <time.h>
#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/program_options.hpp>
#include <Sequence/stateCounter.hpp>
#include "snp.h" 
using namespace std;
using boost::format;
using namespace Sequence;
namespace po = boost::program_options;

int main(int argc, char *argv[]) {
  string fi, fo, ind_ref, chr;
  po::options_description cmdOpts("Allowed options");
  cmdOpts.add_options()
    ("help,h", "produce help message")
    ("in,i", po::value<string>(&fi), "input (simpleSNP) file")
    ("out,o", po::value<string>(&fo), "output (VCF) file")
    ("ref,r", po::value<string>(&ind_ref)->default_value(""), "reference ind")
    ("chr,c", po::value<string>(&chr)->default_value("chrU"), "chr")
  ;
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, cmdOpts), vm);
  po::notify(vm);
  if(vm.count("help") || !vm.count("in") || !vm.count("out")) {
    cout << cmdOpts << endl;
    return 1;
  }

  istream* in = (!vm.count("in") || fi == "" || fi == "-" || fi == "stdin") ?
    &cin : new ifstream( fi.c_str() );
  ostream* out = (!vm.count("out") || fo == "" || fo == "-" || fo == "stdout") ?
    &cout : new ofstream( fo.c_str() );
  
  SNP snp;
  snp.read(*in);

  unsigned nind (snp.get_nind()), npos (snp.get_npos());
  cout << format("%3d inds, %6d loci\n") % nind % npos;
  
  vector<string> data = snp.GetData();
  vector<double> poss = snp.GetPositions();
  unsigned idx_ref = snp.get_label_idx(ind_ref);

  for(unsigned j = 0; j < npos; ++ j) {
    char ref = data[idx_ref][j];
    if(ref == 'N') {
      cerr << format("position %d skipped : ref = N\n") % poss[j];
      continue;
    }
    stateCounter Counts;
    for(unsigned i = 0; i < nind; ++ i) Counts( data[i][j] );
    if(ref != 'A' && Counts.a > 0)
      *out << format("%s\t%.0f\t%s\t%s\t%s\n") % chr % poss[j] % "" % ref % "A";
    if(ref != 'T' && Counts.t > 0)
      *out << format("%s\t%.0f\t%s\t%s\t%s\n") % chr % poss[j] % "" % ref % "T";
    if(ref != 'C' && Counts.c > 0) 
      *out << format("%s\t%.0f\t%s\t%s\t%s\n") % chr % poss[j] % "" % ref % "C";
    if(ref != 'G' && Counts.g > 0)
      *out << format("%s\t%.0f\t%s\t%s\t%s\n") % chr % poss[j] % "" % ref % "G";
  }
  out->flush();

  return EXIT_SUCCESS;
}



