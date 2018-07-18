#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include <time.h>
#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/program_options.hpp>
#include "snp.h" 
using namespace std;
using boost::format;
using namespace Sequence;
namespace po = boost::program_options;

int main(int argc, char *argv[]) {
  string fi, fo, ind1, ind2;
  po::options_description cmdOpts("Allowed options");
  cmdOpts.add_options()
    ("help,h", "produce help message")
    ("in,i", po::value<string>(&fi)->implicit_value(""), "input (SimpleSNP format)")
    ("out,o", po::value<string>(&fo)->implicit_value(""), "output")
    ("ind1,1", po::value<string>(&ind1), "individual 1")
    ("ind2,2", po::value<string>(&ind2), "individual 2")
  ;

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, cmdOpts), vm);
  po::notify(vm);

  if(vm.count("help") || !vm.count("ind1") || !vm.count("ind2")) {
    cout << cmdOpts << endl;
    return 1;
  }
  
  istream* in = (!vm.count("in") || fi == "" || fi == "-" || fi == "stdin") ?
    &cin : new ifstream( fi.c_str() );
  ostream* out = (!vm.count("out") || fo == "" || fo == "-" || fo == "stdout") ?
    &cout : new ofstream( fo.c_str() );
   
  SNP snp;
  snp.read( *in );
  
  unsigned nind (snp.get_nind()), npos (snp.get_npos());
  vector<string> labels (snp.get_labels());

  unsigned idx1 = snp.get_label_idx(ind1);
  unsigned idx2 = snp.get_label_idx(ind2);

  int n_same=0, n_diff=0, n_1n=0, n_2n=0, n_nn=0;
  vector<string> data = snp.GetData();
  for(unsigned j = 0; j < npos; ++j) {
    char s1 = char(std::toupper(data[idx1][j]));
    char s2 = char(std::toupper(data[idx2][j]));
    if(s1 == 'N' && s2 == 'N') {
      n_nn ++;
    } else if(s1 == 'N' && s2 != 'N') {
      n_1n ++;
    } else if(s1 != 'N' && s2 == 'N') {
      n_2n ++;
    } else if(s1 == s2) {
      n_same ++;
    } else {
      n_diff ++;
    }
  }

  *out << format("%3d inds, %6d loci") % nind % npos << endl;
  *out << format("%6d positions: %s  = %s") % n_same % ind1 % ind2 << endl;
  *out << format("%6d positions: %s != %s") % n_diff % ind1 % ind2 << endl;
  *out << format("%6d positions: %s  = N && %s != N") % n_1n % ind1 % ind2 << endl;
  *out << format("%6d positions: %s != N && %s  = N") % n_2n % ind1 % ind2 << endl;
  *out << format("%6d positions: %s  = N && %s  = N") % n_nn % ind1 % ind2 << endl;

  return EXIT_SUCCESS;
}



