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

vector<uint32_t> sample_serial(const uint32_t& n, uint32_t& m) {
  vector<uint32_t> iv;
  if(m > n || m <= 0) m = n;
  double increment = n / m;
  for(uint32_t i = 0; i < m; i ++) 
    iv.push_back(1 + int(i * increment) - 1);
  return iv;
}
vector<uint32_t> sample_random(const uint32_t& n, const uint32_t& m) {
  vector<uint32_t> iv;
  for(uint32_t i = 0; i < m; i ++) 
    iv.push_back(1 + int(n * rand()/(RAND_MAX+1.0)) - 1);
  return iv;
}

int main(int argc, char *argv[]) {
  string fi, fo, ifmt, ofmt;
  uint32_t n_sam;
  clock_t time1 = clock();
  po::options_description cmdOpts("Allowed options");
  cmdOpts.add_options()
    ("help,h", "produce help message")
    ("in,i",  po::value<string>(&fi)->implicit_value(""), "input")
    ("out,o", po::value<string>(&fo)->implicit_value(""), "output")
    ("iformat,x",  po::value<string>(&ifmt)->default_value("snp"), "input format")
    ("oformat,y",  po::value<string>(&ofmt)->default_value("snp"), "output format")
    ("sample,n", po::value<uint32_t>(&n_sam)->default_value(0), "# positions to sample")
  ;

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, cmdOpts), vm);
  po::notify(vm);
  if(vm.count("help")) {
    cout << cmdOpts << endl;
    return 1;
  }

  istream* in = (!vm.count("in") || fi == "" || fi == "-" || fi == "stdin") ?
    &cin : new ifstream( fi.c_str() );
  ostream* out = (!vm.count("out") || fo == "" || fo == "-" || fo == "stdout") ?
    &cout : new ofstream( fo.c_str() );
  
  SNP snp;
  snp.import(*in, ifmt);

  unsigned nind (snp.get_nind()), npos (snp.get_npos());
  cout << format("before sampling: %3d inds, %6d loci\n") % nind % npos;
  n_sam = (n_sam == 0) ? snp.get_npos() : n_sam;

  vector<uint32_t> idxs_ind;
  vector<uint32_t> idxs_pos = sample_serial(npos, n_sam);
  snp.subset(idxs_ind, idxs_pos);
  
  nind = snp.get_nind();
  npos = snp.get_npos();
  cout << format("after sampling: %3d inds, %6d loci\n") % nind % npos;

  snp.write(*out, ofmt);

  cout << right << setw(60) << format("(running time: %.01f minutes)\n") % ( (double)(clock() - time1) / ((double)CLOCKS_PER_SEC * 60) );
  return EXIT_SUCCESS;
}


