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
  string fi, fo, ifmt, ofmt;
  po::options_description cmdOpts("Allowed options");
  int co_maf, co_mis;
  bool co_bia;
  cmdOpts.add_options()
    ("help,h", "produce help message")
    ("in,i",  po::value<string>(&fi)->implicit_value(""), "input")
    ("out,o", po::value<string>(&fo)->implicit_value(""), "output")
    ("iformat,x",  po::value<string>(&ifmt)->default_value("snp"), "input format")
    ("oformat,y",  po::value<string>(&ofmt)->default_value("snp"), "output format")
    ("co_mis,m", po::value<int>(&co_mis)->default_value(-1), "missing data cutoff")
    ("co_maf,a", po::value<int>(&co_maf)->default_value(0), "MAF minimum cutoff")
    ("co_bia,b", po::value<bool>(&co_bia)->default_value(0), "only keep bi-allelic positions")
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

  clock_t time1 = clock();
  SNP snp;
  snp.import(*in, ifmt);
  
  cout << format("before filter: %3d inds, %6d loci\n") % snp.get_nind() % snp.get_npos();
  co_mis = (co_mis == -1) ? snp.get_nind() : co_mis;

  snp.FilterMissing(co_mis);
  cout << format("after filter [co_mis=%3d]: %6d positions\n") % co_mis % snp.get_npos();

  snp.ApplyFreqFilter(co_maf);
  cout << format("after filter [co_maf=%3d]: %6d positions\n") % co_maf % snp.get_npos();

  if(co_bia) {
    snp.RemoveMultiHits();
    cout << format("after filter [co_bia=%3d]: %6d positions\n") % co_bia % snp.get_npos();
  }
  	
  snp.write(*out, ofmt);

  cerr << right << setw(60) << format("(running time: %.01f minutes)\n") % ( (double)(clock() - time1) / ((double)CLOCKS_PER_SEC * 60) );
  return EXIT_SUCCESS;
}



