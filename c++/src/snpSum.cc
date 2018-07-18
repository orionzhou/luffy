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
  cmdOpts.add_options()
    ("help,h", "produce help message")
    ("in,i",  po::value<string>(&fi)->implicit_value(""), "input")
    ("out,o", po::value<string>(&fo)->implicit_value(""), "output")
    ("iformat,x",  po::value<string>(&ifmt)->default_value("snp"), "input format")
    ("oformat,y",  po::value<string>(&ofmt)->default_value("snp"), "output format")
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
   
  unsigned nind (snp.get_nind()), npos (snp.get_npos());
  vector<string> labels (snp.get_labels());
  bool haveOG (snp.get_haveOG());

  *out << nind << " samples (outgroup: " << (haveOG ? "yes" : "no") << ")" << endl;
  *out << npos << " sites" << endl;
  *out << "  " << boost::algorithm::join(labels, ",") << endl;
  *out << "  " << snp.position(0) << "," << snp.position(1) << ",...," << snp.position(npos-1) << endl;
  
  cerr << right << setw(60) << format("(running time: %.01f minutes)\n") % ( (double)(clock() - time1) / ((double)CLOCKS_PER_SEC * 60) );
  return EXIT_SUCCESS;
}
