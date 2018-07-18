#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include <time.h>
#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/program_options.hpp>
#include "ssp.h" 
using namespace std;
using boost::format;
using namespace Sequence;
namespace po = boost::program_options;

int main(int argc, char *argv[]) {
  string fi, fo, og;
  po::options_description cmdOpts("Allowed options");
  cmdOpts.add_options()
    ("help,h", "produce help message")
    ("in,i", po::value<string>(&fi)->implicit_value(""), "input (SimpleSNP format)")
    ("out,o", po::value<string>(&fo)->implicit_value(""), "output")
    ("outgroup,g", po::value<string>(&og), "outgroup")
  ;

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, cmdOpts), vm);
  po::notify(vm);

  if(vm.count("help") || !vm.count("outgroup")) {
    cout << cmdOpts << endl;
    return 1;
  }
  
  istream* in = (!vm.count("in") || fi == "" || fi == "-" || fi == "stdin") ?
    &cin : new ifstream( fi.c_str() );
  ostream* out = (!vm.count("out") || fo == "" || fo == "-" || fo == "stdout") ?
    &cout : new ofstream( fo.c_str() );
   
  SSP ssp;
  ssp.read(*in);

  cerr << "Before SetOG: " << endl;
  cerr << format("%d samples, %6d positions, outgroup[%s]\n") % ssp.nind % ssp.npos % (ssp.outgroup() ? "yes" : "no");
  cerr << format("...now SetOG to %s") % og << endl;

  vector<double> poss = ssp.GetPositions();
  vector<string> data = ssp.GetData();
  unsigned idx_og = ssp.get_label_idx(og);
  
  ssp.labels.erase(ssp.labels.begin() + idx_og);
  ssp.labels.insert(ssp.labels.begin(), og);
  
  string data_og = data[idx_og];
  data.erase(data.begin() + idx_og);
  data.insert(data.begin(), data_og);
  
  ssp.assign(&poss[0], poss.size(), &data[0], data.size());
  if( !ssp.outgroup() ) ssp.set_outgroup(true);
  ssp.nind = data.size();
  
  unsigned n_N = 0, n_fix = 0;
  vector<double> poss_n;
  vector<string> data_n(ssp.nind);
  for(unsigned j = 0; j < ssp.npos; ++ j) {
    char anc = data[0][j];
    if(anc == 'N') {
      n_N ++;
    } else {
      unsigned n_anc = 0;
      for(unsigned i = 1; i < ssp.nind; ++ i) {
        if(data[i][j] == anc) n_anc ++;
      }
      if( n_anc == 0 ) {
        n_fix ++;
      } else {
        poss_n.push_back(poss[j]);
        for(unsigned i = 0; i < ssp.nind; ++ i) 
          data_n[i] += data[i][j];
      }
    }
  }
  ssp.npos = poss_n.size();
  ssp.nind = data_n.size();
  ssp.assign(&poss_n[0], poss_n.size(), &data_n[0], data_n.size());
  
  cerr << format("   %6d positions removed: %s = 'N'\n") % n_N % og;
  cerr << format("   %6d positions removed: fixed in ingroup\n") % n_fix;
  cerr << format("%d samples, %6d positions, outgroup[%s]\n") % ssp.nind % ssp.npos % (ssp.outgroup() ? "yes" : "no");
  ssp.write(*out, "ssp");

  return EXIT_SUCCESS;
}



