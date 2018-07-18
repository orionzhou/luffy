#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include <time.h>
#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/program_options.hpp>
#include "snp.h" 
using boost::format;
using namespace std;
using namespace Sequence;
namespace po = boost::program_options;

vector<string> read_mt_ids(const string& f_id, const string& opt) {
  ifstream fhi(f_id.c_str());
  if(!fhi.is_open()) {
    cerr << format("cannot read file: %s\n") % f_id; 
    exit(1);
  }
  
  string line;
  vector<string> ss;
  int idx = -1;

  getline(fhi, line);
  boost::trim_if(line, boost::is_any_of("\t "));
  boost::split(ss, line, boost::is_any_of("\t "), boost::token_compress_on);
  for(uint32_t i = 0; i < ss.size(); i ++)
    if(ss[i] == opt) idx = i;
  if(idx == -1) {
    cerr << "unknown opt: " << opt << endl;
    exit(1);
  }
  
  vector<string> ids;
  while(fhi.good()) {
    getline(fhi, line);
    boost::split(ss, line, boost::is_any_of("\t "));
    if(ss[idx] == "1") ids.push_back(ss[0]);
  }
  return ids;
}

int main(int argc, char *argv[]) {
  string fi, fo, fa, opt;
  po::options_description cmdOpts("Allowed options");
  cmdOpts.add_options()
    ("help,h", "produce help message")
    ("in,i", po::value<string>(&fi)->implicit_value(""), "input (SimpleSNP format)")
    ("out,o", po::value<string>(&fo)->implicit_value(""), "output")
    ("acc,a", po::value<string>(&fa)->default_value("/home/youngn/zhoup/git/luffy/conf/acc_ids.tbl"), "acc option file")
    ("opt,p", po::value<string>(&opt), "id option")
  ;
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, cmdOpts), vm);
  po::notify(vm);
  if(vm.count("help") || !vm.count("opt")) {
    cout << cmdOpts << endl;
    return 1;
  }

  istream* in = (!vm.count("in") || fi == "" || fi == "-" || fi == "stdin") ?
    &cin : new ifstream( fi.c_str() );
  ostream* out = (!vm.count("out") || fo == "" || fo == "-" || fo == "stdout") ?
    &cout : new ofstream( fo.c_str() );
   
  vector<string> inds = read_mt_ids(fa, opt);
  
  SNP snp;
  snp.read( *in );

  unsigned nind(snp.get_nind()), npos(snp.get_npos());
  cerr << format("before ind_select: %3d inds, %6d loci\n") % nind % npos;
  
  vector<unsigned> idxs_pos;
  vector<unsigned> idxs_ind = snp.get_labels_idx(inds);
  cerr << format("  selected: %s\n") % boost::join(inds, " ");
  
  snp.subset(idxs_ind, idxs_pos);
  snp.RemoveMono();

  nind = snp.get_nind();
  npos = snp.get_npos();
  cerr << format("after ind_select:  %3d inds, %6d loci\n") % nind % npos;

  snp.write(*out, "snp");

  return EXIT_SUCCESS;
}



