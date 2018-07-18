#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include <vector>
#include <map>
#include <time.h>
#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <boost/program_options.hpp>
#include <boost/regex.hpp>
#include <boost/algorithm/string.hpp>
#include <geco.h>
#include "read.h"
#include "ssp.h"
using namespace std;
using namespace geco;
using boost::format;
namespace fs = boost::filesystem;
namespace po = boost::program_options;

int main( int argc, char* argv[] ) {
  string f_in, f_out, f_acc, opt_ind;
  int opt_conf;
  clock_t time1 = clock();
  po::options_description cmdOpts("Allowed options");
  cmdOpts.add_options()
    ("help,h", "produce help message")
    ("f_in,i", po::value<string>(&f_in), "window file")
    ("f_out,o", po::value<string>(&f_out), "output (cov) file")
    ("f_acc,a", po::value<string>(&f_acc)->default_value("/project/youngn/zhoup/Scripts/conf/acc_ids.tbl"), "acc option file")
    ("opt_ind,t", po::value<string>(&opt_ind), "id option")
    ("opt_conf,c", po::value<int>(&opt_conf)->default_value(1), "id_all option")
  ;
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, cmdOpts), vm);
  po::notify(vm);
  if(vm.count("help") || !vm.count("opt_ind") || !vm.count("f_in") || !vm.count("f_id") || !vm.count("f_out")) {
    cout << cmdOpts << endl;
    return 1;
  }

  IntVec idxs;
  StrVec inds;
  idxs = get_acc_idx(f_acc, opt_ind, inds, opt_conf);
  
  string f_cov1, f_cov2;
  if(opt_conf == 1) {
    f_cov1 = "/project/youngn/zhoup/Data/repo/mt_35/30_vnt_acc84/21_coverage/cov1.tbl.gz";
    f_cov2 = "/project/youngn/zhoup/Data/repo/mt_35/30_vnt_acc84/21_coverage/cov2.tbl.gz";
  } else if(opt_conf == 2) {
    f_cov1 = "/project/youngn/zhoup/Data/repo/mt_35/30_vnt_acc288/21_coverage/cov1.tbl.gz";
    f_cov2 = "/project/youngn/zhoup/Data/repo/mt_35/30_vnt_acc288/21_coverage/cov2.tbl.gz";
  } else {
    cout << "unknown opt_conf: " << opt_conf << endl;
    exit(1);
  }
  
  ofstream fho( f_out.c_str() );
  fho << "id\tchr\tbeg\tend\tacc\tcov1_mean\tcov1_max\tcov1_sd\tcov2_mean\tcov2_max\tcov2_sd\n"; 
  
  vector<Loc> lv = read_locs(f_in);
  for(vector<Loc>::iterator it = lv.begin(); it != lv.end(); it ++) {
    string id = it->id, chr = it->chr;
    int beg = it->beg, end = it->end;
    int nrow = end - beg + 1;
    cout << format("%s:%d-%d\n") % chr % beg % end;

    IntVecVec ivv1 = get_raw_cov(chr, beg, end, f_cov1, opt_conf);
    IntVecVec ivv2 = get_raw_cov(chr, beg, end, f_cov2, opt_conf);
    
    for(uint32_t i = 0; i < idxs.size(); i ++) {
      int idx = idxs[i];
      string ind = inds[i];
      IntVec iv1 = ivv1[idx], iv2 = ivv2[idx];
      
      gArray<gScore> covary1(0, nrow, false, false);
      gArray<gScore> covary2(0, nrow, false, false);
      for(int j = 0; j < nrow; j ++) {
        covary1.setValue(j, iv1[j]);
        covary2.setValue(j, iv2[j]);
      }
      
      gScore cov1_max = covary1.getMax()[0];
      gScore cov1_mean = covary1.getMean()[0];
      gScore cov1_sd = covary1.getStdDev()[0];
      gScore cov2_max = covary2.getMax()[0];
      gScore cov2_mean = covary2.getMean()[0];
      gScore cov2_sd = covary2.getStdDev()[0];

      fho << format("%s\t%s\t%d\t%d\t%s\t%.03f\t%d\t%.03f\t%.03f\t%d\t%.03f\n") % id % chr % beg % end % ind % cov1_mean % cov1_max % cov1_sd % cov2_mean % cov2_max % cov2_sd;
    }
  }

  cout << right << setw(60) << format("(running time: %.01f minutes)\n") % ( (double)(clock() - time1) / ((double)CLOCKS_PER_SEC * 60) );
  return EXIT_SUCCESS;
}


