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
#include <api/BamReader.h>
#include <api/BamWriter.h>
#include "read.h"
#include "bam_utils.h" 
using namespace std;
using namespace BamTools;
using boost::format;
namespace fs = boost::filesystem;
namespace po = boost::program_options;

int main(int argc, char *argv[]) {
  string f_bam, f_win, f_out;
  clock_t time1 = clock();
  po::options_description cmdOpts("Allowed options");
  cmdOpts.add_options()
    ("help,h", "produce help message")
    ("f_bam,b", po::value<string>(&f_bam), "input (BAM) file")
    ("f_win,w", po::value<string>(&f_win), "window file")
    ("f_out,o", po::value<string>(&f_out), "output (tabular) file")
  ;
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, cmdOpts), vm);
  po::notify(vm);
  if(vm.count("help") || !vm.count("f_bam") || !vm.count("f_win") || !vm.count("f_out")) {
    cout << cmdOpts << endl;
    exit(1);
  }

  BamReader reader;
  Utilities utils;
  BamRegion region;
  if( !reader.Open(f_bam) ) {
    cout << format("cannot read from %s ...\n") % f_bam;
    exit(1);
  }

  vector<Loc> lv = read_locs(f_win);
  ofstream fho( f_out.c_str() );
  fho << "chr\tbeg\tend\tn_orphan\t\n"; 
  for(vector<Loc>::iterator it = lv.begin(); it != lv.end(); it ++) {
    string id = it->id, chr = it->chr;
    int beg = it->beg, end = it->end;
    string regionStr = str(format("%s:%d-%d") % chr % beg % end);
    cout << regionStr << endl;
    if( !utils.ParseRegionString(regionStr, reader, region) ) {
      cerr << format("invalid region string: %s\n") % regionStr;
      exit(1);
    } 
    reader.SetRegion(region);
    BamAlignment al;
    while( reader.GetNextAlignment(al) ) {
      
    }
  }
  fho.close();
  reader.Close();

  cout << right << setw(60) << format("(running time: %.01f minutes)\n") % ( (double)(clock() - time1) / ((double)CLOCKS_PER_SEC * 60) );
  return EXIT_SUCCESS;
}

