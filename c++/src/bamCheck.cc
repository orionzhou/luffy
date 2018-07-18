#include <iostream>
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
#include "bam_utils.h" 

using namespace std;
using boost::format;
namespace fs = boost::filesystem;
namespace po = boost::program_options;
using namespace BamTools;

int main(int argc, char *argv[]) {
  //define command-line arguments
  string fi;
  clock_t time1 = clock();
  po::options_description cmdOpts("Allowed options");
  cmdOpts.add_options()
    ("help,h", "produce help message")
    ("in,i", po::value<string>(&fi), "input file")
  ;
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, cmdOpts), vm);
  po::notify(vm);
  if(vm.count("help") || !vm.count("in")) {
    cout << cmdOpts << endl;
    return 1;
  }

  Utilities util;
  BamReader2 reader;
  BamWriter writer;
  if( reader.Open(fi) ) {
    cout << format("header read from %s ...\n") % fi;
  } else {
    cout << format("cannot read from %s ...\n") % fi;
    exit(1);
  }

  int cnt_all(0), cnt_unmapped(0), cnt_orphan(0), cnt_orphan_u(0);
  int cnt_mapped(0), cnt_mapped_u(0), cnt_mapped_up(0);
  map<int32_t, uint32_t> m1;
  map<int32_t, uint32_t>::iterator it1;
  pair< map<int32_t, uint32_t>::iterator, bool > p1;

  string name;
  BamAlnVector als, als2;
  BamAlnPairVector alnpairs;
  vector<int32_t> iss;
  int n_total, type;

  while( reader.GetReads(als, name, als2) ) {
    if( !util.FixReads(name, als, type, n_total, iss, alnpairs) ) {
      cout << format("error fixing/paring reads: %s\n") % name;
      util.PrintReads(name, n_total, als);
      cout << endl;
      util.PrintReads(name, n_total, alnpairs);
      exit(0);
    }
//    cout << format("%d alnpairs for %s\n") % alnpairs.size() % name;
    for (int i = 0; i < n_total; i ++) {
      if(i == 0) {
        cnt_all ++;
        if( type == 3 ) {
          cnt_unmapped ++;
        } else if( type == 1 ) {
          cnt_mapped ++;
          if( n_total == 1 ) {
            cnt_mapped_u ++;
            if( iss[i] >= 90 && iss[i] <= 1000 ) {
              p1 = m1.insert( pair<int32_t, uint32_t> (iss[i], 1) );
              if(p1.second == false) p1.first->second ++;
              cnt_mapped_up ++;
            }
          }
        } else {
          cnt_orphan ++;
          if( n_total == 1 ) cnt_orphan_u ++;
        }
      }
    }
  }
  reader.Close();

  cout << format("%d\t%d\t%d\t%d\t%d\t%d\t%d\n") % cnt_all % cnt_unmapped 
    % cnt_orphan % cnt_mapped % cnt_orphan_u % cnt_mapped_u % cnt_mapped_up;

  cout << right << setw(60) << format("(running time: %.01f minutes)\n") % ( (double)(clock() - time1) / ((double)CLOCKS_PER_SEC * 60) );
  return EXIT_SUCCESS;
}
