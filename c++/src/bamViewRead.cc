#include <iostream>
#include <iomanip>
#include <vector>
#include <time.h>
#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <boost/program_options.hpp>
#include <api/BamReader.h>
#include <api/BamWriter.h>
#include "bam_utils.h"
using namespace std;
using boost::format;
namespace fs = boost::filesystem;
namespace po = boost::program_options;
using namespace BamTools;

int main(int argc, char *argv[]) {
  string fi, fo, dirW, regionStr, name;
  int n_more;
  clock_t time1 = clock();
  po::options_description cmdOpts("Allowed options");
  cmdOpts.add_options()
    ("help,h", "produce help message")
    ("in,i", po::value<string>(&fi), "input (BAM) file")
    ("name,n", po::value<string>(&name), "read name")
    ("region,r", po::value<string>(&regionStr), "region")
    ("more,m", po::value<int>(&n_more)->default_value(0), "number more reads to display")
    ("out,o", po::value<string>(&fo)->default_value(""), "output (BAM) file")
  ;

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, cmdOpts), vm);
  po::notify(vm);
  if(vm.count("help") || !vm.count("in") || !vm.count("name")) {
    cout << cmdOpts << endl;
    return 1;
  }

  BamReader reader;
  Utilities util;
  fs::path fi1( fi );
  fs::path fi2( fi + ".bai" );
  if( exists(fi1) ) {
    reader.Open(fi1.string());
    reader.OpenIndex(fi2.string());
    if( reader.HasIndex() ) {
      cout << format("%s opened with index ...\n") % fi;
      if( vm.count("region") ) {
        BamRegion region;
        if( util.ParseRegionString(regionStr, reader, region) ) {
          reader.SetRegion(region);
          cout << format("region set to: %s\n") % regionStr;
        } else {
          cerr << format("invalid region string: %s\n") % regionStr;
          exit(1);
        }
      } else {
        cout << "no region specified: starting from scratch\n";
      }
    } else {
      cout << format("%s opened without index ...\n") % fi;
    }
  } else {
    cerr << format("cannot read from %s ...\n") % fi;
    exit(1);
  }
  RefVector refs = reader.GetReferenceData();

  BamAlignment al, al1, al2;
  bool tag_f = false, tag_s = false;
  uint32_t cnt = 0;
  while( reader.GetNextAlignment(al) ) {
    cnt ++;
    if( al.Name == name && al.IsFirstMate() ) {
      tag_f = true;
      al1 = al;
    }
    if( al.Name == name && al.IsSecondMate() ) {
      tag_s = true;
      al2 = al;
    }
    if(tag_f == true && tag_s == true) break;
  }
  reader.Close();
  if(tag_f != true) {
    cerr << "cannot find first mate\n";
    util.PrintRead(al2, refs);
    exit(1);
  } else if(tag_s != true) {
    cerr << "cannot find second mate\n";
    util.PrintRead(al1, refs);
    exit(1);
  } else {
    util.PrintRead(al1, refs);
    util.PrintRead(al2, refs);
  }
  if( fo != "" ) {
    BamWriter writer;
    writer.Open(fo, reader.GetHeader(), reader.GetReferenceData());
    writer.SaveAlignment(al1);
    writer.SaveAlignment(al2);
    writer.Close();
  }

  cerr << format("%10.0f alignments read\n") % cnt;
  cerr << right << setw(60) << format("(running time: %.01f minutes)\n") % 
    ( (double)(clock() - time1) / ((double)CLOCKS_PER_SEC * 60) );
  return EXIT_SUCCESS;
}


 
