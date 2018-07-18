#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <vector>
#include <map>
#include <time.h>
#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <boost/program_options.hpp>
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
  string fi, fo, regionStr;
  clock_t time1 = clock();
  po::options_description cmdOpts("Allowed options");
  cmdOpts.add_options()
    ("help,h", "produce help message")
    ("in,i", po::value<string>(&fi), "input (BAM) file")
    ("out,o", po::value<string>(&fo), "output file")
    ("region,r", po::value<string>(&regionStr), "region")
  ;
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, cmdOpts), vm);
  po::notify(vm);
  if(vm.count("help") || !vm.count("in") || !vm.count("out")) {
    cout << cmdOpts << endl;
    return 1;
  }

  Utilities util;
  BamReader reader;
  fs::path fi1( fi );
  fs::path fi2( fi + ".bai" );
  if( exists(fi1) ) {
    reader.Open(fi1.string());
    reader.OpenIndex(fi2.string());
    if( reader.HasIndex() ) {
      cerr << format("%s opened with index ...\n") % fi;
      if( vm.count("region") ) {
        BamRegion region;
        if( util.ParseRegionString(regionStr, reader, region) ) {
          reader.SetRegion(region);
          cerr << format("region set to: %s\n") % regionStr;
        } else {
          cerr << format("invalid region string: %s\n") % regionStr;
          exit(1);
        }
      } else {
        cerr << "no region specified: starting from scratch\n";
      }
    } else {
      cerr << format("%s opened without index ...\n") % fi;
    }
  } else {
    cerr << format("cannot read from %s ...\n") % fi;
    exit(1);
  }
  
  RefVector refs = reader.GetReferenceData();
  BamWriter writer;
  writer.Open(fo, reader.GetHeader(), reader.GetReferenceData());

  BamAlignment al;
  string id;
  int type;
  int32_t is;
  vector<CigarOp> cigars;
  while( reader.GetNextAlignment(al) ) {
    if( al.IsDuplicate() || !al.IsMapped() ) continue;
    id = al.Name;
    cigars = al.CigarData;
    type = 0; // type = 1(orphan) 2(abnormal insert size) 3(soft-clip)
    is = al.InsertSize;
    
    if( !al.IsMateMapped() ) {
      type = 1;
    } else if( al.RefID != al.MateRefID || is > 10000 || is < -10000) {
      type = 2;
    } else {
      char typef(cigars[0].Type), typel(cigars[cigars.size()-1].Type);
      uint32_t lenf(cigars[0].Length), lenl(cigars[cigars.size()-1].Length);
      if( (typef == 'S' && lenf >= 10) || (typel == 'S' && lenl >= 10) )
        type = 3;
    }
    if(type > 0) writer.SaveAlignment(al);
  }
  reader.Close();
  writer.Close();

  cerr << right << setw(60) << format("(running time: %.01f minutes)\n") % 
    ( (double)(clock() - time1) / ((double)CLOCKS_PER_SEC * 60) );
  return EXIT_SUCCESS;
}
