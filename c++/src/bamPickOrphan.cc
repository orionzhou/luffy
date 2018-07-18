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

  ofstream fho( fo.c_str() );
  fho << "id\tmate\ttype\tchr\tbeg\tend\tstrd\tcigar\tmq\tis\trg\tseq\n";

  BamAlignment al;
  string id="", chr="", strd="", cigar="", rg="", seq="";
  int mate, type, beg, end, mq, is;
  vector<CigarOp> cigars;
  while( reader.GetNextAlignment(al) ) {
//    if( al.IsDuplicate() ) continue;
    id = al.Name;
    mate = al.IsFirstMate() ? 1 : 2;
    cigars = al.CigarData;
    type = 0; // type = 1(orphan) 2(inter-chromosomal) 3(soft-clip)
    seq = al.QueryBases;
    
    if(al.IsMapped()) chr = refs[al.RefID].RefName;
    beg = al.Position;
    end = al.GetEndPosition();
    strd = al.IsReverseStrand() ? "-" : "+";
    cigar = util.getCigarString(cigars);
    mq = al.MapQuality;
    is = al.InsertSize;
    al.GetTag("RG", rg);
    
    if( (al.IsMapped() && !al.IsMateMapped()) || 
      (!al.IsMapped() && al.IsMateMapped()) ) {
      type = 1;
      if( al.IsMapped() ) {
        fho << format("%d\t%d\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s\n") 
          % id % mate % type % chr % beg % end % strd % cigar % mq % "" % "" % "" ;
      } else {
        fho << format("%d\t%d\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s\n") 
          % id % mate % type % ""  % ""  % ""  % ""   % ""  % "" % "" % rg % seq;
      }
    } else if( al.IsMapped() && al.IsMateMapped() ) {
      if ( al.RefID != al.MateRefID ) {
        type = 2;
        fho << format("%d\t%d\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s\n") 
          % id % mate % type % chr % beg % end % strd % cigar % mq % "" % rg % seq;
      } else if( (cigars[0].Type == 'S' && cigars[0].Length >= 10) ||
  (cigars[cigars.size()-1].Type == 'S' && cigars[cigars.size()-1].Length >= 10) ) {
        type = 3;
        
        fho << format("%d\t%d\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s\n") 
          % id % mate % type % ""  % ""  % ""  % ""   % cigar % "" % "" % rg % seq;

        mate = al.IsFirstMate() ? 2 : 1;
        chr = refs[al.MateRefID].RefName;
        beg = al.MatePosition;
        end = al.MatePosition + al.Length;
        strd = al.IsMateReverseStrand() ? "-" : "+";
        mq = 60;
        is = -al.InsertSize;
        fho << format("%d\t%d\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s\n") 
          % id % mate % type % chr % beg % end % strd % ""  % mq % is % "" % "" ;
      }
    }
  }
  reader.Close();

  cerr << right << setw(60) << format("(running time: %.01f minutes)\n") % 
    ( (double)(clock() - time1) / ((double)CLOCKS_PER_SEC * 60) );
  return EXIT_SUCCESS;
}
