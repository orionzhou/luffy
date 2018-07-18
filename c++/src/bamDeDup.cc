#include <iostream>
#include <iomanip>
#include <vector>
#include <set>
#include <time.h>
#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <boost/program_options.hpp>
#include <api/BamReader.h>
#include <api/BamWriter.h>
using namespace std;
using boost::format;
namespace fs = boost::filesystem;
namespace po = boost::program_options;
using namespace BamTools;

int main(int argc, char *argv[]) {
  string fi, fo;
  clock_t time1 = clock();
  po::options_description cmdOpts("Allowed options");
  cmdOpts.add_options()
    ("help,h", "produce help message")
    ("in,i", po::value<string>(&fi), "input (BAM) file")
    ("out,o", po::value<string>(&fo), "output (BAM) file")
  ;

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, cmdOpts), vm);
  po::notify(vm);
  if(vm.count("help") || !vm.count("in") || !vm.count("out")) {
    cout << cmdOpts << endl;
    return 1;
  }

  BamReader reader;
  fs::path fi1( fi );
  fs::path fi2( fi + ".bai" );
  if( exists(fi1) ) {
    reader.Open(fi1.string());
    reader.OpenIndex(fi2.string());
    if( reader.HasIndex() ) {
      cout << format("%s opened with index ...\n") % fi;
    } else {
      cout << format("%s opened without index ...\n") % fi;
    }
  } else {
    cerr << format("cannot read from %s ...\n") % fi;
    exit(1);
  }

  BamWriter writer;
  writer.SetCompressionMode(BamWriter::Uncompressed);
  writer.Open(fo, reader.GetHeader(), reader.GetReferenceData());

  map<string, BamAlignment> als;
  map<string, BamAlignment>::iterator it;
  BamAlignment al, al2;
  uint32_t cnt1(0), cnt1d(0), cnt2(0), cnt2d(0), cnt3(0);
  while( reader.GetNextAlignment(al) ) {
    if( !al.IsPrimaryAlignment() ) {
      if( !al.IsDuplicate() ) {
        writer.SaveAlignment(al);
      }
    } else {
      if( al.IsMapped() && al.IsMateMapped() ) {
        if( al.IsFirstMate() ) cnt1 ++;
        if( al.IsDuplicate() ) {
          if( al.IsFirstMate() ) cnt1d ++;
        } else {
          writer.SaveAlignment(al);
        }
      } else if ( !al.IsMapped() && !al.IsMateMapped() ) {
        if( al.IsFirstMate() ) cnt3 ++;
        writer.SaveAlignment(al);
      } else {
        if( al.IsFirstMate() ) cnt2 ++;
        it = als.find(al.Name);
        if(it == als.end()) {
          als.insert( pair<string, BamAlignment> (al.Name, al) );
        } else {
          al2 = it->second;
          if( al.IsDuplicate() || al2.IsDuplicate() ) {
            cnt2d ++;
          } else {
            writer.SaveAlignment(al2);
            writer.SaveAlignment(al);
          }
          als.erase( it );
        }
      }
    }
  }
  reader.Close();
  writer.Close();
  
  if( als.size() > 0 ) 
    cerr << "BAM error: " << als.begin()->first << endl;
  cerr << format("%10d paired   / %10d duplicates\n") % cnt1 % cnt1d;
  cerr << format("%10d unpaired / %10d duplicates\n") % cnt2 % cnt2d;
  cerr << format("%10d unmapped\n") % cnt3;
  cerr << right << setw(60) << format("(running time: %.01f minutes)\n") % 
    ( (double)(clock() - time1) / ((double)CLOCKS_PER_SEC * 60) );
  return EXIT_SUCCESS;
}


 
