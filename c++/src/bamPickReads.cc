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

map<string, int> read_ids(const string& fi) {
  ifstream fhi(fi.c_str());
  if(!fhi.is_open()) {
    cout << format("cannot open file: %s") % fi;
    exit(1);
  }
  
  map<string, int> ids;
  string line;
  vector<string> ss;
  while(fhi.good()) {
    getline(fhi, line);
    if(line.length() == 0) continue;
    boost::split(ss, line, boost::is_any_of("\t"));
    ids.insert( pair<string, int>(ss[0], boost::lexical_cast<int>(ss[2])) );
  }
  cout << format("%d ids read from %s\n") % ids.size() % fi;
  return ids;
}

int main(int argc, char *argv[]) {
  string f_bam, f_id, f_out, tag;
  clock_t time1 = clock();
  po::options_description cmdOpts("Allowed options");
  cmdOpts.add_options()
    ("help,h", "produce help message")
    ("f_bam,b", po::value<string>(&f_bam), "input BAM file (sorted by read name)")
    ("f_id,i", po::value<string>(&f_id), "read name file")
    ("f_out,o", po::value<string>(&f_out), "output file")
  ;
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, cmdOpts), vm);
  po::notify(vm);
  if(vm.count("help") || !vm.count("f_id") || !vm.count("f_bam") || !vm.count("f_out")) {
    cout << cmdOpts << endl;
    return 1;
  }
  
  Utilities utils;
  BamReader reader;
  BamWriter writer;
  reader.Open(f_bam); 
  RefVector refs = reader.GetReferenceData();

  map<string, int> sim= read_ids(f_id);
  map<string, int>::iterator sim_it;
  pair< map<string, int>::iterator, bool > sim_p;

  BamAlignment al1, al2;
  string rg, id;
  ofstream fho( f_out.c_str() );
  fho << "id\ttype1\tchr1\tbeg1\tstrand1\tcigar1\tseq1\ttype2\tchr2\tbeg2\tstrand2\tcigar2\tseq2\n";
  int refId = reader.GetReferenceID("chr5");
  while( reader.GetNextAlignment(al1) ) {
    reader.GetNextAlignment(al2);
    if(al1.Name != al2.Name) {
      cout << format("not a pair of reads: %s != %s\n") % al1.Name % al2.Name;
      utils.PrintRead(al1);
      utils.PrintRead(al2);
      exit(0);
    }
    id = al1.Name;
    al1.GetTag("RG", rg);
    sim_it = sim.find(id);
    if(sim_it != sim.end()) {
      if(al2.IsFirstMate() && !al1.IsFirstMate()) swap(al1, al2);
      if(al2.IsFirstMate() || !al1.IsFirstMate()) {
        cout << format("%s: not first + second\n") % id;
        exit(0);
      }
      int beg1 = al1.Position, beg2 = al2.Position;
      bool strand1 = al1.IsReverseStrand(), strand2 = al2.IsReverseStrand();
      string chr1 = refs[al1.RefID].RefName, chr2 = refs[al2.RefID].RefName;
      vector<CigarOp> cigars1 = al1.CigarData, cigars2 = al2.CigarData;
      string cigar1 = utils.getCigarString(cigars1);
      string cigar2 = utils.getCigarString(cigars2);
      string seq1 = al1.QueryBases, seq2 = al2.QueryBases;
//  type: 1(mapped), 2(soft-clipped), 3(unmapped)
      int type1 = (!al1.IsMapped() || al1.RefID != refId) ? 3 :
          ((cigars1[0].Type == 'S' && cigars1[0].Length >= 10) ||
           (cigars1[cigars1.size()-1].Type == 'S' && cigars1[cigars1.size()-1].Length >= 10)) ? 2 : 1;
      int type2 = (!al2.IsMapped() || al2.RefID != refId) ? 3 :
          ((cigars2[0].Type == 'S' && cigars2[0].Length >= 10) ||
           (cigars2[cigars2.size()-1].Type == 'S' && cigars2[cigars2.size()-1].Length >= 10)) ? 2 : 1;
      fho << format("%d\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%s\t%d\t%d\t%s\t%s\n") 
        % id % type1 % chr1 % beg1 % strand1 % cigar1 % seq1 % type2 % chr2 % beg2 % strand2 % cigar2 % seq2;
     
/*      int orphanIsFirst = sim_it->second;
      if(orphanIsFirst) swap(al1, al2);
      int firstMate = al2.IsFirstMate() ? 1 : 2;
      string mateStrand = al2.IsMateReverseStrand() ? "-" : "+";
      fho << format("@%s/%d\n%s\n") % al2.Name % firstMate % al2.QueryBases;
      int mateEnd = al1.GetEndPosition(false, false);
      fho << format("%s\t%s\t%d\t%d\t%d\t%s\n") % mateStrand % "chr5" % mateEnd % al1.MapQuality % 270 % tag */;
    }
  }
  reader.Close();
  cerr << right << setw(60) << format("(running time: %.01f minutes)\n") % 
    ( (double)(clock() - time1) / ((double)CLOCKS_PER_SEC * 60) );
  return EXIT_SUCCESS;
}
