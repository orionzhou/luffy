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

typedef map< string, string > StrStrMap;
StrStrMap get_readgroup(const string& fi, const string& sm) {
  ifstream fhi(fi.c_str());
  if(!fhi.is_open()) {
    cout << format("cannot open file: %s") % fi;
    exit(1);
  }

  StrStrMap rm;
  StrStrMap::iterator it;
  pair< StrStrMap::iterator, bool > ret;
  
  string line;
  vector<string> ss;
  getline(fhi, line);
  while(fhi.good()) {
    getline(fhi, line);
    if(line.length() == 0) continue;
    boost::split(ss, line, boost::is_any_of("\t"));
    if(ss[0] != sm) continue;
    string rgn = ss[1];
    string rgo = ss[7];
    ret = rm.insert( pair<string, string> (rgo, rgn) );
    if(ret.second == false) {
      cerr << format("read group [%s] appeared >1 times\n") % rgo;
      exit(1);
    }
    cout << format("\t%s  ->  %s\n") % rgo % rgn;
  }
  return rm;
}
int main(int argc, char *argv[]) {
  string fi, dirO, sample, fm;
  clock_t time1 = clock();
  po::options_description cmdOpts("Allowed options");
  cmdOpts.add_options()
    ("help,h", "produce help message")
    ("in,i", po::value<string>(&fi), "input (BAM) file")
    ("out,o", po::value<string>(&dirO), "output directory")
    ("sample,s", po::value<string>(&sample), "sample id")
    ("mapping,m", po::value<string>(&fm), "@RG mapping file")
  ;

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, cmdOpts), vm);
  po::notify(vm);
  if(vm.count("help") || !vm.count("in") || !vm.count("out") || !vm.count("sample")) {
    cout << cmdOpts << endl;
    exit(1);
  }
  
//  try { fs::create_directory(dirO); }
//  catch (fs::filesystem_error &e) { cerr << e.what() << endl; }

  Utilities util;
  BamReader2 reader;
  if( reader.Open(fi) ) {
    cout << format("header read from %s ...\n") % fi;
  } else {
    cout << format("cannot read from %s ...\n") % fi;
    exit(1);
  }

  StrStrMap rm = get_readgroup(fm, sample);
  StrStrMap::iterator rm_it;
  
  SamReadGroupDictionary rgd = reader.GetHeader().ReadGroups;
  map<string, ofstream*> fm1;
  map<string, ofstream*> fm2;
  map<string, ofstream*>::iterator fm_it;
  for(SamReadGroupConstIterator it=rgd.ConstBegin(); it!=rgd.ConstEnd(); ++it) {
    string rgo = it->ID;
    rm_it = rm.find(rgo);
    if(rm_it == rm.end()) {
      cerr << format("read group [%s] not found\n") % rgo;
      exit(1);
    }
    string rgn = rm_it->second;
//    cout << format("\t%s  ->  %s\n") % rgo % rgn;
    string fo1_str = str( format("%s/%s.1.fq") % dirO % rgn );
    string fo2_str = str( format("%s/%s.2.fq") % dirO % rgn );
    ofstream* fho1 = new ofstream;
    ofstream* fho2 = new ofstream;
    fm1[rgo] = fho1;
    fm1.find(rgo)->second->open( fo1_str.c_str() );
    fm2[rgo] = fho2;
    fm2.find(rgo)->second->open( fo2_str.c_str() );
  }

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
      exit(1);
    }
    BamAlignment al1, al2;
    if(alnpairs[0].second.IsFirstMate()) {
      al1 = alnpairs[0].second;
      al2 = alnpairs[0].first;
    } else {
      al1 = alnpairs[0].first;
      al2 = alnpairs[0].second;
    }
    if(!al1.IsFirstMate() || al2.IsFirstMate()) {
      cout << "alns both first/second\n";
      util.PrintReads(name, n_total, alnpairs);
      exit(1);
    }

    string rgo;
    al1.GetTag("RG", rgo);
    if(rm.find(rgo) == rm.end()) {
      cout << format("read group[%s] not found\n") % rgo;
      exit(1);
    } else {
      ofstream* fho1 = fm1.find(rgo)->second;
      ofstream* fho2 = fm2.find(rgo)->second;
      string al1_b(al1.QueryBases), al1_q(al1.Qualities), al2_b(al2.QueryBases), al2_q(al2.Qualities);
      if(al1.IsReverseStrand()) {
        util.ReverseComplement(al1_b);
        util.Reverse(al1_q);
      }
      (*fho1) << format("@%s/1\n%s\n+\n%s\n") % name % al1_b % al1_q;
      if(al2.IsReverseStrand()) {
        util.ReverseComplement(al2_b);
        util.Reverse(al2_q);
      }
      (*fho2) << format("@%s/2\n%s\n+\n%s\n") % name % al2_b % al2_q;
    }
  }
  reader.Close();
  for(fm_it = fm1.begin(); fm_it != fm1.end(); fm_it ++) fm_it->second->close();
  for(fm_it = fm2.begin(); fm_it != fm2.end(); fm_it ++) fm_it->second->close();
  
  for(rm_it = rm.begin(); rm_it != rm.end(); rm_it ++) {
    string rgn = rm_it->second;
    string cmd1 = str( format("gzip -f %s/%s.1.fq") % dirO % rgn );
    string cmd2 = str( format("gzip -f %s/%s.2.fq") % dirO % rgn );
    system(cmd1.c_str());
    system(cmd2.c_str());
  }

  cout << right << setw(60) << format("(running time: %.01f minutes)\n") % ( (double)(clock() - time1) / ((double)CLOCKS_PER_SEC * 60) );
  return EXIT_SUCCESS;
}

