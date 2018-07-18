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

struct ISD {
  string id;
  int run;
  double mean, sd, median, rsd, mad, mld, mno, predicted;
  bool tag_hydra;
};
typedef map<string, ISD> ISDMap;
const int getEditDist(const string& tag_md) {
  int mm = 0;
  string::size_type word_pos(0);
  while(word_pos != string::npos) {
    word_pos = tag_md.find("X", word_pos);
    if(word_pos != string::npos) {
      mm ++;
      word_pos += 1;
    }
  }
  return mm;
}
const bool read_isdmap(const string& fi, ISDMap& im, set<string>& ids) {
  ifstream fhi(fi.c_str());
  if(!fhi.is_open()) {
    cout << format("cannot open file: %s") % fi;
    return false;
  }
  string line;
  if(fhi.good()) getline(fhi, line);
  while(fhi.good()) {
    getline(fhi, line);
    if(line.length() == 0) continue;
    vector<string> ss;
    boost::split(ss, line, boost::is_any_of("\t"));
    ISD isd;
    isd.id = ss[0];
    isd.run = boost::lexical_cast<int>(ss[1]);
    isd.tag_hydra = ss[10] == "1" ? true : false;
    if(!isd.tag_hydra) continue;

    isd.mean = boost::lexical_cast<double>(ss[2]);
    isd.sd = boost::lexical_cast<double>(ss[3]);
    isd.median = boost::lexical_cast<double>(ss[4]);
    isd.rsd = boost::lexical_cast<double>(ss[5]);
    isd.mad = boost::lexical_cast<double>(ss[6]);
    isd.mld = boost::lexical_cast<double>(ss[7]);
    isd.mno = boost::lexical_cast<double>(ss[8]);
    isd.predicted = boost::lexical_cast<double>(ss[9]);
    
    ids.insert(isd.id);
    string rg = str(format("%s_%02d") % isd.id % isd.run);
    im.insert( pair<string, ISD> (rg, isd) );
  }
  return true;
}

int main(int argc, char *argv[]) {
  string fi, fo, regionStr;
  double is_mno;
  vector<string> rgs_vec;
  clock_t time1 = clock();
  po::options_description cmdOpts("Allowed options");
  cmdOpts.add_options()
    ("help,h", "produce help message")
    ("fi,i", po::value<string>(&fi), "input (BAM) file")
    ("fo,o", po::value<string>(&fo), "output (BED) file")
    ("is_mno,m", po::value<double>(&is_mno), "max insert size to be considered normal")
    ("rgs,g", po::value< vector<string> >(&rgs_vec), "read group(s) to include")
    ("region,r", po::value<string>(&regionStr), "region")
  ;
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, cmdOpts), vm);
  po::notify(vm);
  if(vm.count("help") || !vm.count("fi") || !vm.count("fo") || 
    !vm.count("is_mno") || !vm.count("rgs")) {
    cout << cmdOpts << endl;
    return 1;
  }

  Utilities util;
  BamReader reader;
  if( reader.Open(fi) ) {
    cout << format("%s opened ...\n") % fi;
  } else {
    cout << format("cannot read from %s ...\n") % fi;
    exit(1);
  }
  if(regionStr != "") {
    BamRegion region;
    if( util.ParseRegionString(regionStr, reader, region) ) {
      reader.SetRegion(region);
      cout << format("region set to: %s\n") % regionStr;
    } else {
      cerr << format("invalid region string: %s\n") % regionStr;
      exit(1);
    }
  } else {
    cout << "region not specified: working on entire BAM\n";
  }
 
  RefVector refs = reader.GetReferenceData();
  
  set<string> rgs;
  set<string>::iterator rgs_it;
  for(vector<string>::iterator it = rgs_vec.begin(); it != rgs_vec.end(); it ++)
    rgs.insert( *it );

  ofstream fho(fo.c_str());
  map<string, int> readh;
  map<string, int>::iterator it;
  pair<map<string, int>::iterator, bool> p;
  BamAlignment al;
  string rg;
  int is, mm;
  while( reader.GetNextAlignment(al) ) {
    if(!al.IsMapped() || !al.IsMateMapped() || al.RefID != al.MateRefID 
      || al.IsDuplicate() || al.MapQuality==0) continue;
    string strand1 = al.IsReverseStrand() ? "-" : "+";
    string strand2 = al.IsMateReverseStrand() ? "-" : "+";
    if(strand1 == strand2) continue;
    is = abs(al.InsertSize);
    al.GetTag("RG", rg);
    al.GetTag("NM", mm);
    
    rgs_it = rgs.find(rg);
    if(rgs_it != rgs.end() && is >= is_mno) {
      p = readh.insert( pair<string, int> (al.Name, mm) );
      if(p.second == false) {
        int mm2 = p.first->second;
        string chr = refs[al.RefID].RefName;
        fho << format("%s\t%d\t%s\t%d\t%d\t%d\t%s\t%d\t%s\t%s\t%d\t%d\n")
          % chr % al.Position % (al.Position+al.Length)
          % chr % al.MatePosition % (al.MatePosition + al.Length)
          % al.Name % 1 % strand1.c_str() % strand2.c_str() % mm % mm2;
        readh.erase(p.first);
      }
    }
  }
  reader.Close();

  cerr << right << setw(60) << format("(running time: %.01f minutes)\n") % 
    ( (double)(clock() - time1) / ((double)CLOCKS_PER_SEC * 60) );
  return EXIT_SUCCESS;
}
