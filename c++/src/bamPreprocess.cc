#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <vector>
#include <map>
#include <time.h>
#include <boost/format.hpp>
#include <boost/program_options.hpp>
#include <boost/regex.hpp>
#include <boost/algorithm/string.hpp>
#include <api/BamReader.h>
#include <api/BamWriter.h>
#include "read.h" 
#include "bam_utils.h" 
using namespace std;
using boost::format;
namespace fs = boost::filesystem;
namespace po = boost::program_options;
using namespace BamTools;

bool changeHeader(const SamHeader& header1, SamHeader& header2) {
  if(header1.HasVersion()) header2.Version = header1.Version;
  if(header1.HasSortOrder()) header2.SortOrder = header1.SortOrder;
  if(header1.HasGroupOrder()) header2.GroupOrder = header1.GroupOrder;
  if(header1.HasReadGroups()) header2.ReadGroups = header1.ReadGroups;
  if(header1.HasPrograms()) header2.Programs = header1.Programs;
  if(header1.HasComments()) header2.Comments = header1.Comments;
  return true;
}
bool changeRefVector(RefVector refs1, map<string, LocMap> refMap, vector<LocMap>& refInfo, RefVector& refs2) {
  RefVector::iterator refIter = refs1.begin();
  uint32_t cnt = 0;
  boost::regex e ("\\.[0-9A-Z]+", boost::regex_constants::icase|boost::regex_constants::perl);
  map<string, uint32_t> refIds;
  for (; refIter != refs1.end(); refIter++) {
    RefData ref1 = *refIter;
    string name1 = ref1.RefName;
    int id1 = distance(refs1.begin(), refIter);
    LocMap locMap;
    locMap.id1 = id1;
    locMap.name1 = name1;
    if(name1 == "chl_Mt") {
      locMap.name2 = "NC_003119";
      locMap.start2 = 1;
      locMap.end2 = ref1.RefLength;
    } else if(name1.find("Mt3.5") != string::npos) {
      locMap.name2 = "chr" + name1.substr(name1.length()-1, 1);
      locMap.start2 = 1;
      locMap.end2 = ref1.RefLength;
    } else {
      name1 = boost::regex_replace(name1, e, "");
      locMap.name1 = name1;
      map<string, LocMap>::iterator it = refMap.find(name1);
      if(it == refMap.end()) {
        cout << format("do not have info for %s from refMap\n") % name1;
        return false;
      } else {
        locMap.name2 = it->second.name2;
        locMap.start2 = it->second.start2;
        locMap.end2 = it->second.end2;
      }
    }
    pair< map<string, uint32_t>::iterator, bool > p;
    p = refIds.insert( pair<string, uint32_t> (locMap.name2, cnt) );
    if(p.second == false) {
      uint32_t idx = p.first->second;
      locMap.id2 = idx;
      uint32_t refLen1 = refs2[idx].RefLength;
      if(locMap.end2 > refLen1) refs2[idx].RefLength = locMap.end2;
    } else {
      locMap.id2 = cnt;
      RefData ref2;
      ref2.RefName = locMap.name2;
      ref2.RefLength = locMap.end2;
      refs2.push_back(ref2);
      cnt ++;
    }
    refInfo.push_back(locMap);
  }
  return true;
}
bool changeAlnLoc(const vector<LocMap>& refInfo, BamAlignment& al) {
  if(al.RefID >= (int)refInfo.size() || al.MateRefID >= (int)refInfo.size()) {
    cout << format("invalid RefID[%d] or MateRefID[%d]\n") % al.RefID % al.MateRefID;
    return false;
  }
  if(al.RefID != -1) {
    LocMap loc1 = refInfo[al.RefID];
    al.RefID = loc1.id2;
    al.Position += loc1.start2 - 1;
  }
  if(al.MateRefID != -1) {
    LocMap loc2 = refInfo[al.MateRefID];
    al.MateRefID = loc2.id2;
    al.MatePosition += loc2.start2 - 1;
  }
  return true;
}

int main(int argc, char *argv[]) {
  string fi, fo, fconfig;
  clock_t time1 = clock();
  po::options_description cmdOpts("Allowed options");
  cmdOpts.add_options()
    ("help,h", "produce help message")
    ("in,i", po::value<string>(&fi), "input file")
    ("out,o", po::value<string>(&fo), "output file prefix")
    ("config,c", po::value<string>(&fconfig)->default_value("/project/youngn/zhoup/Data/repo/mt_35/01_reference/11_ref_mapping.tbl"), "config file")
  ;
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, cmdOpts), vm);
  po::notify(vm);
  if(vm.count("help") || !vm.count("in") || !vm.count("out")) {
    cout << cmdOpts << endl;
    exit(1);
  }

  map<string, LocMap> refMap;
  if( !read_refmap(fconfig, refMap) ) {
    cout << format("error building RefMap from %s\n") % fconfig;
    exit(0);
  } else {
    cout << format("RefMap built from %s\n") % fconfig;
  }

  Utilities util;
  BamReader2 reader;
  BamWriter writer;
  fs::path fo1( fo + ".bam" );
  fs::path fo2( fo + "_stat.txt" );
  fs::path fo3( fo + "_insert.txt" );
  if( reader.Open(fi) ) {
    cout << format("header read from %s ...\n") % fi;
  } else {
    cout << format("cannot read from %s ...\n") % fi;
    exit(1);
  }

  SamHeader header;
  if( !changeHeader(reader.GetHeader(), header) ) {
    cout << "cannot change BamHeader" << endl;
    exit(0);
  } else {
    cout << "BamHeader extracted " << endl;
  }
  RefVector refs;
  vector<LocMap> refInfo;
  if( !changeRefVector(reader.GetReferenceData(), refMap, refInfo, refs) ) {
    cout << format("error generating new RefVector\n");
    exit(0);
  } else {
    cout << format("%d refSeqs changed to %d\n") % (int)refInfo.size() % (int)refs.size();
  }
  writer.Open( fo1.string(), header, refs );

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
    for (int i = 0; i < n_total; i ++) {
      BamAlignment al1(alnpairs[i].first), al2(alnpairs[i].second);
      if(type != 3 && n_total == 1) {
        if(al1.IsMapped()) al1.MapQuality = 40;
        if(al2.IsMapped()) al2.MapQuality = 40;
      }
      if(type == 1 || type == 2) {
        al1.SetIsProperPair(true);
        al2.SetIsProperPair(true);
      }
      al1.SetIsPrimaryAlignment(true);
      al2.SetIsPrimaryAlignment(true);
      al1.AddTag("PZ", "i", i+1);
      al2.AddTag("PZ", "i", i+1);
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
      if(type != 3) { 
        if(!changeAlnLoc(refInfo, al1) || !changeAlnLoc(refInfo, al2)) {
          cout << format("error changing location for %s") % name;
          exit(0);
        }
      }
      al1.RemoveTag("NH");
      al1.RemoveTag("IH");
      al1.RemoveTag("YU");
      al2.RemoveTag("NH");
      al2.RemoveTag("IH");
      al2.RemoveTag("YU");
      al1.AddTag("X0", "i", n_total);
      al2.AddTag("X0", "i", n_total);
      writer.SaveAlignment(al1);
      writer.SaveAlignment(al2);
    }
  }
  reader.Close();
  writer.Close();

  ofstream fho2( fo2.string().c_str() );
  fho2 << format("%d\t%d\t%d\t%d\t%d\t%d\t%d\n") % cnt_all % cnt_unmapped 
    % cnt_orphan % cnt_mapped % cnt_orphan_u % cnt_mapped_u % cnt_mapped_up;
  fho2.close();

  ofstream fho3( fo3.string().c_str() );
  for( it1=m1.begin(); it1!=m1.end(); it1++ ) 
    fho3 << format("%d\t%d\n") % it1->first % it1->second;
  fho3.close();

  cout << right << setw(60) << format("(running time: %.01f minutes)\n") % ( (double)(clock() - time1) / ((double)CLOCKS_PER_SEC * 60) );
  return EXIT_SUCCESS;
}

