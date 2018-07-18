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
#include <geco.h>
using namespace std;
using namespace geco;
using boost::format;
namespace fs = boost::filesystem;
namespace po = boost::program_options;

template <class T> 
ostream & operator << (ostream & out, const gArray<T>& array) {
  size_t lfield;
  if(typeid(T) == typeid(gScore) ){
    lfield=7;
    cout.precision(1);
    cout << fixed;
  } else {
    lfield=5;
  }
  string sep(4+lfield*array.getSize(),'-');
  out << "Pos:";
  for (gPos i=0;i<array.getSize();i++)
    out << setw(lfield) << i;
  out << endl << sep.c_str() << endl;
  out << "Val:";
  for (gPos i=0;i<array.getSize();i++) {
    if(array.isNA(i)) out << setw(lfield) << "*";
    else out << setw(lfield) << array[i];
  }
  out << endl;
  return out;
}
bool get_geco_transcript(const gLocalTxtSequenceRetriever& ret, const gSequenceRetriever& seqret, Transcript& t, gTranscript& gtr) {
  gElementSequenceMode mode(gRef);
  gPos cdsBeg(t.cdsBeg), cdsEnd(cdsEnd);
  gArray<gPos> exonBegs(0, t.locG_exon.size(), false, false);
  gArray<gPos> exonEnds(0, t.locG_exon.size(), false, false);
 
  string locStr_exon = "";
  for(int i = 0; i < int(t.locG_exon.size()); i ++) {
    LocPair lp = t.locG_exon[i];
    locStr_exon += str( format(" %d-%d") % lp.first % lp.second );
    exonBegs.setValue( i, lp.first - 1 );
    exonEnds.setValue( i, lp.second );
  }
  gTranscript gtr2(t.chr, exonBegs, exonEnds, t.forward, ret, mode, 0, 0, cdsBeg, cdsEnd);
  gtr = gtr2;

  LocPairVec locL_intron;
  for(int i = 0; i < int(gtr.getIntronCount()); i ++) {
    gElementInterval bd = gtr.getIntronInterval(i);
    locL_intron.push_back( LocPair (bd.getStart()+1, bd.getEnd()) );
  }
  t.locL_intron = locL_intron;

  t.locLs["cds"] =  t.locL_cds;
  t.locLs["intron"] = t.locL_intron;
  t.locLs["utr5"] = t.locL_utr5;
  t.locLs["utr3"] = t.locL_utr3;

  t.seq = seqret.getSequence(t.chr, t.beg-1, t.end);
  boost::to_upper(t.seq);
  return true;
}
string locVec2String(const LocPairVec& loc) {
  vector<string> strVec;
  for(LocPairVec::const_iterator it = loc.begin(); it != loc.end(); it ++) {
    string locStr = str( format("%d-%d") % it->first % it->second );
    strVec.push_back( locStr );
  }
  return boost::algorithm::join(strVec, ",");
}
int binary_search(vector<int>& vec, const int& key) {
  if(adjacent_find(vec.begin(), vec.end(), greater<int>()) != vec.end()) {
    cout << format("vector not sorted!\n");
    exit(1);
  }
  int mid, lb=0, ub=vec.size()-1, flag=1;
  for(mid = (lb+ub)/2; lb <= ub; mid=(lb+ub)/2) {
    if(vec[mid] == key) {
      flag = 0;
      break;
    } else if(vec[mid] > key) {
      ub = mid - 1;
    } else {
      lb = mid + 1;
    }
  }
  if(vec[mid] > key) mid --;
  return mid;
}
vector<string> get_region_from_pos(const Transcript& t, const gArray<gPos>& poss) {
  map<int, string> typeMap; 
  vector<int> begs;
  map<string, LocPairVec>::const_iterator it;
  LocPairVec::const_iterator it2;
  for(it = t.locLs.begin(); it != t.locLs.end(); it ++) {
    string type = it->first;
    LocPairVec lpVec = it->second;
    for(it2 = lpVec.begin(); it2 != lpVec.end(); it2++) {
      int beg = it2->first;
      begs.push_back( beg );
      typeMap.insert( pair<int, string> (beg, type) );
    }
  }
  sort(begs.begin(), begs.end());

  vector<string> types;
  for(gSize i = 0; i < poss.getSize(); i ++) {
    int idx = binary_search(begs, poss[i]+1);
    types.push_back( typeMap[begs[idx]] );
  }
  return types;
}
string get_subseq_from_loc(const string& seq_ori, const LocPairVec& loc, bool forward) {
  vector<string> subseqVec;
  string seq = seq_ori;
  if( !forward ) {
    gSequence gSeq(seq_ori);
    seq = gSeq.getReverseComplement();
  }

  for(LocPairVec::const_iterator it = loc.begin(); it < loc.end(); it ++) 
    subseqVec.push_back( seq.substr(it->first-1, it->second-it->first+1) );
  return boost::algorithm::join(subseqVec, "");
}


