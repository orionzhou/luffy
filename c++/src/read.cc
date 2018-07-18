#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <vector>
#include <map>
#include <time.h>
#include <fstream>
#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <boost/program_options.hpp>
#include <boost/regex.hpp>
#include <boost/assign/std/vector.hpp>
#include <boost/algorithm/string.hpp>
#include "read.h" 
using boost::format;
using namespace std;
using namespace boost::assign;
namespace fs = boost::filesystem;
namespace po = boost::program_options;

bool read_refmap(const string& fi, map<string, LocMap>& refMap) {
  ifstream fhi(fi.c_str());
  if(!fhi.is_open()) {
    cout << format("cannot open file: %s") % fi;
    return false;
  }
  string line;
  pair< map<string, LocMap>::iterator, bool > p;
  while(fhi.good()) {
    getline(fhi, line);
    if(line.length() == 0) continue;
    vector<string> ss;
    boost::split(ss, line, boost::is_any_of("\t"));
    LocMap locMap;
    locMap.name1 = ss[0];
    locMap.name2 = ss[1];
    locMap.start2 = atoi(ss[2].c_str());
    locMap.end2 = atoi(ss[3].c_str());
    p = refMap.insert( pair<string, LocMap> (ss[0], locMap) );
    if(p.second == false) {
      cout << format(">2 entries [%s] in locMap file\n") % ss[0];
      return false;
    }
  }
  return true;
}
bool compare_loc_pair(LocPair locA, LocPair locB) { return (locA.first < locB.first); }
LocPairVec get_loc_from_string(const string& locStr, const boost::regex expression) {
  LocPairVec lpVec;
  boost::match_flag_type flags = boost::match_default;
  boost::match_results<string::const_iterator> what;
  string::const_iterator beg = locStr.begin(), end = locStr.end();
  while( regex_search(beg, end, what, expression, flags) ) {
    int loc_beg = atoi( string(what[1].first, what[1].second).c_str() );
    int loc_end = atoi( string(what[2].first, what[2].second).c_str() );
    if(loc_beg > loc_end) swap(loc_beg, loc_end);
    lpVec.push_back( pair<int, int> (loc_beg, loc_end) );
    beg = what[0].second;
  }
  sort(lpVec.begin(), lpVec.end(), compare_loc_pair);
  return lpVec;
}
LocPairVec loc_global_2_local(const uint32_t& beg, const uint32_t& end, const bool& forward, const LocPairVec& locG) {
  LocPairVec locL;
  if(forward) {
    for(LocPairVec::const_iterator it = locG.begin(); it != locG.end(); it ++)
      locL.push_back( LocPair (it->first - beg + 1, it->second - beg + 1) );
  } else {
    for(LocPairVec::const_reverse_iterator it = locG.rbegin(); it != locG.rend(); it ++)
      locL.push_back( LocPair (end - it->second + 1, end - it->first + 1) );
  }
  return locL;
}
bool read_gtb(const string& f_gtb, map <string, Transcript>& tm) {
  ifstream fhi(f_gtb.c_str());
  if(!fhi.is_open()) { cout << format("cannot open file: %s") % f_gtb; return false; }

  string line;
  pair< map<string, Transcript>::iterator, bool > p;
  if(fhi.good()) getline(fhi, line);
  boost::regex expression("([0-9]+)\\.\\.([0-9]+)");
  
  while(fhi.good()) {
    getline(fhi, line);
    if(line.length() == 0) continue;
    Transcript t;
    vector<string> ss;
    boost::split(ss, line, boost::is_any_of("\t"));
    t.id = ss[1];
    t.id_gene = ss[0];
    t.chr = ss[2];
    t.forward = ss[3] == "-1" ? false : true;
    t.frame = atoi(ss[4].c_str());
    t.source = ss[9];
    t.type = ss[10];
    t.conf = ss[11];
    t.note = ss[12];

    t.locG_cds = get_loc_from_string(ss[5], expression);
    t.locG_exon = get_loc_from_string(ss[6], expression);
    t.locG_utr5 = get_loc_from_string(ss[7], expression);
    t.locG_utr3 = get_loc_from_string(ss[8], expression);

    t.beg = t.locG_exon.at(0).first;
    t.end = t.locG_exon.at(t.locG_exon.size()-1).second;
    t.cdsBeg = t.locG_cds.at(0).first;
    t.cdsEnd = t.locG_cds.at(t.locG_cds.size()-1).second;
    
    t.locL_cds = loc_global_2_local(t.beg, t.end, t.forward, t.locG_cds);
    t.locL_exon = loc_global_2_local(t.beg, t.end, t.forward, t.locG_exon);
    t.locL_utr5 = loc_global_2_local(t.beg, t.end, t.forward, t.locG_utr5);
    t.locL_utr3 = loc_global_2_local(t.beg, t.end, t.forward, t.locG_utr3);

    p = tm.insert( pair<string, Transcript> (t.id, t) );
    if(p.second == false) {
      cout << format(">2 trnascripts named [%s] gtb file\n") % t.id;
      return false;
    }
  }
  cout << format("%d transcript models read from %s\n") % tm.size() % f_gtb;
  return true;
}
set<string> read_ids(const string& fi) {
  ifstream fhi(fi.c_str());
  if(!fhi.is_open()) {
    cout << format("cannot open file: %s") % fi;
    exit(1);
  }
  set<string> ids;
  string line;
  uint32_t i = 0;
  while(fhi.good()) {
    getline(fhi, line);
    if(++i == 1 && line == "id") continue;
    if(line.length() == 0) continue;
    ids.insert(line);
  }
  cout << format("%d ids read from %s\n") % ids.size() % fi;
  return ids;
}
vector<Loc> read_locs(const string& fi) {
  ifstream fhi(fi.c_str());
  if(!fhi.is_open()) {
    cout << format("cannot open file: %s") % fi;
    exit(1);
  }
  vector<Loc> lv;
  string line;
  vector<string> ss;
  
  getline(fhi, line);
  boost::split(ss, line, boost::is_any_of("\t"));
  int idx_c, idx_b, idx_e, idx_i = -1;
  for(uint32_t i = 0; i < ss.size(); i ++) {
    if(ss[i] == "id") idx_i = i;
    if(ss[i] == "chr") idx_c = i;
    if(ss[i] == "beg") idx_b = i;
    if(ss[i] == "end") idx_e = i;
  }

  while(fhi.good()) {
    getline(fhi, line);
    if(line.length() == 0) continue;
    boost::split(ss, line, boost::is_any_of("\t"));
    Loc loc;
    loc.chr = ss[idx_c];
    loc.beg = atoi(ss[idx_b].c_str());
    loc.end = atoi(ss[idx_e].c_str());
    loc.length = loc.end - loc.beg + 1;
    if(idx_i != -1) loc.id = ss[idx_i];
    lv.push_back(loc);
  }
  return lv;
}
IntVecVec get_raw_cov(const string& chr, const int& beg, const int& end, const string& f_cov, const int& opt_conf) {
  string cmd = str( format("tabix %s %s:%d-%d") % f_cov % chr % beg % end );
    
  int ncol, nrow = end - beg + 1, col_offset = 3;
  if(opt_conf == 1) {
    ncol = 84;
  } else if(opt_conf == 2) {
    ncol = 288;
  } else {
    cout << "unknown opt_conf: " << opt_conf << endl;
    exit(1);
  }
  
  int MAX_BUFSIZE = 25535;
  char buffer[MAX_BUFSIZE];
  FILE *stream = popen(cmd.c_str(), "r");
  vector<uint32_t> res;
  int row = 0;
  while( fgets(buffer, MAX_BUFSIZE, stream) != NULL ) {
    string line(buffer);
    if(line.substr(line.length()-1, 1) == "\n") line.erase(line.length()-1, 1);
    if(line.length() == 0) continue;
    vector<string> ss;
    boost::split(ss, line, boost::is_any_of("\t"));
    if(int(ss.size()) != ncol + col_offset) {
      cout << format("not %d+%d columns: %d\n%s\n") % ncol % col_offset % ss.size() % line;
      exit(1);
    }
    for(uint32_t j = col_offset; j < ss.size(); j ++) {
      string value = ss[j];
      uint32_t cov = value.empty() ? 0 : boost::lexical_cast<uint32_t>(value);
      res.push_back(cov);
    }
    row ++;
  }
  pclose(stream);
  if(int(res.size()) != ncol * nrow) {
    cout << format("vector has %d elements != %d * %d\n") % res.size() % ncol % nrow;
    exit(1);
  }
  
  IntVecVec ivv;
  for(int i = 0; i < ncol; i ++) {
    IntVec iv;
    for(int j = 0; j < nrow; j ++)
      iv.push_back(res[j*ncol+i]);
    ivv.push_back(iv);
  }
  return ivv;
}
StrVec get_cov(const string& chr, const int& beg, const int& end, const string& f_cov, const int& opt_conf) {
  StrVec sv;
  IntVecVec ivv = get_raw_cov(chr, beg, end, f_cov, opt_conf);
  for(uint32_t i = 0; i < ivv.size(); i ++) {
    IntVec iv = ivv[i];
    string cov_str;
    for(uint32_t j = 0; j < iv.size(); j ++) {
      uint32_t cov = iv[j];
      cov_str += cov>=2 ? "2" : boost::lexical_cast<string>(cov);
    }
    sv.push_back( cov_str );
  }
  return sv;
}
IntStrMap get_snp(const string& chr, const int& beg, const int& end, const string& f_snp, const int& opt_conf) {
  int ncol, col_offset = 4;
  if(opt_conf == 1) {
    ncol = 84;
  } else if(opt_conf == 2) {
    ncol = 288;
  } else {
    cout << "unknown opt_conf: " << opt_conf << endl;
    exit(1);
  }
  
  IntStrMap snpMap;
  string cmd = str( format("tabix %s %s:%d-%d") % f_snp % chr % beg % end );
  int MAX_BUFSIZE = 25535;
  char buffer[MAX_BUFSIZE];
  FILE *stream = popen(cmd.c_str(), "r");
  pair< map<uint32_t, string>::iterator, bool > p;
  while( fgets(buffer, MAX_BUFSIZE, stream) != NULL ) {
    string line(buffer);
    if(line.substr(line.length()-1, 1) == "\n") line.erase(line.length()-1, 1);
    if(line.length() == 0) continue;
    vector<string> ss;
    boost::split(ss, line, boost::is_any_of("\t"));
    if( int(ss.size()) != ncol + col_offset ) {
      cout << format("not %d+%d columns: %d\n%s\n") % ncol % col_offset % ss.size() % line;
      exit(1);
    }
    uint32_t pos = boost::lexical_cast<uint32_t>(ss[1]);
    string snp = ss[3];
    for(uint32_t j = col_offset; j < ss.size(); j ++) 
      snp += ss[j];
    p = snpMap.insert( pair< uint32_t, string > (pos, snp) );
    if(p.second == false) {
      cout << format("location[%s:%d] appears >1 times\n") % chr % pos;
      exit(1);
    }
  }
  pclose(stream);
  return snpMap;
}


