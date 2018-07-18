#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <fstream>
#include <vector>
#include <map>
#include <cmath>
#include <time.h>
#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <boost/regex.hpp>
#include <boost/program_options.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include "location.h"
using namespace std;
using boost::format;
namespace fs = boost::filesystem;
namespace po = boost::program_options;

struct BlastRecord {
  string qId;
  uint32_t qBeg;
  uint32_t qEnd;
  string qSrd;
  
  string tId;
  uint32_t tBeg;
  uint32_t tEnd;
  string tSrd;
  
  float ident;
  double e;
  float score;
  uint32_t len;
  uint32_t match;
  uint32_t misMatch;
  uint32_t gap;
};
typedef vector<BlastRecord> BlastRecords;

BlastRecord make_blast_record(const vector<string>& ss) {   
  BlastRecord br;
  br.qId = ss[0];
  br.tId = ss[1];
  br.ident = boost::lexical_cast<float>(ss[2]);
  br.len = boost::lexical_cast<uint32_t>(ss[3]);
  br.misMatch = boost::lexical_cast<uint32_t>(ss[4]);
  br.gap = boost::lexical_cast<uint32_t>(ss[5]);
  br.qBeg = boost::lexical_cast<uint32_t>(ss[6]);
  br.qEnd = boost::lexical_cast<uint32_t>(ss[7]);
  br.tBeg = boost::lexical_cast<uint32_t>(ss[8]);
  br.tEnd = boost::lexical_cast<uint32_t>(ss[9]);
  br.qSrd = "+";
  br.tSrd = "+";
  if(br.tBeg > br.tEnd) {
    swap(br.tBeg, br.tEnd);
    br.qSrd = "-";
  }
  if(br.qBeg > br.qEnd)
    cerr << format("qBeg[%d] > qEnd[%d]\n") % br.qBeg % br.qEnd;
  br.e = boost::lexical_cast<double>(ss[10]);
  br.score = boost::lexical_cast<float>(ss[11]);
  return br;
}

void blast_tiling(const BlastRecords& brs, string& qId, ofstream& fho, const unsigned& len_min) {
  LocVec lv1;
  int i = 0;
  for(BlastRecords::const_iterator it = brs.begin(); it != brs.end(); it++) {
    Location loc;
    loc.beg = it->qBeg;
    loc.end = it->qEnd;
    loc.score = it->score;
    loc.idxs.insert( i++ );
    lv1.push_back(loc);
  }

  LocVec lv2 = tiling(lv1, true);
  for(LocVec::const_iterator it = lv2.begin(); it != lv2.end(); it++) {
    Location loc = *it;
    uint32_t qBeg(loc.beg), qEnd(loc.end);
    uint32_t qLen = qEnd - qBeg + 1;
    if(qLen < len_min) continue;
    int idx = loc.idx_ext;
    BlastRecord br = brs[idx];
    string qSrd = br.qSrd;
    string tSrd = br.tSrd;
    float ident(br.ident), score(br.score);
    double e(br.e);
    string tId = br.tId;
    uint32_t tBeg = round( br.tBeg + (qBeg-br.qBeg)*(br.tEnd-br.tBeg)/(br.qEnd-br.qBeg) );
    uint32_t tEnd = round( br.tBeg + (qEnd-br.qBeg)*(br.tEnd-br.tBeg)/(br.qEnd-br.qBeg) );
    uint32_t tLen = tEnd - tBeg + 1;
    fho << format("%s\t%d\t%d\t%s\t%d\t%s\t%d\t%d\t%s\t%d\t%g\t%g\t%g\n") 
      % qId % qBeg % qEnd % qSrd % qLen % tId % tBeg % tEnd % tSrd % tLen 
      % ident % e % score;
  }
  cout << qId << endl;
}

int main( int argc, char* argv[] ) {
  string fi, fo;
  unsigned len_min;
  clock_t time1 = clock();
  po::options_description cmdOpts("Allowed options");
  cmdOpts.add_options()
    ("help,h", "produce help message")
    ("in,i", po::value<string>(&fi), "input file")
    ("out,o", po::value<string>(&fo), "output file")
    ("min,m", po::value<unsigned>(&len_min)->default_value(100), "mininum tiling length")
  ;
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, cmdOpts), vm);
  po::notify(vm);
  if(vm.count("help") || !vm.count("in") || !vm.count("out")) {
    cout << cmdOpts << endl;
    return 1;
  }

  ifstream fhi(fi.c_str());
  if(!fhi.is_open()) { cout << format("cannot read: %s\n") % fi; return false; }
  ofstream fho(fo.c_str());
  if(!fho.is_open()) { cout << format("cannot write: %s\n") % fo; return false; }
  fho << "qId\tqBeg\tqEnd\tqSrd\tqLen"
    << "\ttId\ttBeg\ttEnd\ttSrd\ttLen"
    << "\tident\te\tscore" << endl;
  
  BlastRecords brs;
  string qId_p = "";
  string line;
  while(fhi.good()) {
    getline(fhi, line);
    if(line.length() == 0) continue;
    boost::erase_all(line, " ");
    
    vector<string> ss;
    boost::split(ss, line, boost::is_any_of("\t"));
    BlastRecord br = make_blast_record(ss);
    if(br.qId == qId_p) {
      brs.push_back(br);
    } else {
      if(qId_p != "") {
        blast_tiling(brs, qId_p, fho, len_min);
      }
      brs.clear();
      brs.push_back(br);
      qId_p = br.qId;
    }
  }
  if(brs.size() > 0) {
    blast_tiling(brs, qId_p, fho, len_min);
  }

  cout << right << setw(60) << format("(running time: %.01f minutes)\n") % ( (double)(clock() - time1) / ((double)CLOCKS_PER_SEC * 60) );
  return EXIT_SUCCESS;
}
