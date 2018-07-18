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

struct GalRecord {
  string id;
  string tId;
  uint32_t tBeg;
  uint32_t tEnd;
  uint32_t tSize;
  string tSrd;
  
  string qId;
  uint32_t qBeg;
  uint32_t qEnd;
  uint32_t qSize;
  string qSrd;
  
  uint32_t ali;
  uint32_t mat;
  uint32_t mis;
  uint32_t qN;
  uint32_t tN;
  float ident;
  float score;

  LocVec qLoc;
  LocVec tLoc;
};
typedef vector<GalRecord> GalRecords;

void check_gal_record(const GalRecord& br) {
  if(br.qBeg > br.qEnd)
    cerr << format("%s: qBeg[%d] > qEnd[%d]\n") % br.qId % br.qBeg % br.qEnd;
  if(br.tBeg > br.tEnd)
    cerr << format("%s: tBeg[%d] > tEnd[%d]\n") % br.tId % br.tBeg % br.tEnd;
  if(br.score <= 0)
    cerr << format("%s:%d-%d score:%g\n") % br.qId % br.qBeg % br.qEnd % br.score;
  if(br.qLoc.size() != br.tLoc.size())
    cerr << format("unequal blocks: %s:%d-%d\n") % br.qId % br.qBeg % br.qEnd; 
  if(locVecLen(br.qLoc) != locVecLen(br.tLoc))
    cerr << format("unequal alnlen: %s:%d-%d\n") % br.qId % br.qBeg % br.qEnd; 
}
GalRecord make_gal_record(const vector<string>& ss) {   
  using boost::lexical_cast;
  
  GalRecord br;
  br.id = ss[0];

  br.tId = ss[1];
  br.tBeg = lexical_cast<uint32_t>(ss[2]);
  br.tEnd = lexical_cast<uint32_t>(ss[3]);
  br.tSrd = ss[4];
  br.tSize = lexical_cast<uint32_t>(ss[5]);
  
  br.qId = ss[6];
  br.qBeg = lexical_cast<uint32_t>(ss[7]);
  br.qEnd = lexical_cast<uint32_t>(ss[8]);
  br.qSrd = ss[9];
  br.qSize = lexical_cast<uint32_t>(ss[10]);
   
  br.ali = ss[11].empty() ? 0 : lexical_cast<uint32_t>(ss[11]);
  br.mat = ss[12].empty() ? 0 : lexical_cast<uint32_t>(ss[12]);
  br.mis = ss[13].empty() ? 0 : lexical_cast<uint32_t>(ss[13]);
  br.qN  = ss[14].empty() ? 0 : lexical_cast<uint32_t>(ss[14]);
  br.tN  = ss[15].empty() ? 0 : lexical_cast<uint32_t>(ss[15]);
  br.ident = ss[16].empty() ? 0 : lexical_cast<float>(ss[16]);
  br.score = ss[17].empty() ? lexical_cast<float>(br.mat) :  lexical_cast<float>(ss[17]);

  br.qLoc = locStr2Vec(ss[18]);
  br.tLoc = locStr2Vec(ss[19]);

  check_gal_record(br);
  return br;
}

void gal_tiling(const GalRecords& grs, string& qId, ofstream& fho, const unsigned& len_min) {
  LocVec lv1;
  int i = 0;
  for(GalRecords::const_iterator it = grs.begin(); it != grs.end(); it++) {
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
    
    GalRecord gr = grs[ loc.idx_ext ];
    string id(gr.id), qSrd("+");
    string tSrd = gr.qSrd == gr.tSrd ? "+" : "-";
    uint32_t tSize(gr.tSize), qSize(gr.qSize);
    float ident(gr.ident), score(gr.score);
     
    string tId = gr.tId;
    uint32_t rqb(qBeg-gr.qBeg+1), rqe(qEnd-gr.qEnd+1);
    Location rqloc;
    rqloc.beg = rqb;
    rqloc.end = rqe;
    LocVec rqLoc;
    rqLoc.push_back(rqloc);
    
    LocVec qLoc = loc_intersect(gr.qLoc, rqLoc);
    rqb = qLoc.begin()->beg;
    rqe = qLoc.rbegin()->end;
    qBeg = gr.qBeg + rqb - 1;
    qEnd = gr.qEnd + rqe - 1;
    
    uint32_t rtb = coordTransform(rqb, gr.qLoc, "+", gr.tLoc, "+"); //need fill
    uint32_t rte = coordTransform(rqe, gr.qLoc, "+", gr.tLoc, "+");
    Location rtloc;
    rtloc.beg = rtb;
    rtloc.end = rte;
    LocVec rtLoc;
    rtLoc.push_back(rtloc);

    LocVec tLoc = loc_intersect(gr.qLoc, rtLoc);
    uint32_t tBeg = tSrd=="+" ? gr.tBeg+rtb-1 : gr.tEnd-rte+1;
    uint32_t tEnd = tSrd=="+" ? gr.tBeg+rte-1 : gr.tEnd-rtb+1;
    fho << format("%s\t%s\t%d\t%d\t%s\t%d\t%s\t%d\t%d\t%s\t%d \
      \t\t\t\t\t\t%g\t%g\t%s\t%s\n") 
      % id % tId % tBeg % tEnd % tSrd % tSize
           % qId % qBeg % qEnd % qSrd % qSize 
      % ident % score % locVec2Str(qLoc) % locVec2Str(tLoc);
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
    ("min,m", po::value<unsigned>(&len_min)->default_value(1), "mininum tiling length")
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
  fho << "id\ttId\ttBeg\ttEnd\ttSrd\ttSize"
    << "\tqId\tqBeg\tqEnd\tqSrd\tqSize"
    << "\tali\tmat\tmis\tqN\ttN\tident\tscore\tqLoc\ttLoc" << endl;
  
  GalRecords grs;
  string qId_p = "";
  string line;
  while(fhi.good()) {
    getline(fhi, line);
    if(line.length() == 0) continue;
    boost::erase_all(line, " ");
    
    vector<string> ss;
    boost::split(ss, line, boost::is_any_of("\t"));
    if(ss[0] == "id") continue;
    
    GalRecord gr = make_gal_record(ss);
    if(gr.qId == qId_p) {
      grs.push_back(gr);
    } else {
      if(qId_p != "") {
        gal_tiling(grs, qId_p, fho, len_min);
      }
      grs.clear();
      grs.push_back(gr);
      qId_p = gr.qId;
    }
  }
  if(grs.size() > 0) {
    gal_tiling(grs, qId_p, fho, len_min);
  }

  cout << right << setw(60) << format("(running time: %.01f minutes)\n") % ( (double)(clock() - time1) / ((double)CLOCKS_PER_SEC * 60) );
  return EXIT_SUCCESS;
}
