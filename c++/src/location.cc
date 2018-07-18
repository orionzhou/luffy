#include <fstream>
#include <iostream>
#include <algorithm>
#include <vector>
#include <set>
#include <map>
#include <string>
#include <assert.h>
#include <limits>

#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include "location.h"
using namespace std;
using boost::format;
using boost::lexical_cast;

uint32_t locVecLen (const LocVec& locv) {
  uint32_t len = 0;
  for(vector<Location>::const_iterator it = locv.begin(); it != locv.end(); it++) {
    Location loc = *it;
    if(loc.beg > loc.end)
      cerr << format("error: beg[%d] > end[%d]\n") % loc.beg % loc.end;
    len += loc.end - loc.beg + 1;
  }
  return len;
}
LocVec locStr2Vec (const string& locS) {
  LocVec locV;
  vector<string> ss, ss1;
  boost::split(ss, locS, boost::is_any_of(","));
  for(vector<string>::const_iterator it = ss.begin(); it != ss.end(); it++) {
    string locS1 = *it;
    boost::split(ss1, locS1, boost::is_any_of("-_"));
    if(ss1.size() != 2)
      cerr << format("locstr not 2 parts: %s\n") % locS1;
    Location loc;
    loc.beg = lexical_cast<uint32_t>(ss1[0]);
    loc.end = lexical_cast<uint32_t>(ss1[1]);
    locV.push_back( loc );
  }
  return locV;
}
string locVec2Str (const LocVec& locV) {
  vector<string> ss;
  for(LocVec::const_iterator it = locV.begin(); it != locV.end(); it++) {
    Location loc = *it;
    string ss1 = str( format("%d-%d") % loc.beg % loc.end );
    ss.push_back(ss1);
  }
  return boost::algorithm::join(ss, ",");
}

bool compare_loc (const Location& l1, const Location& l2) {
  if(l1.beg < l2.beg) return true;
  if(l1.beg > l2.beg) return false;

  if(l1.end < l2.end) return true;
  if(l1.end > l2.end) return false;
  return false;
}
void print_locm (const map<uint32_t, Location>& locm) {
  cout << "----------" << endl;
  for(map<uint32_t, Location>::const_iterator it = locm.begin(); it != locm.end(); it++) {
    Location loc = it->second;
    vector<string> idxs;
    for(IntSet::iterator its = loc.idxs.begin(); its != loc.idxs.end(); its++)
      idxs.push_back( boost::lexical_cast<string>(*its) );
    string idx_str = boost::algorithm::join(idxs, ",");
    cout << format("%10d\t%10d\t%g\t[%d] %s\n") % loc.beg % loc.end % 
      loc.score % loc.idx_ext % idx_str;
  }
  cout << "----------" << endl;
}
void print_locv (const vector<Location>& locv) {
  cout << "----------" << endl;
  for(vector<Location>::const_iterator it = locv.begin(); it != locv.end(); it++) {
    Location loc = *it;
    vector<string> idxs;
    for(IntSet::iterator its = loc.idxs.begin(); its != loc.idxs.end(); its++)
      idxs.push_back( boost::lexical_cast<string>(*its) );
    string idx_str = boost::algorithm::join(idxs, ",");
    cout << format("%10d\t%10d\t%g\t[%d] %s\n") % loc.beg % loc.end % 
      loc.score % loc.idx_ext % idx_str;
  }
  cout << "----------" << endl;
}


void insert_loc (const Location& loc, set<uint32_t>& begs, map<uint32_t, Location>& locm) {
  uint32_t beg(loc.beg), end(loc.end);
  set<uint32_t>::iterator its;
  map<uint32_t, Location>::iterator itmp, itms;
  its = begs.find(beg);
  //cout << loc.beg << " " << loc.end << " " << *(loc.idxs.begin()) << endl;
  if(loc.beg > loc.end) exit(EXIT_FAILURE);
  if(its == begs.end()) {
    its = begs.upper_bound(beg);
    bool flag_beg = (its == begs.begin());
    bool flag_end = (its == begs.end());
    uint32_t end_p(0), beg_s(numeric_limits<uint32_t>::max());
    if(flag_beg && !flag_end) {
      itms = locm.find( *its );
      beg_s = itms->second.beg;
    } else if(!flag_beg && flag_end) {
      itmp = locm.end();
      itmp --;
      end_p = itmp->second.end;
    } else if(!flag_beg && !flag_end) {
      itms = locm.find( *its );
      itmp = itms;
      itmp --;
      beg_s = itms->second.beg;
      end_p = itmp->second.end;
    }
    if(end_p < beg && end < beg_s) {
      Location locn = loc;
      begs.insert( locn.beg );
      locm.insert( pair<uint32_t, Location> (locn.beg, locn) );
    } else if(end_p >= beg && end < beg_s) {
      if(end_p < end) {
        Location locn1(loc), locn2(loc);
        locn1.end = end_p;
        locn2.beg = end_p + 1;
        insert_loc(locn1, begs, locm);
        insert_loc(locn2, begs, locm);
      } else {
        Location locn1(loc), locn2(itmp->second);
        itmp->second.end = beg - 1;
        locn1.idxs.insert(itmp->second.idxs.begin(), itmp->second.idxs.end());
        insert_loc(locn1, begs, locm);
        if(end_p > end) {
          locn2.beg = end + 1;
          insert_loc(locn2, begs, locm);
        }
      }
    } else if(end_p < beg && end >= beg_s) {
      Location locn1(loc), locn2(loc);
      locn1.end = beg_s - 1;
      locn2.beg = beg_s;
      insert_loc(locn1, begs, locm);
      insert_loc(locn2, begs, locm);
    } else { // end_p >= beg && end >= beg_s
      Location locn1(loc), locn2(loc);
      locn1.end = end_p;
      locn2.beg = beg_s;
      insert_loc(locn1, begs, locm);
      insert_loc(locn2, begs, locm);
      if(end_p + 1 < beg_s) {
        Location locn3 = loc;
        locn3.beg = end_p + 1;
        locn3.end = beg_s - 1;
        insert_loc(locn3, begs, locm);
      }
    }
  } else {
    itms = locm.find(beg);
    if(end < itms->second.end) {
      Location locn = itms->second;
      itms->second.end = end;
      itms->second.idxs.insert( loc.idxs.begin(), loc.idxs.end() );
      locn.beg = end + 1;
      insert_loc(locn, begs, locm);
    } else {
      itms->second.idxs.insert( loc.idxs.begin(), loc.idxs.end() );
      if(end > itms->second.end) {
        Location locn = loc;
        locn.beg = itms->second.end + 1;
        insert_loc(locn, begs, locm);
      }
    }
  }
}
LocVec loc_disjoin(const LocVec& lv) {
  LocMap lm;
  set<uint32_t> begs;
  for(LocVec::const_iterator it = lv.begin(); it != lv.end(); it++) {
    Location loc = *it;
    insert_loc(loc, begs, lm);
  }
  LocVec lvo;
  return lvo;
}
LocVec loc_reduce(const LocVec& lv) {
  LocVec lvo;
  return lvo;
}
LocVec loc_intersect(const LocVec& lv1, const LocVec& lv2) {
  LocVec lv;
  return lv;
}
LocVec loc_union(const LocVec& lv1, const LocVec& lv2) {
  LocVec lv;
  return lv;
}

uint32_t coordTransform(uint32_t& pos, const LocVec& lvi, const string& srdi, const LocVec& lvo, const string& srdo) {
  return pos;
}

LocVec tiling(const LocVec& lvi, const bool& flag_max) {
  LocVec lv = lvi;
  //sort(lv.begin(), lv.end(), compare_loc);
  
  set<uint32_t> begs;
  set<uint32_t>::iterator its;
  map<uint32_t, Location> locm;
  map<uint32_t, Location>::iterator itm;
  for (LocVec::iterator it = lv.begin(); it != lv.end(); it++)
    insert_loc( *it, begs, locm );
  
  LocVec lvo;
  vector<bool> tags;
  int idx_ext_p = -1;
  uint32_t end_p = 0;
  for(map<uint32_t, Location>::const_iterator it = locm.begin(); it != locm.end(); it++) {
    Location loc = it->second;
    IntSet idxs = loc.idxs;
    
    int idx_ext = -1;
    double score_ext = flag_max ? 
      numeric_limits<double>::min() : numeric_limits<double>::max();
    for (set<int>::iterator it = idxs.begin(); it != idxs.end(); it ++) {
      int idx = *it;
      Location locO = lvi[idx];
      double score = locO.score;
      if( (flag_max && score>score_ext) || (!flag_max && score<score_ext) ) {
        idx_ext = idx;
        score_ext = score;
      }
    }
    loc.idx_ext = idx_ext;
    loc.score = score_ext;
    if(idx_ext == -1)
      cerr<<"error idx_ext\t"<<lvi[idx_ext].beg<<"\t"<<lvi[idx_ext].end<<endl;

    if( loc.beg == end_p + 1 && loc.idx_ext == idx_ext_p ) {
      lvo.back().end = loc.end;
      lvo.back().idxs.insert( loc.idxs.begin(), loc.idxs.end() );
    } else {
      lvo.push_back( loc );
    }
    idx_ext_p = loc.idx_ext;
    end_p = loc.end;
  }
  //print_locv(lvo);
  return lvo;
}

