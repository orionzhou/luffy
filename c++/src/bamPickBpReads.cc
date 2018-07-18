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

struct BreakPoint {
  int id_p;
  string chr, type_p, ins;
  int beg, end, beg_r, end_r;
  int size_d, size_i;
  int n_ind_all, n_ind, n_reads, n_reads_uniq;
  set<string> inds;
  int id_c;
  string chr_l, chr_r, strand_l, strand_r, type_c, acc;
  int pos_l, pos_r, len, n_acc;
};
struct ReadInfo {
  int id;
  string ind, name;
  bool isFirst;
  int matePos;
  bool isReverse;
  vector<CigarOp> cigar_ops;
};

bool read_reads(const string& fi, map<string, ReadInfo>& rim, vector<string>& pres) {
  ifstream fhi(fi.c_str());
  if(!fhi.is_open()) {
    cout << format("cannot open file: %s") % fi;
    return false;
  }

  pair< map<string, ReadInfo>::iterator, bool > rim_p;
  map<string, int> pm;
  map<string, int>::iterator pm_it;
  pair< map<string, int>::iterator, bool > pm_p;

  boost::regex e("(\\d+)([MIDNSHP])");
  boost::match_results<string::const_iterator> m;
  boost::match_flag_type flags = boost::match_default;
  string line;
  getline(fhi, line);
  while(fhi.good()) {
    getline(fhi, line);
    if(line.length() == 0) continue;
    vector<string> ss;
    boost::split(ss, line, boost::is_any_of("\t"));
    ReadInfo ri;
    ri.id = atoi(ss[0].c_str());
    ri.ind = ss[1];
    ri.name = ss[2];
    ri.isFirst = ss[3] == "1" ? true : false;
    ri.matePos = atoi(ss[4].c_str());
    ri.isReverse = ss[6] == "-" ? true : false;

    string cigar = ss[5];
    string::const_iterator start = cigar.begin();
    string::const_iterator end = cigar.end();
    vector<CigarOp> cigar_ops;
    while(boost::regex_search(start, end, m, e, flags)) {
      uint32_t length = atoi(string(m[1].first, m[1].second).c_str());
      char type = *(m[2].first);
      CigarOp co (type, length);
      cigar_ops.push_back(co);
      start = m[0].second;
    }
    ri.cigar_ops = cigar_ops;

    string key = ri.name + "." + ri.ind + "." + ss[0];
    rim_p = rim.insert( pair<string, ReadInfo> (key, ri) );
    if(rim_p.second == false) {
      cout << format("read[%s] is duplicated in %s\n") % key % fi;
      return false;
    }
    pm_p = pm.insert( pair<string, int> (ri.ind, 1) );
  }

  for(pm_it = pm.begin(); pm_it != pm.end(); pm_it ++)
    pres.push_back(pm_it->first);
  cout << format("%d records read\n") % rim.size();
  cout << format("%d samples\n") % pres.size();
  return true;
}
bool read_breakpoints(const string& fi, vector<BreakPoint>& bv, set<string>& inds_all) {
  ifstream fhi(fi.c_str());
  if(!fhi.is_open()) {
    cout << format("cannot open file: %s") % fi;
    exit(0);
  }
  
  string line;
  getline(fhi, line);

  while(fhi.good()) {
    getline(fhi, line);
    if(line.length() == 0) continue;
    vector<string> ss;
    boost::split(ss, line, boost::is_any_of("\t"));
    BreakPoint bp;
    bp.id_p = atoi(ss[0].c_str());
    bp.chr = ss[1];
    bp.beg = atoi(ss[2].c_str());
    bp.end = atoi(ss[3].c_str());
    bp.beg_r = atoi(ss[4].c_str());
    bp.end_r = atoi(ss[5].c_str());
    bp.type_p = ss[6];
    bp.size_d = atoi(ss[7].c_str());
    bp.ins = ss[8];
    bp.size_i = atoi(ss[9].c_str());
    bp.n_ind_all = atoi(ss[10].c_str());
    bp.n_ind = atoi(ss[11].c_str());
    bp.n_reads = atoi(ss[12].c_str());
    bp.n_reads_uniq = atoi(ss[13].c_str());
    
    set<string> inds;
    vector<string> ps, pps;
    boost::split(ps, ss[14], boost::is_any_of(" "));
    for(unsigned i = 0; i < ps.size(); i ++) {
      boost::split(pps, ps[i], boost::is_any_of(":"));
      string ind = pps[0];
      inds.insert(ind);
      if(inds_all.find(ind) == inds_all.end()) inds_all.insert(ind);
    }
    bp.inds = inds;
    bp.id_c = atoi(ss[15].c_str());
    bp.chr_l = ss[16];
    bp.pos_l = atoi(ss[17].c_str());
    bp.strand_l = ss[18];
    bp.chr_r = ss[19];
    bp.pos_r = atoi(ss[20].c_str());
    bp.strand_r = ss[21];
    bp.type_c = ss[22];
    bp.acc = ss[23];
    bp.len = atoi(ss[24].c_str());
    bp.n_acc = atoi(ss[25].c_str());

    bv.push_back(bp);
  }
  cout << format("%d break_points in %d inds\n") % bv.size() % inds_all.size();
  return true;
}

int main(int argc, char *argv[]) {
  string d_bam, f_bp, f_out;
  clock_t time1 = clock();
  po::options_description cmdOpts("Allowed options");
  cmdOpts.add_options()
    ("help,h", "produce help message")
    ("d_bam,b", po::value<string>(&d_bam), "BAM directory")
    ("f_bp,p", po::value<string>(&f_bp), "break point file")
    ("f_out,o", po::value<string>(&f_out), "output file")
  ;
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, cmdOpts), vm);
  po::notify(vm);
  if(vm.count("help") || !vm.count("d_bam") || !vm.count("f_bp") || !vm.count("f_out")) {
    cout << cmdOpts << endl;
    return 1;
  }

  Utilities utils;
  BamWriter writer;
  map<string, BamReader> readers;
  map<string, BamReader>::iterator readers_it;
  ofstream fho( f_out.c_str() );
  fho << "ind\tid_pindel\ttype\tbeg\tcigar\tid\tfirst\tseq\n";
 
  vector<BreakPoint> bv;
  set<string> inds_all;
  read_breakpoints(f_bp, bv, inds_all);
  for(set<string>::iterator it = inds_all.begin(); it != inds_all.end(); it ++) {
    string ind = *it;
    BamReader reader;
    string f_bam = d_bam + "/" + ind + ".bam";
    if(!reader.Open(f_bam) || !reader.OpenIndex(f_bam+".bai") || !reader.HasIndex()) {
      cout << "error opening " << f_bam << endl;
      exit(1);
    }
//    writer.Open( fo, reader.GetHeader(), reader.GetReferenceData() );
    cout << "working on " << ind << endl;
    
    for(unsigned i = 0; i < bv.size(); i ++) {
      BreakPoint bp = bv[i];
      string regionStr = str( format("%s:%d..%d") % bp.chr % (bp.beg-1000) % (bp.end+1000) );

      BamRegion region;
      if(!utils.ParseRegionString(regionStr, reader, region)) { cout << "invalid string: " << regionStr << endl; }
      if(!reader.SetRegion(region)) { cout << "cannot set region\n"; }
      set<string> inds = bp.inds;
      
      BamAlignment al;
      while( reader.GetNextAlignment(al) ) {
/*        rim_it = rim.find(str(format("%s.%s.%d") % al.Name % ind % it2->id));
        if(rim_it != rim.end()) {
          cnt ++;
          ReadInfo ri = rim_it->second;
          assert(ri.ind == ind);
          if(ri.id != it2->id) {
            cout << format("bp_id: %d != %d at %s %s\n") % it2->id % ri.id % ind % al.Name;
            exit(1);
          }
          if(al.IsMapped() && !al.IsMateMapped()) {
            assert(al.RefID==chrId);
            al.SetIsMateMapped(true);
            al.MateRefID = chrId;
            al.MatePosition = ri.matePos;
            if(!ri.isReverse) al.SetIsMateReverseStrand(true);
          } else {
            assert(!al.IsMapped() && al.IsMateMapped() && al.MateRefID == chrId);
            al.SetIsMapped(true);
            al.RefID = chrId;
            al.Position = ri.matePos;
            if(!ri.isReverse) {
              util.Reverse(al.Qualities);
              util.ReverseComplement(al.QueryBases);
              al.SetIsReverseStrand(true);
            }
            al.MapQuality = 40;
            al.CigarData = ri.cigar_ops;
          }
          al.SetIsProperPair(true);
          al.InsertSize = abs(al.Position - al.MatePosition) + 90;
          al.AddTag("AC", "Z", ind);
          al.AddTag("BP", "i", it2->id);
          al.AddTag("BT", "i", 1);
          writer.SaveAlignment(al); */
        int type = 0;
        if(!al.IsMapped() && al.IsMateMapped()) {
          type = 1;
        } else {
          int beg = al.Position, end = al.GetEndPosition();
          if( (beg < bp.beg && beg+90>bp.beg) || (end > bp.end && end-90 < bp.end) ) type = 2;
        }
        if(type > 0) {
          fho << format("%s\t%d\t%d\t%d\t%s\t%s\t%d\t%s\n") % ind % bp.id_p 
            % type % al.Position % utils.getCigarString(al.CigarData) % al.Name % al.IsFirstMate() % al.QueryBases;
        }
      }
      if(i < 0) break;
    }
  }
  
  cout << right << setw(60) << format("(running time: %.01f minutes)\n") % ( (double)(clock() - time1) / ((double)CLOCKS_PER_SEC * 60) );
  return EXIT_SUCCESS;
}
