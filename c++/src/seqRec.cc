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
#include "geco_sequences.h"
#include "read.h"
#include "seqRec.h"
#include "ssp.h"
using namespace std;
using namespace geco;
using boost::format;
namespace fs = boost::filesystem;
namespace po = boost::program_options;

int main( int argc, char* argv[] ) {
  int opt_conf;
  string opt_ind, f_acc;
  string chr, f_id, f_out, f_gtb, d_refseq;
  clock_t time1 = clock();
  po::options_description cmdOpts("Allowed options");
  cmdOpts.add_options()
    ("help,h", "produce help message")
    ("f_acc,a", po::value<string>(&f_acc)->default_value("/project/youngn/zhoup/Scripts/conf/acc_ids.tbl"), "acc option file")
    ("opt_ind,t", po::value<string>(&opt_ind), "ind option")
    ("opt_conf,c", po::value<int>(&opt_conf)->default_value(1), "id_all tion")
    ("f_id,i", po::value<string>(&f_id), "transcript id file")
    ("f_out,o", po::value<string>(&f_out), "output file prefix")
    ("f_gtb,g", po::value<string>(&f_gtb)->default_value("/project/youngn/zhoup/Data/genome/mt_35/10_model/66_final.gtb"), "gene model table file")
    ("d_refseq,s", po::value<string>(&d_refseq)->default_value("/project/youngn/zhoup/Data/db/geco/mt_35/"), "refseq directory")
  ;
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, cmdOpts), vm);
  po::notify(vm);
  if(vm.count("help") || !vm.count("f_id") || !vm.count("f_out")) {
    cout << cmdOpts << endl;
    return 1;
  }

  string f_cov, f_snp;
  StrVec inds;
  if(opt_conf == 1) {
    f_cov = "/project/youngn/zhoup/Data/repo/mt_35/30_vnt_acc84/21_coverage/cov2.tbl.gz";
    f_snp = "/project/youngn/zhoup/Data/repo/mt_35/30_vnt_acc84/11_snps/11_merged.tbl.gz";
    inds = get_acc_ids(f_acc, "acc84"); 
  } else if(opt_conf == 2) {
    f_cov = "/project/youngn/zhoup/Data/repo/mt_35/30_vnt_acc288/21_coverage/cov2.tbl.gz";
    f_snp = "/project/youngn/zhoup/Data/repo/mt_35/30_vnt_acc288/11_snps/11_merged.tbl.gz";
    inds = get_acc_ids(f_acc, "acc288"); 
  } else {
    cout << "unknown opt_conf: " << opt_conf << endl;
    exit(1);
  }

  gLocalTxtSequenceRetriever ret(d_refseq.c_str(), "");
  gSequenceRetriever seqret(ret);

  map<string, Transcript> tm;
  if( !read_gtb(f_gtb, tm) ) { cout << format("error reading gtb file: %s\n") % f_gtb; }

  set<string> ids = read_ids(f_id);
  ofstream fho, fhl;
  fho.open((f_out+".tbl").c_str());
  fho << format("id\tacc\ttype\tlength\tnum_N\tnum_snp\tseq\n");
  fhl.open((f_out+".log").c_str());

  int counter = 0;
  for(set<string>::iterator it0 = ids.begin(); it0 != ids.end(); it0 ++) {
    string id = *it0;
    map<string, Transcript>::iterator it = tm.find(id);
    if(it == tm.end()) { cout << format("cannot find %s\n") % id; exit(1); }
    Transcript t = it->second;
    
    gTranscript gtr;
    get_geco_transcript(ret, seqret, t, gtr);
  
//    cout << format("%s\t%s\t%d\t%d\t%s\t%d\n") % id % t.chr % t.beg % t.end % t.forward % t.seq.length();
    fhl << format("%s\t%s\t%d\t%d\t%s\t%d\n") % id % t.chr % t.beg % t.end % t.forward % t.seq.length();

    StrVec covs = get_cov(t.chr, t.beg, t.end, f_cov, opt_conf);
    vector<gVariations> varVec;
    for(uint32_t i = 0; i < inds.size(); i ++) {
      gVariations var(t.chr);
      for(int pos = t.beg; pos <= t.end; pos ++) {
        string cov = covs[i].substr(pos-t.beg, 1);
        if(cov == "0")
          var.addSubstitution(t.chr, gReferenceInterval(pos-1, pos), "N");
      }
      varVec.push_back(var);
    }

    IntStrMap snpMap = get_snp(t.chr, t.beg, t.end, f_snp, opt_conf);
    for(IntStrMap::iterator it = snpMap.begin(); it != snpMap.end(); it ++) {
      uint32_t pos = it->first;
      uint32_t posL = pos - t.beg + 1;
      string snp = it->second;
      string ref = snp.substr(0, 1);
      if(ref != t.seq.substr(posL-1, 1)) {
        cout << format("  %d[%d]: ref[%s] != %s\n") % pos % posL % ref % t.seq.substr(posL-1, 1);
        exit(1);
      }

      for(uint32_t i = 0; i < inds.size(); i ++) {
        string alt = snp.substr(i+1, 1);
        string cov = covs[i].substr(posL-1,1);
        if(alt == "N") {
          if(cov == "1") {
//              fhl << format("\t1\t%d[%d]\t%s\tref[%s]\talt[%s]\tcov[%s]\n") % inds[i] % pos % posL % ref % alt % cov;
            varVec[i].addSubstitution(t.chr, gReferenceInterval(pos-1, pos), alt);
          }
        } else {
          if(cov == "0") {
            if(ref == alt) {
              fhl << format("\t2\t%d[%d]\t%s\tref[%s]\talt[%s]\tcov[%s]\n") % inds[i] % pos % posL % ref % alt % cov;
            } else {
              fhl << format("\t3\t%d[%d]\t%s\tref[%s]\talt[%s]\tcov[%s]\n") % inds[i] % pos % posL % ref % alt % cov;
            }
          }
          if(alt != ref) { 
            varVec[i].addSubstitution(t.chr, gReferenceInterval(pos-1, pos), alt);
          }
        }
      }
    }

    for(unsigned i = 0; i < inds.size(); i ++) {
      gVariations var = varVec[i];
      gElement gel(t.chr, t.beg-1, t.end, t.forward, seqret, gElm, var);
      string seqStr = gel.getSequence();
      boost::to_upper(seqStr);

      gArray<gPos> possG(0, var.getCount(), true);
      gArray<gPos> possL(0, var.getCount(), true);
      for(gSize j = 0; j < var.getCount(); j ++)
        possG.setValue(j, var.getVariationStart(j));
      if(var.getCount() > 0) 
        possL = gel.getElementPositionsFromReference(possG);
      vector<string> types = get_region_from_pos(t, possL);
//        cout << possG << endl;
//        cout << possL << endl;
//        cout << boost::join(types, "  ") << endl;

      map<string, LocPairVec>::iterator it;
      map<string, int> num_Ns, num_snps;
      for(it = t.locLs.begin(); it != t.locLs.end(); it ++) {
        string type = it->first;
        num_Ns[type] = 0;
        num_snps[type] = 0;
      }
      for(gSize j = 0; j < var.getCount(); j ++) {
        string type = types[j];
        if(var.getVariationSequence(j) == "N")
          num_Ns[type] ++;
        else
          num_snps[type] ++;
      }

      for(it = t.locLs.begin(); it != t.locLs.end(); it ++) {
        string type = it->first;
        LocPairVec loc = it->second;
        string seq = get_subseq_from_loc(seqStr, loc, true);

        fho << format("%s\t%s\t%s\t%d\t%d\t%d\t%s\n") % id % inds[i] % type % seq.length() % num_Ns[type] % num_snps[type] % seq;
      }
    }
    if( ++counter % 100 == 0 )
      cout << right << setw(50) << format("%d records: %.01f minutes\n") % counter % ( (double)(clock() - time1) / ((double)CLOCKS_PER_SEC * 60) );
  }
  fho.close();

  return EXIT_SUCCESS;
}



