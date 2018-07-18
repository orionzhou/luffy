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
#include "sam.h"

using namespace std;
using boost::format;
namespace fs = boost::filesystem;
namespace po = boost::program_options;

typedef struct {
  int beg, end;
  samfile_t *in;
} tmpstruct_t;

// callback for bam_fetch()
static int fetch_func(const bam1_t *b, void *data) {
  bam_plbuf_t *buf = (bam_plbuf_t*)data;
  bam_plbuf_push(b, buf);
  return 0;
}
// callback for bam_plbuf_init()
static int pileup_func(uint32_t tid, uint32_t pos, int n, const bam_pileup1_t *pl, void *data) {
  tmpstruct_t *tmp = (tmpstruct_t*)data;
  if ((int)pos >= tmp->beg && (int)pos < tmp->end)
    printf("%s\t%d\t%d\n", tmp->in->header->target_name[tid], pos + 1, n);
  return 0;
}

int main(int argc, char *argv[]) {
  string fi, regionStr;
  po::options_description cmdOpts("Allowed options");
  cmdOpts.add_options()
    ("help,h", "produce help message")
    ("in,i", po::value<string>(&fi), "input file")
    ("region,r", po::value<string>(&regionStr)->default_value(""), "region string")
  ;
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, cmdOpts), vm);
  po::notify(vm);
  if(vm.count("help") || !vm.count("in")) {
    cout << cmdOpts << endl;
    return 1;
  }

  tmpstruct_t tmp;
  tmp.beg = 0; tmp.end = 0x7fffffff;
  tmp.in = samopen(fi.c_str(), "rb", 0);
  if (tmp.in == 0) {
    cerr << format("Fail to open BAM file %s\n") % fi;
    return 1;
  }
  if (regionStr == "") { // if a region is not specified
    sampileup(tmp.in, -1, pileup_func, &tmp);
  } else {
    int ref;
    bam_index_t *idx;
    bam_plbuf_t *buf;
    idx = bam_index_load( fi.c_str() ); // load BAM index
    if (idx == 0) {
      cerr << "BAM indexing file is not available\n";
      return 1;
    }
    bam_parse_region(tmp.in->header, regionStr.c_str(), &ref, &tmp.beg, &tmp.end); // parse the region
    if (ref < 0) {
      cerr << format("Invalid region %s\n") % regionStr;
      return 1;
    }
    buf = bam_plbuf_init(pileup_func, &tmp); // initialize pileup
    bam_fetch(tmp.in->x.bam, idx, ref, tmp.beg, tmp.end, buf, fetch_func);
    bam_plbuf_push(0, buf); // finalize pileup
    bam_index_destroy(idx);
    bam_plbuf_destroy(buf);
  }
  samclose(tmp.in);
  
  return EXIT_SUCCESS;
}
