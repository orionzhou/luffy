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
using namespace std;
using boost::format;
namespace fs = boost::filesystem;
namespace po = boost::program_options;

typedef pair< int, int > LocPair;
typedef vector< pair<int, int> > LocPairVec;
typedef vector<uint32_t> IntVec;
typedef vector<string> StrVec;
typedef map< uint32_t, string > IntStrMap;
typedef map< string, string > StrStrMap;
typedef map< string, vector<uint32_t> > StrVecMap;
typedef vector< vector<uint32_t> > IntVecVec;

template< class T >
struct next {
    T seed;
    next( T seed ) : seed(seed) { }
    T operator()() {
        return seed++;
    }
};
struct Loc {
    string id, chr;
    uint32_t beg, end, length;
    bool forward;
    string type, note;
};
struct LocMap {
    string name1, name2;
    uint32_t id1, start1, end1, id2, start2, end2;
};
struct Transcript {
    string id, id_gene, chr;
    bool forward;
    int frame, beg, end, cdsBeg, cdsEnd;
    LocPairVec locG_cds, locG_exon, locG_utr5, locG_utr3, locG_intron;
    LocPairVec locL_cds, locL_exon, locL_utr5, locL_utr3, locL_intron;
    map<string, LocPairVec> locLs;
    string source, type, conf, note, seq;
};

bool read_refmap(const string& fi, map<string, LocMap>& refMap);
bool read_gtb(const string& f_gtb, map<string, Transcript>& tm);
set<string> read_ids(const string& fi);
vector<Loc> read_locs(const string& fi);
IntVecVec get_raw_cov(const string& chr, const int& beg, const int& end, const string& f_cov, const int& opt_conf);
StrVec get_cov(const string& chr, const int& beg, const int& end, const string& f_cov, const int& opt_conf);
IntStrMap get_snp(const string& chr, const int& beg, const int& end, const string& f_snp, const int& opt_conf);


