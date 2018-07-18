#include <string>
#include <vector>
using namespace std;

typedef vector<int> IntVec;
typedef set<int> IntSet;
struct Location {
    uint32_t beg;
    uint32_t end;
    double score;
    int idx_ext;
    IntSet idxs;
};
typedef vector< Location > LocVec;
typedef map< uint32_t, Location > LocMap;

uint32_t locVecLen (const LocVec& locv);
LocVec locStr2Vec (const string& locS); 
string locVec2Str (const LocVec& locV);

LocVec loc_disjoin(const LocVec& lv);
LocVec loc_reduce(const LocVec& lv);
LocVec loc_intersect(const LocVec& lv1, const LocVec& lv2);
LocVec loc_union(const LocVec& lv1, const LocVec& lv2);

bool compare_loc (const Location& l1, const Location& l2);
uint32_t coordTransform(uint32_t& pos, const LocVec& lvi, const string& srdi, const LocVec& lvo, const string& srdo);
LocVec tiling (const LocVec& lvi,  const bool& flag_max = false);

