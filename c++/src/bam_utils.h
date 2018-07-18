#ifndef BAM_UTILS_H
#define BAM_UTILS_H

#include <api/BamAux.h>
#include <api/BamReader.h>
#include <string>
#include <vector>
using namespace std;

namespace BamTools {

typedef vector<BamAlignment> BamAlnVector;
typedef pair<BamAlignment, BamAlignment> BamAlnPair;
typedef vector< pair<BamAlignment, BamAlignment> > BamAlnPairVector;
typedef vector< pair<string, int> > RefInfo;

class BamReader;

class Utilities {
  public:
    static bool ParseRegionString(const string& regionString, 
      const BamReader& reader, BamRegion& region);
    bool PairReads(const string& name, const BamAlnVector& als, 
      int& n_total, BamAlnPairVector& alpairs);
    bool FixReads(const string& name, BamAlnVector& als, 
      int& type, int& n_total, vector<int32_t>& iss, BamAlnPairVector& alnpairs);
    static void PrintReads(const string& name, const int& n_total, 
      BamAlnVector& als);
    static void PrintReads(const string& name, const int& n_total, 
      BamAlnPairVector& alpairs);
    static void PrintRead(const BamAlignment& al, const RefVector& refs);
    static void Reverse(string& sequence);
    static void ReverseComplement(string& sequence);
    static string getCigarString(const vector<CigarOp> cigars); 
  protected:
    static bool CheckReadPair(const BamAlignment& al1, const BamAlignment& al2, 
      int32_t& is);
};

class BamReader2 : public BamReader {
  public:
    bool GetReads(BamAlnVector& als, string& name, BamAlnVector& als2);
};

}


#endif
