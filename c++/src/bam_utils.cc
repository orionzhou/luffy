#include <cstdlib>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <vector>
#include <map>
#include <string>
#include <assert.h>
#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>
#include <api/BamReader.h>
#include "bam_utils.h"
using namespace std;
using boost::format;
using namespace BamTools;

namespace BamTools {
const char REVCOMP_LOOKUP[] = {'T', 0, 'G', 'H', 0, 0, 
  'C', 'D', 0, 0, 0, 0, 'K', 'N', 0, 0, 0, 'Y', 'W', 'A', 'A', 'B', 
  'S', 'X', 'R', 0 };
}

bool Utilities::ParseRegionString(const std::string& regionString, 
    const BamReader& reader, BamRegion& region) {
  if(regionString.empty()) return false;
  size_t foundFirstColon = regionString.find(':');
  string chrS, chrE;
  int posS, posE;
  if(foundFirstColon == string::npos) {
    chrS = regionString;
    posS = 0;
    chrE = regionString;
    posE = -1;
  } else {
    chrS = regionString.substr(0, foundFirstColon);
    size_t foundRangeDots = regionString.find("..", foundFirstColon + 1);
    if(foundRangeDots == string::npos) {
      posS = atoi(regionString.substr(foundFirstColon + 1).c_str());
      chrE = chrS;
        posE = -1;
    } else {
      posS = atoi(regionString.substr(foundFirstColon + 1, 
        foundRangeDots - foundFirstColon - 1).c_str() );
      size_t foundSecondColon = regionString.find(":", foundRangeDots + 1);
       if(foundSecondColon == string::npos) {
        chrE = chrS;
        posE = atoi(regionString.substr(foundRangeDots + 2).c_str());
       } else {
        chrE = regionString.substr(foundRangeDots + 2, 
            foundSecondColon - (foundRangeDots + 2));
        posE = atoi(regionString.substr(foundSecondColon + 1).c_str());
      }
    }
  }

  const RefVector references = reader.GetReferenceData();

  int refIdS = reader.GetReferenceID(chrS);
  if(refIdS == (int)references.size()) return false;
  const RefData& refS = references.at(refIdS);
  if(posS > refS.RefLength) return false;

  int refIdE = reader.GetReferenceID(chrE);
  if(refIdE == (int)references.size()) return false;
  const RefData& refE = references.at(refIdE);
  if(posE > refE.RefLength) return false;

  if(posE == -1) posE = refE.RefLength;

  region.LeftRefID = refIdS;
  region.LeftPosition = posS;
  region.RightRefID = refIdE;
  region.RightPosition = posE;
  return true;
}
bool Utilities::CheckReadPair(const BamAlignment& al1, 
    const BamAlignment& al2, int32_t& is) {
  string name1 = al1.Name, name2 = al2.Name; 
  bool isfirst1 = al1.IsFirstMate(), isfirst2 = al2.IsFirstMate();
  bool ismapped1 = al1.IsMapped(), ismapped2 = al2.IsMapped();
  bool ismatemapped1 = al1.IsMateMapped(), ismatemapped2 = al2.IsMateMapped();
  int is1 = al1.InsertSize, is2 = al2.InsertSize;
  int pos1 = al1.Position, pos2 = al2.Position;
  int matepos1 = al1.MatePosition, matepos2 = al2.MatePosition;
  string rg1, rg2;
  al1.GetTag("RG", rg1);
  al2.GetTag("RG", rg2);

  if(!ismapped1 && pos1 != -1) pos1 = -1;
  if(!ismapped2 && pos2 != -1) pos2 = -1;

  if(rg1 != rg2) {
    cout << format("reads not belong to same RG: %s != %s\n") % rg1 % rg2;
    return false;
  }

  if(name1 != name2) {
    cout << format("read name not identical: %s != %s\n") % name1 % name2;
    return false;
  }

  if(ismapped1 != ismatemapped2 || ismapped2 != ismatemapped1) {
    cout << format("mapping tag inconsistency: %s\n %d-%d %d-%d\n") % name1 % ismapped1 % ismatemapped2 % ismapped2 % ismatemapped1;
    return false;
  }

  if( (isfirst1 && isfirst2) || (!isfirst1 && !isfirst2) ) {
    cout << format("reads not paired: %s\n") % name1;
    return false;
  }

  if(abs(is1) != abs(is2)) {
    cout << format("insert size not equal: %s [%d %d]\n") % name1 % is1 % is2;
    return false;
  }

  if(pos1 != matepos2 || pos2 != matepos1) {
    cout << format("positions not match: [%d - %d] : [%d - %d]\n") % pos1 % matepos1 % pos2 % matepos2;
    return false;
  }

  is = abs(is1);
  return true;
}
void Utilities::PrintReads(const string& name, const int& n_total, BamAlnVector& als) {
  cout << format("%s\t%d\n") % name % n_total;
  for (BamAlnVector::iterator it = als.begin(); it != als.end(); it++) {
    BamAlignment al(*it);
    string rg;
    al.GetTag("RG", rg);
    cout << format("\t%d[%s] (%d-%d) [%d:%d %d] %d\n") % al.IsFirstMate() % rg % al.IsMapped() % al.IsMateMapped() % al.RefID % al.Position % al.MatePosition % al.InsertSize;
  }
}
void Utilities::PrintReads(const string& name, const int& n_total, BamAlnPairVector& alpairs) {
  cout << format("%s\t%d\n") % name % n_total;
  for (BamAlnPairVector::iterator it = alpairs.begin(); it != alpairs.end(); it ++) {
    BamAlignment al1 = it->first, al2 = it->second;
    string rg1, rg2;
    al1.GetTag("RG", rg1);
    al2.GetTag("RG", rg2);
    cout << format("\t%d[%s] (%d-%d) [%d:%d %d] %d\n") % al1.IsFirstMate() % rg1 % al1.IsMapped() % al1.IsMateMapped() % al1.RefID % al1.Position % al1.MatePosition % al1.InsertSize;
    cout << format("\t%d[%s] (%d-%d) [%d:%d %d] %d\n") % al2.IsFirstMate() % rg2 % al2.IsMapped() % al2.IsMateMapped() % al2.RefID % al2.Position % al2.MatePosition % al2.InsertSize;
  }
}
void Utilities::PrintRead(const BamAlignment& al, const RefVector& refs) {
  string mate1 = al.IsFirstMate() ? "mate1" : "mate2";
  string qc = al.IsFailedQC() ? "failedQC" : "passedQC";
  string dup = al.IsDuplicate() ? "dup" : "!dup";
  string proper = al.IsProperPair() ? "proper" : "!proper";
  string mapped1 = al.IsMapped() ? "mapped" : "!mapped";
  string mapped2 = al.IsMateMapped() ? "mapped" : "!mapped";
  string strand1 = al.IsReverseStrand() ? "-" : "+";
  string strand2 = al.IsMateReverseStrand() ? "-" : "+";
  string chr1 = al.RefID == -1 ? "?" : refs[al.RefID].RefName;
  string chr2 = al.MateRefID == -1 ? "?" : refs[al.MateRefID].RefName;
  cout << format("%s\t%s\t%s\t%s\t%s\tIS=%d\n") % al.Name % mate1
    % qc % dup % proper % al.InsertSize;
  cout << format("  me\t%s\t%s:%d[%s] (%d)\n") % mapped1
    % chr1 % al.Position % strand1 % al.MapQuality;
  cout << format("  buddy\t%s\t%s:%d[%s]\n") % mapped2
    % chr2 % al.MatePosition % strand2;
}
bool isPair1 (const BamAlignment& a1, const BamAlignment& a2) {
  BamAlignment al1 = a1, al2 = a2;
  if((al1.IsFirstMate() && al2.IsFirstMate()) || (!al1.IsFirstMate() && !al2.IsFirstMate()))
    return false;
  if(!al1.IsMapped() && !al2.IsMapped()) {
    return true;
  } else if(al1.IsMapped() && al2.IsMapped()) {
    if(al1.RefID == al2.RefID && al1.MatePosition == al2.Position && al1.Position == al2.MatePosition && abs(al1.InsertSize) == abs(al2.InsertSize)) 
      return true;
    else
      return false;
  } else {
    return false;
  }
}
bool isPair2 (const BamAlignment& a1, const BamAlignment& a2) {
  BamAlignment al1 = a1, al2 = a2;
  if((al1.IsFirstMate() && al2.IsFirstMate()) || (!al1.IsFirstMate() && !al2.IsFirstMate()))
    return false;
  if(!a1.IsMapped() && a2.IsMapped()) {
    al1 = a2;
    al2 = a1;
  }
  if(!al1.IsMapped() && !al2.IsMapped()) {
    return true;
  } else if(al1.IsMapped() && !al2.IsMapped()) { 
    if(al1.RefID == al2.RefID && al1.Position == al2.MatePosition) 
      return true;
    else 
      return false;
  } else {
    return false;
  }
}
bool Utilities::PairReads(const string& name, const BamAlnVector& als_o, int& n_total, BamAlnPairVector& alpairs) {
  BamAlnVector als = als_o;
  if( (int)als.size() % 2 != 0 ) {
    cout << "mate1 and mate2 not half to half" << endl;
    return false;
  }
  n_total = (int)als.size() / 2;
  if(n_total == 1) {
    alpairs.push_back( BamAlnPair (als[0], als[1]) );
  } else {
    BamAlnVector::iterator it1, it2;
    for (int i=0; i<n_total; i++) {
      for(it1 = als.begin(); it1 != als.end(); it1 ++ )
        if( (*it1).IsMapped() ) break;
      if(it1 == als.end()) {
        cout << format("%d pairs for %s, yet none is mapped\n") % n_total % name;
        return false;
      }
      BamAlignment al1 = *it1;
      BamAlnVector als1;
      als1.push_back(al1);
      als.erase(it1);
      it2 = search(als.begin(), als.end(), als1.begin(), als1.begin()+1, isPair1);
      if(it2 == als.end()) it2 = search(als.begin(), als.end(), als1.begin(), als1.begin()+1, isPair2);
      if(it2 == als.end()) {
        cout << format("cannot find pair\n");
        return false;
      } else {
        BamAlignment al2 = *it2;
        als.erase(it2);
        alpairs.push_back( BamAlnPair (al1, al2) );
      }
    }
    if( !als.empty() || (int)alpairs.size() != n_total ) return false;
  }
  return true;
}
bool Utilities::FixReads(const string& name, BamAlnVector& als, int& type, int& n_total, vector<int32_t>& iss, BamAlnPairVector& alpairs) {
  vector<int> types;
  if(!alpairs.empty()) alpairs.clear();
  if(!iss.empty()) iss.clear();

  string rg;
  if((int)als.size() > 2) {
    vector<int> idxs_rm;
    als[0].GetTag("RG", rg);
    for(int i=1; i<(int)als.size(); i++) { 
      string rg2;
      als[i].GetTag("RG", rg2);
      if( rg != rg2 ) idxs_rm.push_back(i);
    }
    for(vector<int>::reverse_iterator rit = idxs_rm.rbegin(); 
      rit < idxs_rm.rend(); rit ++)
      als.erase(als.begin() + (*rit));
  }
  if( !this->PairReads(name, als, n_total, alpairs) ) return false;

  for(BamAlnPairVector::iterator it = alpairs.begin(); 
    it != alpairs.end(); it ++) {
    BamAlignment &al1 = it->first, &al2 = it->second;
    bool ismapped1(al1.IsMapped()), ismatemapped1(al1.IsMateMapped());
    bool ismapped2(al2.IsMapped()), ismatemapped2(al2.IsMateMapped());
    int type1 = (ismapped1 && ismatemapped1) ? 
      1 : (!ismapped1 && !ismatemapped1) ? 3 : 2;
    int type2 = (ismapped2 && ismatemapped2) ? 
      1 : (!ismapped2 && !ismatemapped2) ? 3 : 2;
    if(type1 == type2) {
      if(type1 == 1 && al1.InsertSize == al2.InsertSize)
        al2.InsertSize = - al1.InsertSize;
      types.push_back(type1);
    } else if ( type1==1 && type2==2 ) {
      types.push_back(2);
      if(al1.IsMateMapped()) al1.SetIsMateMapped(false);
      if(al1.MatePosition != -1) al1.MatePosition = -1;
      if(al1.InsertSize != 0) al1.InsertSize = 0;
    } else if ( type1==2 && type2==1 ) {
      types.push_back(2);
      if(al2.IsMateMapped()) al2.SetIsMateMapped(false);
      if(al2.MatePosition != -1) al2.MatePosition = -1;
      if(al2.InsertSize != 0) al2.InsertSize = 0;
    }
  }
  type = *min_element(types.begin(), types.end());

  //remove inconsistent/redundant alnPairs
  vector<int> idxs_rm;
  for(int i=0; i<n_total; i++) 
  if( (types[i] < type) || (i > 0 && type == 3) )
    idxs_rm.push_back(i);
  for(vector<int>::reverse_iterator rit = idxs_rm.rbegin(); 
    rit < idxs_rm.rend(); rit ++) {
    int idx = *rit;
    alpairs.erase(alpairs.begin() + idx);
    types.erase(types.begin() + idx);
  }
  assert(types.size() == alpairs.size());
  n_total = (int)alpairs.size();

  //checking paired read, store insert_size into iss
  for (BamAlnPairVector::iterator it = alpairs.begin(); 
    it != alpairs.end(); it ++) {
    int32_t is;
    if( !this->CheckReadPair(it->first, it->second, is) ) return false;
    iss.push_back(is);
  }
  return true;
}
void Utilities::Reverse(string& sequence) {
  reverse(sequence.begin(), sequence.end());
}
void Utilities::ReverseComplement(string& sequence) {
  size_t seqLen = sequence.length();
  for(size_t i = 0; i < seqLen; ++i) 
  sequence.replace(i, 1, 1, REVCOMP_LOOKUP[(int)sequence.at(i) - 65]);
  Reverse(sequence);
}
string Utilities::getCigarString(const vector<CigarOp> cigars) {
  string cigarStr = "";
  for(unsigned i = 0; i < cigars.size(); i ++) 
  cigarStr += str( format("%d%s") % cigars[i].Length % cigars[i].Type );
  return cigarStr;
}

bool BamReader2::GetReads(BamAlnVector& als, string& name, BamAlnVector& als2) {
  if(!als.empty()) als.clear();
  if(!als2.empty()) {
    assert ((int) als2.size() == 1);
    als.push_back(als2[0]);
    name = als2[0].Name;
    als2.clear();
  }
  BamAlignment al;
  bool tag_end = true, tag_dn = false;
  while( this->GetNextAlignment(al) ) {
    if(tag_end == true) tag_end = false;
    if(name.empty()) name = al.Name;
    if(name != al.Name) {
      als2.push_back(al);
      tag_dn = true;
      break;
    } else {
      als.push_back(al);
    }
  }
  if(tag_end == true) {
    return false;
  } else if(tag_dn == true) {
    assert((int)als2.size() == 1);
    return true;
  } else {
    assert(als2.empty());
    return true;
  }
}


