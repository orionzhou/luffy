#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <map>
#include <Sequence/RedundancyCom95.hpp>
using namespace std;
using namespace Sequence;
using std::vector;
using std::string;
using boost::format;

int main(int argc, char *argv[]) {
  const char chs = {"A", "T", "C", "G"};
  Sequence::RedundancyCom95 rc;

  for( unsigned i = 0; i < 4; ++i ) {
    for( unsigned j = 0; j < 4; ++j ) {
      for( unsigned k = 0; k < 4; ++k ) {
        string codon( { dna_alphabet[i], dna_alphabet[j], 
          dna_alphabet[k] } );
        cout << format("%s: %d %d %d   %d %d %d\n") % codon 
          % rc.FirstNon(codon) % rc.First2S(codon) + rc.First2V(codon) % -1
          % rc.ThirdNon(codon) % rc.Third2S(codon) + rc.Third2V(codon) % rc.ThirdFour(codon);
      }
    }
  }
}
