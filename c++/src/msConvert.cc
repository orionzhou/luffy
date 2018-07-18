#include <fstream>
#include <iostream>
#include <stdio.h>
#include <boost/lexical_cast.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <vector>

#include <Sequence/SimParams.hpp>
#include <Sequence/SimData.hpp>
#include <Sequence/PolySIM.hpp>

using namespace Sequence;
using namespace std;
namespace fs = boost::filesystem;

int main(int argc, char *argv[]) {
  fs::path dirW(argv[0]);
  fs::path f01 = dirW / "01.txt";
  std::cout << "opening " << f01.string().c_str() << " ...\n";
  if ( !exists( f01 ) ) {
    std::cout << "  Not found!\n";
    return 1;
  } 
  fs::path f02 = dirW / "02.txt";
  fs::fstream fh1( f01 );

  SimParams p;
  fh1 >> p;
  printf ("Reps: %3d  Individuals: %5d\n", p.runs(), p.totsam());

  SimData d(p.totsam());
//  ios_base::sync_with_stdio(true);
  for( unsigned i=1; i<=p.runs(); i++ ) {
    d.read( fh1 );
    PolySIM P(&d);
    printf ("ss: %3d  ThetaW:%5.02f  ThetaPi:%5.02f\n", P.NumPoly(), P.ThetaW(), P.ThetaPi());
  }
  fh1.close();
}

