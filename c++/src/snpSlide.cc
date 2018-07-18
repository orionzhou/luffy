#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include <time.h>
#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/program_options.hpp>
#include <Sequence/PolySites.hpp>
#include <Sequence/PolyTableSlice.hpp>
#include "snp.h" 
using namespace std;
using boost::format;
using namespace Sequence;
namespace po = boost::program_options;

int main(int argc, char *argv[]) {
  string fi, fo;
  clock_t time1 = clock();
  po::options_description cmdOpts("Allowed options");
  cmdOpts.add_options()
    ("help,h", "produce help message")
    ("in,i", po::value<string>(&fi)->implicit_value(""), "input (SimpleSNP format)")
    ("out,o", po::value<string>(&fo)->implicit_value(""), "output")
  ;

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, cmdOpts), vm);
  po::notify(vm);

  if(vm.count("help") || !vm.count("in") || !vm.count("out")) {
    cout << cmdOpts << endl;
    return 1;
  }

  istream* in = (!vm.count("in") || fi == "" || fi == "-" || fi == "stdin") ?
    &cin : new ifstream( fi.c_str() );
  ostream* out = (!vm.count("out") || fo == "" || fo == "-" || fo == "stdout") ?
    &cout : new ofstream( fo.c_str() );

  SNP snp;
  snp.read(*in);
  
  cout << format("%3d inds, %6d loci\n") % snp.get_nind() % snp.get_npos();
  PolySites psi(snp.sbegin(), snp.send());
//  PolySNP psp(&psi);
//  cout << "ThetaPi: " << psp.ThetaPi() << endl;
//  cout << "ThetaW: " << psp.ThetaW() << endl;
//  cout << "Tajima's D: " << psp.TajimasD() << endl;
  double chrLen = *(psi.pbegin()+snp.get_npos()-1) + 1;
  unsigned winSize = 100000; 
  unsigned winStep = 100000; 
  PolyTableSlice<PolySites> pts(psi.sbegin(), psi.send(), winSize, winStep, chrLen);

  for (unsigned i = 0; i < pts.size(); i++) {
    PolySites psi1(pts[i]);
    //PolySNP psp1(&psi);
    //double thetaW = winCalc.ThetaW() / winSize;
    if(psi1.numsites() > 0) {
      *out << format("%.0f\t%.0f") % *(psi1.pbegin()) % *(psi1.pbegin()+psi1.numsites()-1) << endl;
    } else {
      *out << format("\t") << endl;
    }
  }
  cout << format("(running time: %.01f minutes)\n") % ( (double)(clock() - time1) / ((double)CLOCKS_PER_SEC * 60) );
	return EXIT_SUCCESS;
}


