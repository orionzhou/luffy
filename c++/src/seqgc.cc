#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include <time.h>
#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>
using namespace std;
using boost::format;
namespace fs = boost::filesystem;
namespace po = boost::program_options;

float get_gc( string seq ) {
  int gcCount = 0;
  int totalBp = 0;
  for(string::const_iterator pos = seq.begin(); pos != seq.end(); ++ pos) {
    switch (toupper(*pos)) {
      case 'A' :
        totalBp ++;
        break;
      case 'C' :
        gcCount ++;
        totalBp ++;
        break;
      case 'G' :
        gcCount ++;
        totalBp ++;
        break;
      case 'T' :
        totalBp ++;
        break;
      default:
        break;
    }
  }
  return totalBp < 20 ? -1 : (float)gcCount / (float)totalBp;
}

int main( int argc, char* argv[] ) {
  string fi, fo;
  clock_t time1 = clock();
  po::options_description cmdOpts("Allowed options");
  cmdOpts.add_options()
    ("help,h", "produce help message")
    ("in,i", po::value<string>(&fi), "input file (fasta)")
    ("out,o", po::value<string>(&fo), "output file (BED)")
  ;
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, cmdOpts), vm);
  po::notify(vm);
  if(vm.count("help") || !vm.count("in") || !vm.count("out")) {
    cout << cmdOpts << endl;
    return 1;
  }

  ifstream fhi(fi.c_str());
  if(!fhi.is_open()) {
    cout << format("cannot open file: %s") % fi;
    return false;
  }
  ofstream fho(fo.c_str());
  if(!fho.is_open()) {
    cout << format("cannot open file: %s") % fo;
    return false;
  }

  string line, name, seq;
  while(fhi.good()) {
    getline(fhi, line);
    if(line.empty() || line[0] == '>') {
      if(!name.empty()) {
        fho << format("%s\t%.3f\n") % name % get_gc(seq);
        name.clear();
      } 
      if(!line.empty()) {
        name = line.substr(1);
      }
      seq.clear();
    } else if(!name.empty()) {
      if(line.find(' ') != string::npos) {
        name.clear();
        seq.clear();
      } else {
        seq += line;
      }
    }
  }

  if(!name.empty()) {
    fho << format("%s\t%.3f\n") % name % get_gc(seq);
  }
  cout << right << setw(60) << format("(running time: %.01f minutes)\n") % ( (double)(clock() - time1) / ((double)CLOCKS_PER_SEC * 60) );
  return EXIT_SUCCESS;
}


