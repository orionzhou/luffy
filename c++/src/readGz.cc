#include <cstdlib>
#include <fstream>
#include <iostream>
#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <boost/program_options.hpp>
#include <boost/regex.hpp>
#include <boost/assign/std/vector.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/progress.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
using namespace std;
using boost::format;
namespace fs = boost::filesystem;
namespace io = boost::iostreams;
namespace po = boost::program_options;

int main(int argc, char *argv[]) {
  string fi, fo;
  po::options_description cmdOpts("Allowed options");
  cmdOpts.add_options()
    ("help,h", "produce help message")
    ("in,i", po::value<string>(&fi), "input (simpleSNP) file")
    ("out,o", po::value<string>(&fo), "output file")
  ;
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, cmdOpts), vm);
  po::notify(vm);
  if(vm.count("help") || !vm.count("in") || !vm.count("out")) {
    cout << cmdOpts << endl;
    return 1;
  }

  ifstream fh( fi.c_str(), ios_base::in | ios_base::binary );
  io::filtering_stream<io::input> in;
  in.push(io::gzip_decompressor());
  in.push(fh);
  
//  streamsize ss;
  string line;
  uint32_t cnt = 0;
  while ( getline(in, line) ) {
    cnt ++;
    if(cnt > 10) break;
    cout << format("%4d: %s\n") % cnt % line;
  }
}

