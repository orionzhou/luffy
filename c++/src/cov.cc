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

int main( int argc, char* argv[] ) {
  string fi, fo;
  clock_t time1 = clock();
  po::options_description cmdOpts("Allowed options");
  cmdOpts.add_options()
    ("help,h", "produce help message")
    ("in,i", po::value<string>(&fi)->default_value("/tmp"), "input file")
    ("out,o", po::value<string>(&fo), "output file")
  ;
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, cmdOpts), vm);
  po::notify(vm);
  if(vm.count("help") || !vm.count("in") || !vm.count("out")) {
    cout << cmdOpts << endl;
    return 1;
  }

  string cmd = "ls -l";
  int MAX_BUFSIZE = 25535;
  char buffer[MAX_BUFSIZE];
  FILE *stream = popen(cmd.c_str(), "r");
  while( fgets(buffer, MAX_BUFSIZE, stream) != NULL ) {
    string line(buffer);
    if(line.substr(line.length()-1, 1) == "\n") line.erase(line.length()-1, 1);
    if(line.length() == 0) continue;
    vector<string> ss;
    boost::split(ss, line, boost::is_any_of(" "));
    cout << ss[ss.size()-2] << "\t" << ss[ss.size()-1] << endl;
  }
  pclose(stream);

  cout << right << setw(60) << format("(running time: %.01f minutes)\n") % ( (double)(clock() - time1) / ((double)CLOCKS_PER_SEC * 60) );
  return EXIT_SUCCESS;
}
