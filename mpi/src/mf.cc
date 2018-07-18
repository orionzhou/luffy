#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>
#include <iomanip>
#include <math.h>
#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <boost/program_options.hpp>
#include <boost/regex.hpp>
#include <boost/assign/std/vector.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include "mpi.h"
#include "task_division.h"

using namespace std;
using boost::format;
using namespace boost::assign;
namespace fs = boost::filesystem;
namespace po = boost::program_options;

struct FqRecord {
  string id;
  string prefix;
  bool gzipped;
};
typedef vector<FqRecord> FqRecordVector;
FqRecordVector get_fq_files ( string fi ) {
  FqRecordVector frv;
  ifstream fhi( fi.c_str() );
  if(!fhi.is_open()) {
    cout << format("cannot open file: %s") % fi;
    exit(1);
  } else {
    string line;
    while(fhi.good()) {
      getline(fhi, line);
      if(line.length() == 0) continue;
      vector<string> ss;
      boost::split(ss, line, boost::is_any_of("\t"));
      if(ss[0] == "id") continue;
      FqRecord fr;
      fr.id = ss[0];
      fr.prefix = ss[1];
      fr.gzipped = boost::lexical_cast<bool>(ss[2]);
      frv.push_back( fr );
    }
  }
  return frv;
}
void run_bwa (int task_beg, int task_end, FqRecordVector frv, const string& dirW) {
  string genomedb = "/scratch1/zhoup/Data/db/bwa/mt_35";
  for(int i = task_beg; i <= task_end; i ++) {
    FqRecord fr = frv[i];
    string f01_1 = dirW + "06_reads/" + fr.id + "_" + fr.prefix + ".1.fq";
    string f01_2 = dirW + "06_reads/" + fr.id + "_" + fr.prefix + ".2.fq";
    if(fr.gzipped) {
      f01_1 += ".gz";
      f01_2 += ".gz";
    }
    
    vector<string> cmds;
    string f02_1 = dirW + "11_pipe_bwa/03_bwa/" + fr.id + "_" + fr.prefix + ".1.sai";
    string f02_2 = dirW + "11_pipe_bwa/03_bwa/" + fr.id + "_" + fr.prefix + ".2.sai";
    string f03 = dirW + "11_pipe_bwa/03_bwa/" + fr.id + "_" + fr.prefix + ".bam";
    string f04 = dirW + "11_pipe_bwa/06_pos_sorted/" + fr.id + "_" + fr.prefix;
    cmds.push_back( str(format("bwa aln -n 8 -t 4 %s %s > %s") % genomedb % f01_1 % f02_1) );
    cmds.push_back( str(format("bwa aln -n 8 -t 4 %s %s > %s") % genomedb % f01_2 % f02_2) );
    cmds.push_back( str(format("bwa sampe %s %s %s %s %s | samtools view -Sb - > %s") % genomedb % f02_1 % f02_2 % f01_1 % f01_2 % f03) );
   
    for(vector<string>::iterator it = cmds.begin(); it != cmds.end(); it++) {
//      if(it - cmds.begin() < 3) continue;
      cout << it->c_str() << endl;
      system(it->c_str());
    }
  }
}

int main(int argc, char *argv[]) {
  int row_beg, row_end;
  po::options_description cmdOpts("Allowed options");
  cmdOpts.add_options()
    ("help,h", "produce help message")
    ("row_beg,b", po::value<int>(&row_beg), "begin row")
    ("row_end,e", po::value<int>(&row_end), "end   row")
  ;
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, cmdOpts), vm);
  po::notify(vm);
  if(vm.count("help") || !vm.count("row_beg") || !vm.count("row_end")) {
    cout << cmdOpts << endl;
    return 1;
  }

  int np, myrank;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &np);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Status status;

  string dirW = "/scratch1/zhoup/Data/repo/mt_35/";
  string fq_info = dirW + "06_reads/01.tbl";
  FqRecordVector frv = get_fq_files(fq_info);
  int nt = row_end - row_beg + 1;

  if (myrank == 0) {
    TaskDivVector tdv = task_divide ( nt, 1, np-1 );
    cout << format("%d procs, %d tasks\n") % (np-1) % nt;
    for(TaskDivVector::iterator it = tdv.begin(); it < tdv.end(); it ++) {
      int proc = it->proc, task_beg = it->task_beg, task_end = it->task_end;
      int task[2] = {task_beg, task_end};
      cout << format("proc[%d]: sending task[%d - %d] to proc[%d]\n") % myrank % task_beg % task_end % proc;
      MPI_Send(task, 2, MPI_INT, proc, 1, MPI_COMM_WORLD);
    }

    clock_t time1 = clock();
    for(TaskDivVector::iterator it = tdv.begin(); it < tdv.end(); it ++) {
      int proc = it->proc, task_beg = it->task_beg, task_end = it->task_end;
      int task[2];
      MPI_Recv(&task, 2, MPI_INT, proc, 1, MPI_COMM_WORLD, &status);
      cout << format("proc[%d]: receiving task[%d - %d] from proc[%d]\n") % myrank % task[0] % task[1] % proc;
    }
    cout << right << setw(60) << format("(running time: %.01f minutes)\n") % ( (double)(clock() - time1) / ((double)CLOCKS_PER_SEC * 60) );
  }  else {
    int task[2];
    MPI_Recv(&task, 2, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);
    cout << format("proc[%d]: receiving task[%d - %d] from proc[0]\n") % myrank % task[0] % task[1];
    int task_beg = task[0] + row_beg - 1;
    int task_end = task[1] + row_beg - 1;
    run_bwa(task_beg, task_end, frv, dirW);
    MPI_Send(&task, 2, MPI_INT, 0, 1, MPI_COMM_WORLD);
  }
  MPI_Finalize();
  return EXIT_SUCCESS;
}



