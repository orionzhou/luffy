#include <stdio.h>
#include <iostream>
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

template< class T >
struct next {
  T seed;
  next( T seed ) : seed(seed) { }
  T operator()() {
    return seed++;
  }
};

int main (int argc, char *argv[]) {
  int nt = 10, np, myrank;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &np);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Status status;

  if (myrank == 0) {
    vector<int> tasks;
    vector<int> tasks_done;
    push_back(tasks).repeat_fun(nt, next<int>(1));
    
    TaskDivVector tdv = task_divide ( nt, 1, np-1 );
    for(TaskDivVector::iterator it = tdv.begin(); it < tdv.end(); it ++) {
      int proc = it->proc, task_beg = it->task_beg, task_end = it->task_end;
      
      int task[2] = {task_beg, task_end};
      cout << format("proc[%d]: sending task[%d - %d] to proc[%d]\n") % myrank % task_beg % task_end % proc;
      MPI_Send(task, 2, MPI_INT, proc, 1, MPI_COMM_WORLD);
    }

    for(TaskDivVector::iterator it = tdv.begin(); it < tdv.end(); it ++) {
      int proc = it->proc, task_beg = it->task_beg, task_end = it->task_end;

      int task[2];
      MPI_Recv(&task, 2, MPI_INT, proc, 1, MPI_COMM_WORLD, &status);
      cout << format("proc[%d]: receiving task[%d - %d] from proc[%d]\n") % myrank % task[0] % task[1] % proc;
    }
  }  else {
    int task[2];
    MPI_Recv(&task, 2, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);
    cout << format("proc[%d]: receiving task[%d - %d] from proc[0]\n") % myrank % task[0] % task[1];
    MPI_Send(&task, 2, MPI_INT, 0, 1, MPI_COMM_WORLD);
  }

  MPI_Finalize();
  return EXIT_SUCCESS;
}



