# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <ctime>
# include <cmath>
using namespace std;

struct TaskDiv {
  int proc;
  int task_beg;
  int task_end;
  int task_number;
};
typedef vector<TaskDiv> TaskDivVector;

int i4_div_rounded ( int a, int b );
int i4_huge ( void );
int i4_sign ( int i );
TaskDivVector task_divide ( int task_number, int proc_first, int proc_last );
void timestamp ( void );

int i4_div_rounded ( int a, int b ) {
  int a_abs, b_abs, c, c_abs, c_s;
  if ( a == 0 ) {
    c_abs = i4_huge ( );
    c_s = i4_sign ( b );
  } else {
    a_abs = abs ( a );
    b_abs = abs ( b );
    c_s = i4_sign ( a ) * i4_sign ( b );
    c_abs = a_abs / b_abs;
    if ( ( 2 * c_abs + 1 ) * b_abs < 2 * a_abs ) {
      c_abs = c_abs + 1;
    }
  }
  c = c_s * c_abs;
  return c;
}
int i4_sign ( int i ){
  int value;
  if ( i < 0 ) {
    value = -1;
  } else {
    value = 1;
  }
  return value;
}
int i4_huge ( void ) {
  return 2147483647;
}

TaskDivVector task_divide ( int task_number, int proc_first, int proc_last ) {
  int i_hi, i_lo, proc, proc_number, proc_remain, task_proc, task_remain;
  proc_number = proc_last + 1 - proc_first;

  TaskDivVector tdv;
  i_hi = 0;
  task_remain = task_number;
  proc_remain = proc_number;

  for ( proc = proc_first; proc <= proc_last; proc++ ) {
    task_proc = i4_div_rounded ( task_remain, proc_remain );

    proc_remain = proc_remain - 1;
    task_remain = task_remain - task_proc;

    i_lo = i_hi + 1;
    i_hi = i_hi + task_proc;
    
    TaskDiv td;
    td.proc = proc;
    td.task_beg = i_lo;
    td.task_end = i_hi;
    td.task_number = task_proc;
    tdv.push_back( td );
  }
  return tdv;
}

void timestamp ( void ) {
# define TIME_SIZE 40
  static char time_buffer[TIME_SIZE];
  const struct tm *tm;
  size_t len;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );
  len = strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );
  cout << time_buffer << "\n";
  return;
# undef TIME_SIZE
}
