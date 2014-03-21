#include "timer.h"
#include <iostream>
#include <cstdio>

int Timer::print(int rank) {
  if (rank < 0) {
    printf("Time for %-20s%10.2f s\n", (name + ":").c_str(), data);
  } else {
    printf("Processor %2d Time for %-20s%10.2f s\n", rank, (name + ":").c_str(), data);    
  }
  return 0;
}

double Timer::time() {
  return data;
}

int Timer::start() {
  gettimeofday(&t_start, NULL);
  return 0;
}

int Timer::pause() {
  timeval t_end;
  gettimeofday(&t_end, NULL);
  data += (double)(t_end.tv_sec-t_start.tv_sec)
      +((double)(t_end.tv_usec-t_start.tv_usec))/1000000.;
  return 0;
}

