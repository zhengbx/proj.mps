#ifndef TIME
#define TIME

#include "sys/time.h"
#include <string>

class Timer {
  double data;
  std::string name;
  timeval t_start;
public:
  Timer(std::string _name = ""): data(0), name(_name) {};
  int print(int rank = -1);
  int start();
  int pause();
  double time();
};

#endif
