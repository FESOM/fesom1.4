#include "Stopwatch.h"

using namespace std;

namespace AWI
{
   Stopwatch::Stopwatch(std::string name_)
   : name(name_)
   {
      gettimeofday(&start_time, NULL);
   }
   
   long Stopwatch::elapsed_milliseconds()
   {
      struct timeval end_time;
      gettimeofday(&end_time, NULL);
      long elapsed_seconds  = end_time.tv_sec  - start_time.tv_sec;
      long elapsed_useconds = end_time.tv_usec - start_time.tv_usec;
      long elapsed_milliseconds = ((elapsed_seconds) * 1000 + elapsed_useconds/1000.0) + 0.5;
      return elapsed_milliseconds;
   }
}