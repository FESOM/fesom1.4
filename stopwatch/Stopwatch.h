#ifndef Stopwatch_41D20995_1B6F_4712_A0EB_48252AF2A3BD
#define Stopwatch_41D20995_1B6F_4712_A0EB_48252AF2A3BD

#include <iostream>
#include <string>
#include <sys/time.h>
#include <stdio.h>

namespace AWI
{   
   class Stopwatch
   {
   public:      
      Stopwatch(std::string name);
      
      long elapsed_milliseconds();

   private:
      std::string name;
      struct timeval start_time;
   };
   
}

#endif
