#ifndef StopwatchManager_921DC723_33FF_4121_BFD8_5A0444583C25
#define StopwatchManager_921DC723_33FF_4121_BFD8_5A0444583C25

#include <iostream>
#include <string>
#include <vector>
#include <utility>
#include <sstream>

#include "Stopwatch.h"



namespace AWI
{
   class Stopwatch;
   
   class StopwatchManager
   {
   public:      
      StopwatchManager();
      void start(std::string name);
      void stop(std::string name);
      void stop();
      void print(std::string name);
      void dump(std::ostream &outstream);
      
   private:
      static void print(std::string name, long elapsed_milliseconds, std::ostream &outstream);

   private:
      std::vector< std::pair<std::string, Stopwatch> > watches;
      std::stringstream output_buffer;
   };
   
}

#endif
