#include "StopwatchManager.h"
#include <iomanip>

using namespace std;

namespace AWI
{
   StopwatchManager::StopwatchManager()
   : watches(), output_buffer()
   {
   	// http://www.yolinux.com/TUTORIALS/LinuxTutorialMixingFortranAndC.html
   	setbuf(stdout,NULL); /* set output to unbuffered */
   }
   

   void StopwatchManager::start(std::string name)
   {
      // [] does not work to insert an element as we do not have a default constructor for Stopwatch
      // http://stackoverflow.com/a/695663
      watches.push_back(make_pair(name, Stopwatch(name)));
   }
   

   void StopwatchManager::stop(std::string name)
   {
      vector< std::pair<std::string, Stopwatch> >::iterator it;
      for(it = watches.begin(); it != watches.end(); it++)
      {
         if(name == it->first)
         {
            watches.erase(it);
            long elapsed_milliseconds = it->second.elapsed_milliseconds();
            StopwatchManager::print(name, elapsed_milliseconds, cout);
            StopwatchManager::print(name, elapsed_milliseconds, output_buffer);
            return;
         }
      }
   }

   
   void StopwatchManager::stop()
   {
      if(watches.rbegin() != watches.rend())
      {
         string lastname = watches.rbegin()->first;
         stop(lastname);
      }
   }
   
   
   void StopwatchManager::print(std::string name)
   {
      vector< std::pair<std::string, Stopwatch> >::iterator it;
      for(it = watches.begin(); it != watches.end(); it++)
      {
         if(name == it->first)
         {
            long elapsed_milliseconds = it->second.elapsed_milliseconds();
            StopwatchManager::print(name, elapsed_milliseconds, output_buffer);
            return;
         }
      }
   }

   
   void StopwatchManager::print(std::string name, long elapsed_milliseconds, std::ostream &outstream)
   {
      outstream<<setiosflags(ios::fixed)<<setprecision(3)<<"elapsed seconds for '"<<name<<"': "<<(elapsed_milliseconds/1000.0)<<endl;
   }
   
   
   void StopwatchManager::dump(std::ostream &outstream)
   {
      outstream << output_buffer.rdbuf();
   }

   
}