#ifndef OutputSchedulerManager_0F498594_2624_4073_A395_BE6167ADDE05
#define OutputSchedulerManager_0F498594_2624_4073_A395_BE6167ADDE05

#include <iostream>
#include <string>
#include <vector>
#include <utility>
#include <sstream>
#include <map>
#include "OutputScheduler.h"

namespace YAML
{
   class Node;
}

namespace AWI
{
   class OutputScheduler;
   
   class OutputSchedulerManager
   {
   public:      
      OutputSchedulerManager();
      /*
       returns true if the given time is past the start time for the given schedule and the given time is a multiple of the step from the given schedule
       */
      bool shouldOutput(std::string name, int istep, double dsec, int day_step, int month_step, int year, bool is_last_day_in_month);
      void scheduleInfo(std::string name, std::string& unit, unsigned int& first, unsigned int& rate);
      void print();
      void printDefaultSchedules();
      
   private:
      static OutputScheduler parseSchedule(const YAML::Node& schedule);
      static std::string defaultSchedules();

   private:
      std::map<std::string, OutputScheduler> schedules;
   };
   
}

#endif
