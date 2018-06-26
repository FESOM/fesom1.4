#include "OutputSchedulerFC.h"
#include "OutputSchedulerManager.h"
#include <stdexcept>
#include <cstring>

static AWI::OutputSchedulerManager schedulerManager;

void should_output(const char* name, int *istep, double *dsec, int *day_step, int *month_step, int *year, bool *is_last_day_in_month, bool *should_write_output)
{
   std::string n(name);
   // with ifort, we get a 0 for false and 255 for true from fortran here (for a one byte fortran bool, i.e. logical*1)
   // make sure we have a real bool, i.e. we can do val==true (which would evaluate to false if val is 255)
   bool is_last_day_in_month_ = *is_last_day_in_month != false;
   *should_write_output = schedulerManager.shouldOutput(n, *istep, *dsec, *day_step, *month_step, *year, is_last_day_in_month_);
}

void schedule_info(const char* name, char* unit_chars, int* first, int* rate)
{
   std::string name_(name);
   std::string unit_string;
   unsigned int first_, rate_;
   schedulerManager.scheduleInfo(name_, unit_string, first_, rate_);
   
   // copy our string to the given character array
   
   size_t char_array_size = sizeof(unit_chars) / sizeof(unit_chars[0]);
   if(unit_string.size() > char_array_size)
   {
      std::stringstream exceptionMessage;
      exceptionMessage << __FILE__ << ":" << __LINE__ <<" given char array size is too small: "<<char_array_size<<" required: "<<unit_string.size();
      throw std::runtime_error(exceptionMessage.str());
   }
   strcpy(unit_chars, unit_string.c_str());
   for(int i = 0; i < sizeof(unit_chars) / sizeof(unit_chars[0]); i++)
   // fill any empty items of the char array with the NUL symbol
   for(int i = 0; i < char_array_size-unit_string.size(); i++)
      unit_chars[unit_string.size()+i] = (char)32;
   
   *first = (int)first_;
   *rate = (int)rate_;
}

void print_schedules()
{
   schedulerManager.print();
}

void print_default_schedules()
{
   schedulerManager.printDefaultSchedules();
}
