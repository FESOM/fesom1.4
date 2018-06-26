#ifndef OutputSchedulerFC_B84C3704_CD0B_4AF6_AA81_B70EBE5C0399
#define OutputSchedulerFC_B84C3704_CD0B_4AF6_AA81_B70EBE5C0399

#include "OutputSchedulerFCMacros.h"
extern "C" {
   void should_output(const char* name, int *istep, double *dsec, int *day_step, int *month_step, int *year, bool *is_last_day_in_month, bool *should_write_output);
   void schedule_info(const char* name, char* unit, int* first, int* rate);
   void print_schedules();
   void print_default_schedules();
}

#endif
