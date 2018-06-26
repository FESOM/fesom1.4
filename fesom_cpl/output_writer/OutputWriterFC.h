#ifndef OutputWriter_DC5415FF_54AC_416A_89A7_F6B8C63C9C10
#define OutputWriter_DC5415FF_54AC_416A_89A7_F6B8C63C9C10

#include "OutputWriterFCMacros.h"
extern "C" {
   void add_ocean_levels_wetrange(const int *level, const int *index, const int *size);
   void add_ocean_depth(const double *depth);
   void set_output_file(const char* filepath, const char* var_name, const int *var_size, const char* var_description, const char* var_unit, const char* dim_name, const char* global_attribute_name, const char* global_attribute_txt, const char* units_txt, const char* calendar_txt);
   void set_ocean_levels_output_file(const char* filepath, const char* var_name, const int *nod2Dcount, const int *levelcount, const char* var_description, const char* var_unit, const char* dim_name, const char* global_attribute_name, const char* global_attribute_txt, const char* units_txt, const char* calendar_txt);
   void begin_append_output(const char* var_name, const int *var_size, const float *var_data, const double *ysec);
   void end_append_output(const char* var_name);
}

#endif
