#include "OutputWriterFC.h"
#include "OutputWriterManager.h"
#include <string>

static AWI::OutputWriterManager writerManager;


void add_ocean_levels_wetrange(const int *level, const int *index, const int *size)
{
   writerManager.addOceanLevelsWetrange(*level-1, *index-1, *size);
}

void add_ocean_depth(const double *depth)
{
   writerManager.addOceanDepth(*depth);
}


void set_output_file(const char* filepath, const char* var_name, const int *var_size, const char* var_description, const char* var_unit, const char* dim_name, const char* global_attribute_name, const char* global_attribute_txt, const char* units_txt, const char* calendar_txt)
{
   std::string filepath_(filepath);
   std::string var_name_(var_name);
   std::string var_description_(var_description);
   std::string var_unit_(var_unit);
   std::string dim_name_(dim_name);
   std::string global_attribute_name_(global_attribute_name);
   std::string global_attribute_txt_(global_attribute_txt);
   std::string units_txt_(units_txt);
   std::string calendar_txt_(calendar_txt);

   writerManager.setOutputFile(filepath_, var_name_, *var_size, var_description_, var_unit_, dim_name_, global_attribute_name_, global_attribute_txt_, units_txt_, calendar_txt_);
}


void set_ocean_levels_output_file(const char* filepath, const char* var_name, const int *nod2Dcount, const int *levelcount, const char* var_description, const char* var_unit, const char* dim_name, const char* global_attribute_name, const char* global_attribute_txt, const char* units_txt, const char* calendar_txt)
{
   std::string filepath_(filepath);
   std::string var_name_(var_name);
   std::string var_description_(var_description);
   std::string var_unit_(var_unit);
   std::string dim_name_(dim_name);
   std::string global_attribute_name_(global_attribute_name);
   std::string global_attribute_txt_(global_attribute_txt);
   std::string units_txt_(units_txt);
   std::string calendar_txt_(calendar_txt);

   writerManager.setOceanLevelsOutputFile(filepath_, var_name_, *nod2Dcount, *levelcount, var_description_, var_unit_, dim_name_, global_attribute_name_, global_attribute_txt_, units_txt_, calendar_txt_);
}


void begin_append_output(const char* var_name, const int *var_size, const float *var_data, const double *ysec)
{
   std::string var_name_(var_name);
   writerManager.beginAppendOutput(var_name_, *var_size, var_data, *ysec);
}


void end_append_output(const char* var_name)
{
   std::string var_name_(var_name);
   writerManager.endAppendOutput(var_name_);
}
