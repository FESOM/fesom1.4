#ifndef OutputWriter_A9246E43_0DB9_47F6_8A23_A7083DAB96BF
#define OutputWriter_A9246E43_0DB9_47F6_8A23_A7083DAB96BF

#include <iostream>
#include <string>
#include <vector>
#include <utility>
#include <sstream>
#include <map>
#include "OutputWriterThing.h"
#include <memory>
#include <netcdf>

namespace netCDF
{
   class NcDim;
}

namespace AWI
{
   
   class OutputWriter : public OutputWriterThing
   {
   public:      
      OutputWriter(const std::string filepath, const std::string var_name, const int var_size, const std::string var_description, const std::string var_unit, const std::string dim_name, const std::string global_attribute_name, std::string global_attribute_txt, std::string units_txt, std::string calendar_txt);
      void appendOutput(const int var_size, const float *var_data, const double ysec);
      OutputWriter(const OutputWriter&) = delete;
      void operator=(const OutputWriter&) = delete;
      
   private:
      static void create_output_file(netCDF::NcFile &ncf, std::string varname, std::string dimname, int dimsize, std::string description_txt, std::string unit_txt, std::string global_attribute_name, std::string global_attribute_txt, std::string units_txt, std::string calendar_txt);
   
   private:
      const std::string filepath;
      const std::string var_name;
      int save_count;
      std::unique_ptr<netCDF::NcFile> netcdf_file;
   };
   
}

#endif
