#ifndef OceanLevelsOutputWriter_DC84C5BC_F849_4220_BD7C_F45B54D720B1
#define OceanLevelsOutputWriter_DC84C5BC_F849_4220_BD7C_F45B54D720B1

#include <string>
#include <vector>
#include "OutputWriterThing.h"
#include <memory>
#include <netcdf>

namespace netCDF
{
   class NcDim;
   class NcVar;
}

namespace AWI
{
   
   class OceanLevelsOutputWriter : public OutputWriterThing
   {
   public:      
      OceanLevelsOutputWriter(const std::string filepath, const std::string var_name, const int nod2Dcount, const int levelcount, const std::vector<double> * const depths, const std::vector< std::vector< std::pair<size_t,size_t> > > * const wetranges, const std::string var_description, const std::string var_unit, const std::string dim_name, const std::string global_attribute_name, std::string global_attribute_txt, std::string units_txt, std::string calendar_txt);
      void appendOutput(const int size, const float *data, const double ysec);
      OceanLevelsOutputWriter(const OceanLevelsOutputWriter&) = delete;
      void operator=(const OceanLevelsOutputWriter&) = delete;

   private:
      void create_output_file(netCDF::NcFile &ncf, std::string varname, std::string dimname, int dimsize_x, int dimsize_y, std::string description_txt, std::string unit_txt, std::string global_attribute_name, std::string global_attribute_txt, std::string units_txt, std::string calendar_txt);
      static void writeFragmentToNetcdf(netCDF::NcVar &var, const size_t save_count, const float *var_data, const int levelindex, const int columnindex, const int columncount, const int dataindex);
   
   private:
      const std::string filepath;
      const std::string var_name;
      const std::vector<double> *depths;
      const std::vector< std::vector< std::pair<size_t,size_t> > > *wetranges;
      std::unique_ptr<netCDF::NcFile> netcdf_file;
   };
   
}

#endif
