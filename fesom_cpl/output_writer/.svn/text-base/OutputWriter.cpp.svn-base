#include "OutputWriter.h"
#include "OutputWriterUtils.h"
#include <iomanip>
#include <sstream>
#include <cstdlib>

using namespace std;
using namespace netCDF;
using namespace netCDF::exceptions;

namespace AWI
{
   OutputWriter::OutputWriter(const std::string filepath_, const std::string var_name_, const int var_size, const std::string var_description, const std::string var_unit, const std::string dim_name, const std::string global_attribute_name, std::string global_attribute_txt, std::string units_txt, std::string calendar_txt)
   : filepath(filepath_), var_name(var_name_), save_count(0)
   {
      // cache the netcdf file handle to try to trick the netcdf library to allow access via independent threads (this seems to require a thread safe hdf5 library linked to netcdf if Ncfile::nc4)
      netcdf_file = std::unique_ptr<netCDF::NcFile>(new NcFile(filepath, NcFile::replace, NcFile::classic)); // as we currently use no compression: use classic format to not rely on having a thread-safe hdf5
      OutputWriter::create_output_file(*netcdf_file, var_name, dim_name, var_size, var_description, var_unit, global_attribute_name, global_attribute_txt, units_txt, calendar_txt);

      //netcdf_file->sync();
      int retCode = nc_sync(netcdf_file->getId());
      if(retCode!=NC_NOERR)
      {
         std::stringstream exceptionMessage;
         exceptionMessage << __FILE__ << ":" << __LINE__ << " netcdf error";
         throw std::runtime_error(exceptionMessage.str());
      }
   }
   
   
   // fesom data stored in a 1D netcdf array, suitable for fesom 2D output and sparse/compact 3D output with only wet nodes getting stored
   void OutputWriter::appendOutput(const int var_size, const float *var_data, const double ysec)
   {
      vector<size_t> index(1, save_count);
      NcVar time_var = netcdf_file->getVar("time");
      time_var.putVar(index, ysec);
      
      vector<size_t> start;
      start.push_back(save_count);
      start.push_back(0);
      vector<size_t> count;
      count.push_back(1);
      count.push_back(var_size);
      NcVar var = netcdf_file->getVar(var_name);
      var.putVar(start, count, &var_data[0]);

      //netcdf_file->sync();
      int retCode = nc_sync(netcdf_file->getId());
      if(retCode!=NC_NOERR)
      {
         std::stringstream exceptionMessage;
         exceptionMessage << __FILE__ << ":" << __LINE__ << " netcdf error";
         throw std::runtime_error(exceptionMessage.str());
      }

      save_count ++;
   }
   
   
   void OutputWriter::create_output_file(netCDF::NcFile &ncf, std::string varname, std::string dimname, int dimsize, std::string description_txt, std::string unit_txt, std::string global_attribute_name, std::string global_attribute_txt, std::string units_txt, std::string calendar_txt)
   {
      NcDim dim_rec = OutputWriterUtils::create_shared_header(ncf, units_txt, calendar_txt);
      
      NcDim dimid_2d = ncf.addDim(dimname, dimsize);
      vector<NcDim> dim_2d;
      dim_2d.push_back(dim_rec);
      dim_2d.push_back(dimid_2d);
      NcVar var = ncf.addVar(varname, ncFloat, dim_2d);
      
      var.putAtt("description", description_txt);
      var.putAtt("units", unit_txt);
      var.putAtt("grid_type", "unstructured");
      var.putAtt("_FillValue", var.getType(), 1.e+30f);
     // could not figure out how to do this via the C++ API (note: NcItem::idGlobal is NC_GLOBAL), so use plain C API
      nc_put_att_text(ncf.getId(), NC_GLOBAL, global_attribute_name.c_str(), global_attribute_txt.size(), global_attribute_txt.c_str());

      //      netcdf_file->enddef();
      int retCode = nc_enddef(ncf.getId());
      if(retCode!=NC_NOERR)
      {
         std::stringstream exceptionMessage;
         exceptionMessage << __FILE__ << ":" << __LINE__ << " netcdf error";
         throw std::runtime_error(exceptionMessage.str());
      }
   }

}
