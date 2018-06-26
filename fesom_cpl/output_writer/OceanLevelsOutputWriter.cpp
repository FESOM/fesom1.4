#include "OceanLevelsOutputWriter.h"
#include "OutputWriterUtils.h"
#include <iomanip>
#include <sstream>
#include <cstdlib>
#include <cassert>

using namespace std;
using namespace netCDF;
using namespace netCDF::exceptions;

namespace AWI
{
   OceanLevelsOutputWriter::OceanLevelsOutputWriter(const std::string filepath_, const std::string var_name_, const int nod2Dcount, const int levelcount, const std::vector<double> * const depths_, const std::vector< std::vector< std::pair<size_t,size_t> > > * const wetranges_, const std::string var_description, const std::string var_unit, const std::string dim_name, const std::string global_attribute_name, std::string global_attribute_txt, std::string units_txt, std::string calendar_txt)
   : filepath(filepath_), var_name(var_name_), depths(depths_), wetranges(wetranges_)
   {
      if(depths->size() != levelcount || wetranges->size() != levelcount)
      {
         std::stringstream exceptionMessage;
         exceptionMessage << __FILE__ << ":" << __LINE__ <<" different number of ocean levels: levelcount: "<<levelcount<<", depths->size(): "<<depths->size()<<", wetranges->size(): "<<wetranges->size();
         throw std::runtime_error(exceptionMessage.str());
      }
      // cache the netcdf file handle to try to trick the netcdf library to allow access via independent threads (this seems to require a thread safe hdf5 library linked to netcdf if Ncfile::nc4)
      netcdf_file = std::unique_ptr<netCDF::NcFile>(new NcFile(filepath, NcFile::replace, NcFile::classic)); // as we currently use no compression: use classic format to not rely on having a thread-safe hdf5
      OceanLevelsOutputWriter::create_output_file(*netcdf_file, var_name, dim_name, nod2Dcount, levelcount, var_description, var_unit, global_attribute_name, global_attribute_txt, units_txt, calendar_txt);

      //netcdf_file->sync();
      int retCode = nc_sync(netcdf_file->getId());
      if(retCode!=NC_NOERR)
      {
         std::stringstream exceptionMessage;
         exceptionMessage << __FILE__ << ":" << __LINE__ << " netcdf error";
         throw std::runtime_error(exceptionMessage.str());
      }
   }
   
   
   // fesom data stored in a 2D netcdf array, suitable to store fesom 3D output in separate ocean levels
   void OceanLevelsOutputWriter::appendOutput(const int size, const float *data, const double ysec)
   {
      NcVar time_var = netcdf_file->getVar("time");
      const size_t save_count = time_var.getDim(0).getSize();
      NcVar var = netcdf_file->getVar(var_name);
      vector<NcDim> dims = var.getDims();
      
      vector<size_t> timepos = {save_count};
      time_var.putVar(timepos, ysec);
      
      vector<size_t>::size_type dataindex1d = 0;
      for(int l = 0; l < wetranges->size(); l++) {
         for(int r = 0; r < (*wetranges)[l].size(); r++) {
            OceanLevelsOutputWriter::writeFragmentToNetcdf(var, save_count, data
                                                           , l // levelindex
                                                           , (*wetranges)[l][r].first // columnindex
                                                           , (*wetranges)[l][r].second // columncount
                                                           , dataindex1d // dataindex
            );
            dataindex1d += (*wetranges)[l][r].second;
         }
      }
      
      //netcdf_file->sync();
      int retCode = nc_sync(netcdf_file->getId());
      if(retCode!=NC_NOERR)
      {
         std::stringstream exceptionMessage;
         exceptionMessage << __FILE__ << ":" << __LINE__ << " netcdf error";
         throw std::runtime_error(exceptionMessage.str());
      }
   }

   
   void OceanLevelsOutputWriter::writeFragmentToNetcdf(netCDF::NcVar &var, const size_t save_count, const float *var_data, const int levelindex, const int columnindex, const int columncount, const int dataindex)
   {
      assert(var.getDimCount() == 3);
      assert(save_count < var.getDim(0).getSize());
      assert(levelindex < var.getDim(1).getSize());
      assert(columnindex < var.getDim(2).getSize());

      vector<size_t> start(3); // 3 dims, i.e. var.getDimCount()
      vector<size_t> count(3);
      // time
      start[0] = save_count;
      count[0] = 1;
      // down
      start[1] = levelindex;
      count[1] = 1;
      // right
      start[2] = columnindex;
      count[2] = columncount; // number of variables to be written from the data array, i.e. this value + dataindex may not exceed the array bounds
      var.putVar(start, count, &var_data[dataindex]);
      
      // missing values are filled automatically with netcdf fill values
   }
   
   
   void OceanLevelsOutputWriter::create_output_file(netCDF::NcFile &ncf, std::string varname, std::string dimname, int dimsize_x, int dimsize_y, std::string description_txt, std::string unit_txt, std::string global_attribute_name, std::string global_attribute_txt, std::string units_txt, std::string calendar_txt)
   {
      NcDim dim_rec = OutputWriterUtils::create_shared_header(ncf, units_txt, calendar_txt);
      
      string dimname_depth = "depth";
      
      NcDim dim_nod2d = ncf.addDim(dimname, dimsize_x);
      NcDim dim_depth = ncf.addDim(dimname_depth, dimsize_y);
      NcVar depth_var = ncf.addVar(dimname_depth, ncDouble, dim_depth);
      depth_var.putAtt("units", "m");
      depth_var.putAtt("long_name", "depth");
      depth_var.putAtt("positive", "down");
      depth_var.putAtt("axis", "Z");

      vector<NcDim> dim_xd;
      dim_xd.push_back(dim_rec);
      dim_xd.push_back(dim_depth);
      dim_xd.push_back(dim_nod2d);
      NcVar var = ncf.addVar(varname, ncFloat, dim_xd);
      
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

      depth_var.putVar(&(*depths)[0]);
   }
}
