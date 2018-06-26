#include "OutputWriterUtils.h"
#include <netcdf>

using namespace std;
using namespace netCDF;

namespace AWI
{
   
   netCDF::NcDim OutputWriterUtils::create_shared_header(netCDF::NcFile &ncf, std::string units_txt, std::string calendar_txt)
   {
      NcDim dim_rec = ncf.addDim("time");
      
		// NetCDF Climate and Forecast (CF) Metadata
      NcVar time_var = ncf.addVar("time", ncDouble, dim_rec);
		time_var.putAtt("long_name", "time");
		time_var.putAtt("units", units_txt);
		time_var.putAtt("calendar", calendar_txt);
      
      return dim_rec;
   }

}
