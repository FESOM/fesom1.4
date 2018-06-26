#ifndef OutputWriterThing_54E9D86B_2A70_409B_9A98_6D3C66DDBF34
#define OutputWriterThing_54E9D86B_2A70_409B_9A98_6D3C66DDBF34


namespace netCDF
{
   class NcFile;
   class NcDim;
}

namespace AWI
{
   
   class OutputWriterThing
   {
   public:      
      virtual void appendOutput(const int size, const float *data, const double ysec) = 0;
   };
   
}

#endif
