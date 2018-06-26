#ifndef OutputWriterManager_9A2312DD_D4B5_4779_9D40_02D62F35DBAE
#define OutputWriterManager_9A2312DD_D4B5_4779_9D40_02D62F35DBAE

#include <string>
#include <map>
#include <memory>
#include <vector>
#include <thread>
#include "OutputWriterThing.h"


namespace AWI
{
   
   class OutputWriterManager
   {
      public:
      void addOceanLevelsWetrange(const int level, const int index, const int size);
      void addOceanDepth(const double d);
      void setOutputFile(const std::string filepath, const std::string var_name, const int var_size, const std::string var_description, const std::string var_unit, const std::string dim_name, const std::string global_attribute_name, std::string global_attribute_txt, std::string units_txt, std::string calendar_txt);
      void setOceanLevelsOutputFile(const std::string filepath, const std::string var_name, const int nod2Dcount, const int levelcount, const std::string var_description, const std::string var_unit, const std::string dim_name, const std::string global_attribute_name, std::string global_attribute_txt, std::string units_txt, std::string calendar_txt);
      void beginAppendOutput(const std::string var_name, const int var_size, const float *var_data, const double ysec);
      void endAppendOutput(const std::string var_name);
      
   private:
      std::vector< std::vector< std::pair<size_t,size_t> > > wetranges;
      std::vector<double> depths;
      std::map< std::string, std::unique_ptr<OutputWriterThing> > writers;
      std::map<std::string, std::thread*> writer_threads;
   };
   
}

#endif
