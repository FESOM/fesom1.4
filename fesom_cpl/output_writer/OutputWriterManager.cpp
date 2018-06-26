#include "OutputWriterManager.h"
#include <stdexcept>
#include "OutputWriter.h"
#include "OceanLevelsOutputWriter.h"

using namespace std;

namespace AWI
{
   
   void OutputWriterManager::addOceanLevelsWetrange(const int level, const int index, const int size)
   {
      if(wetranges.size() < level+1)
         wetranges.resize(level+1);
      
      wetranges[level].emplace_back(index, size);
   }

   
   void OutputWriterManager::addOceanDepth(const double d)
   {
      depths.push_back(d);
   }
   
   
   void OutputWriterManager::setOutputFile(const std::string filepath, const std::string var_name, const int var_size, const std::string var_description, const std::string var_unit, const std::string dim_name, const std::string global_attribute_name, std::string global_attribute_txt, std::string units_txt, std::string calendar_txt)
   {
      map<string, unique_ptr<OutputWriterThing> >::iterator it = writers.find(var_name);
      if(it != writers.end())
         writers.erase(it);
      
      writers[var_name] = unique_ptr<OutputWriter>(new OutputWriter(filepath, var_name, var_size, var_description, var_unit, dim_name, global_attribute_name, global_attribute_txt, units_txt, calendar_txt));
   }

   
   void OutputWriterManager::setOceanLevelsOutputFile(const std::string filepath, const std::string var_name, const int nod2Dcount, const int levelcount, const std::string var_description, const std::string var_unit, const std::string dim_name, const std::string global_attribute_name, std::string global_attribute_txt, std::string units_txt, std::string calendar_txt)
   {
      map<string, unique_ptr<OutputWriterThing> >::iterator it = writers.find(var_name);
      if(it != writers.end())
         writers.erase(it);

      writers[var_name] = unique_ptr<OceanLevelsOutputWriter>(new OceanLevelsOutputWriter(filepath, var_name, nod2Dcount, levelcount, &depths, &wetranges, var_description, var_unit, dim_name, global_attribute_name, global_attribute_txt, units_txt, calendar_txt));
   }

   
   void OutputWriterManager::beginAppendOutput(const std::string var_name, const int var_size, const float *var_data, const double ysec)
   {
      map<string, unique_ptr<OutputWriterThing> >::iterator it = writers.find(var_name);
      if(it != writers.end())
      {
         unique_ptr<OutputWriterThing> &w = it->second;         
         writer_threads[var_name] = new thread(&OutputWriterThing::appendOutput, w.get(), var_size, var_data, ysec);
      }
      else
      {
         std::stringstream exceptionMessage;
         exceptionMessage << __FILE__ << ":" << __LINE__ <<" unknown output writer: "<<var_name;
         throw std::runtime_error(exceptionMessage.str());
      }

   }
   
   
   void OutputWriterManager::endAppendOutput(const std::string var_name)
   {
      map<string, thread*>::iterator it = writer_threads.find(var_name);
      if(it != writer_threads.end())
      {
         thread *t = it->second;
         t->join();
         writer_threads.erase(it);
         delete t;
      }
   }
   
}
