#include "OutputSchedulerManager.h"
#include <iomanip>
#include <sstream>
#include <cstdlib>
#include "yaml-cpp/yaml.h"

using namespace std;

namespace AWI
{
   OutputSchedulerManager::OutputSchedulerManager()
   : schedules()
   {
      string env_key = "AWI_FESOM_YAML";
      string yaml_txt;
      char const* tmp = getenv(env_key.c_str());
      if(tmp == NULL)
         yaml_txt = OutputSchedulerManager::defaultSchedules();
      else
         yaml_txt = std::string(tmp);
      
      istringstream yaml_stream(yaml_txt);
      YAML::Parser parser(yaml_stream);
      YAML::Node yaml;
      if(!parser.GetNextDocument(yaml))
         return;
      
      const YAML::Node& schedules_yaml = yaml["output_schedules"];
      
      
      for(unsigned i=0;i<schedules_yaml.size();i++)
      {
         const YAML::Node& entry = schedules_yaml[i];
         OutputScheduler schedule = OutputSchedulerManager::parseSchedule(entry);
         const YAML::Node& entry_vars = entry["vars"];
         for(unsigned ii=0;ii<entry_vars.size();ii++)
         {
            string key;
            entry_vars[ii] >> key;
            std::pair<std::map<std::string, OutputScheduler>::iterator,bool> return_pair;
            return_pair = schedules.insert(map<string, OutputScheduler>::value_type(key, schedule));
            if(return_pair.second == false)
            {
               std::stringstream exceptionMessage;
               exceptionMessage << __FILE__ << ":" << __LINE__ <<" duplicate output schedule: "<<key;
               throw std::runtime_error(exceptionMessage.str());
            }
         }
      }
   }
   
   
   bool OutputSchedulerManager::shouldOutput(std::string name, int istep, double dsec, int day_step, int month_step, int year, bool is_last_day_in_month)
   {
      bool should_output = false;
   
      bool is_last_second_in_day = (int)dsec==3600*24;
      bool is_last_month_in_year = month_step % 12 == 0;
      
      map<string, OutputScheduler>::iterator it = schedules.find(name);
      if(it != schedules.end())
      {
         OutputScheduler &s = it->second;
         
         if(s.unit == "y")
         {
            if(is_last_month_in_year && is_last_day_in_month && is_last_second_in_day)
               if(s.first <= year && year% s.rate == 0)
                  should_output = true;
         }
         else if(s.unit == "m")
         {
            if(is_last_day_in_month && is_last_second_in_day)
               if(s.first <= month_step && month_step% s.rate == 0)
                  should_output = true;
         }
         else if(s.unit == "d")
         {
            if(is_last_second_in_day)
               if(s.first <= day_step && day_step% s.rate == 0)
                  should_output = true;
         }
         else if(s.unit == "s")
         {
            if(s.first <= istep && istep% s.rate == 0)
               should_output = true;
         }
         else
         {
            std::stringstream exceptionMessage;
            exceptionMessage << __FILE__ << ":" << __LINE__ <<" unknown unit: "<<s.unit;
            throw std::runtime_error(exceptionMessage.str());
         }
      }
      else
      {
         std::stringstream exceptionMessage;
         exceptionMessage << __FILE__ << ":" << __LINE__ <<" unknown output schedule: "<<name;
         throw std::runtime_error(exceptionMessage.str());
      }
      
      return should_output;
   }
   
   
   void OutputSchedulerManager::scheduleInfo(std::string name, std::string& unit, unsigned int& first, unsigned int& rate)
   {
      map<string, OutputScheduler>::iterator it = schedules.find(name);
      if(it != schedules.end())
      {
         OutputScheduler &s = it->second;
         unit = s.unit;
         first = s.first;
         rate = s.rate;
      }
   }


   void OutputSchedulerManager::print()
   {
      std::stringstream msg;
      for (map<string, OutputScheduler>::iterator it = schedules.begin(); it != schedules.end(); ++it)
      {
         const string &n = it->first;
         const OutputScheduler &s = it->second;
         
         msg<<"name: \""<<n<<"\" unit: "<<s.unit<<" first: "<<s.first<<" rate: "<<s.rate<<"\n";
      }
      cout<<msg.str();
   }

   
   void OutputSchedulerManager::printDefaultSchedules()
   {
      cout<<OutputSchedulerManager::defaultSchedules()<<"\n";
   }
   
   
   std::string OutputSchedulerManager::defaultSchedules()
   {
      return "---\n"
      "output_schedules:\n"
      "  - vars: [uo, vo, wo,sisnthick,u2o,v2o,rho,urho,vrho,uv,mlotst,uto,vto,uso,vso,thdgr,thdgrsn,uhice,vhice,uhsnow,vhsnow,flice,tair,shum,uwind,vwind,prlq,prsn,runoff,evap,lwrd,rsdo,hfds,wnet,olat,osen,olwout,virtual_salt,relax_salt,tauuo,tauvo,sitimefrac,w2o,wso,wto,omldamax,zossq,sisnmass,sidmassth,sidmasssi,sidmasstranx,sidmasstrany,volo,opottemptend,pbo,soga,thetaoga,tos,tso,sos,siarean,siareas,siextentn,siextents,sivoln,sivols,evs,sidmassevapsubl,sifllatstop,sistrxdtop,sistrydtop,sistrxubot,sistryubot,wfo,fsitherm,sisnconc]\n"
      "    unit: m\n"
      "    rate: 1\n"
      "  - vars: [zos,thetao,so,sic,sithick,siu,siv,sispeed,sivol]\n"
      "    unit: d\n"
      "    rate: 1\n"
      "  - vars: [restart]\n"
      "    unit: m\n"
      "    first: 72\n"
      "    rate: 1";
   }

   
   OutputScheduler OutputSchedulerManager::parseSchedule(const YAML::Node& schedule)
   {
      // fetch the mandatory items
      string unit;
      schedule["unit"] >> unit;
      int rate;
      schedule["rate"] >> rate;
      
      int first = 1;
      // fetch the optional items
      if(const YAML::Node *item = schedule.FindValue("first")) {
         *item >> first;
      }
      
      return OutputScheduler(unit, first, rate);
   }

}
