#include "StopwatchManagerHelper.h"
#include "StopwatchManager.h"
#include <string>
#include <sstream>

using namespace std;

AWI::StopwatchManager manager;

// when calling this from fortran, make sure to terminate the string with a null char:
// call stopwatch_start_id(id, 'watch_name'//CHAR(0))
void stopwatch_start_id(int *identifier, const char* n)
{
   stringstream str;
   str << "[" << (*identifier) << "] " << n;
   manager.start(str.str());
}

// when calling this from fortran, make sure to terminate the string with a null char:
// call stopwatch_start('watch_name'//CHAR(0))
void stopwatch_start(const char* n)
{
   string str(n);
   manager.start(str);   
}


void stopwatch_stop()
{
   manager.stop();
}


void stopwatch_print_all()
{
   manager.dump(std::cout);
}