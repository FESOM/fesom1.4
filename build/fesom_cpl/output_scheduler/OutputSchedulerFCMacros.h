#ifndef OutputSchedulerFCMacros_HEADER_INCLUDED
#define OutputSchedulerFCMacros_HEADER_INCLUDED

/* Mangling for Fortran global symbols without underscores. */
#define OutputSchedulerFCMacros_GLOBAL(name,NAME) name##_

/* Mangling for Fortran global symbols with underscores. */
#define OutputSchedulerFCMacros_GLOBAL_(name,NAME) name##_

/* Mangling for Fortran module symbols without underscores. */
#define OutputSchedulerFCMacros_MODULE(mod_name,name, mod_NAME,NAME) mod_name##_mp_##name##_

/* Mangling for Fortran module symbols with underscores. */
#define OutputSchedulerFCMacros_MODULE_(mod_name,name, mod_NAME,NAME) mod_name##_mp_##name##_

/*--------------------------------------------------------------------------*/
/* Mangle some symbols automatically.                                       */
#define should_output OutputSchedulerFCMacros_GLOBAL_(should_output, SHOULD_OUTPUT)
#define print_schedules OutputSchedulerFCMacros_GLOBAL_(print_schedules, PRINT_SCHEDULES)
#define schedule_info OutputSchedulerFCMacros_GLOBAL_(schedule_info, SCHEDULE_INFO)
#define print_default_schedules OutputSchedulerFCMacros_GLOBAL_(print_default_schedules, PRINT_DEFAULT_SCHEDULES)

#endif
