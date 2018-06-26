#ifndef FC_HEADER_INCLUDED
#define FC_HEADER_INCLUDED

/* Mangling for Fortran global symbols without underscores. */
#define FC_GLOBAL(name,NAME) name##_

/* Mangling for Fortran global symbols with underscores. */
#define FC_GLOBAL_(name,NAME) name##_

/* Mangling for Fortran module symbols without underscores. */
#define FC_MODULE(mod_name,name, mod_NAME,NAME) mod_name##_mp_##name##_

/* Mangling for Fortran module symbols with underscores. */
#define FC_MODULE_(mod_name,name, mod_NAME,NAME) mod_name##_mp_##name##_

/*--------------------------------------------------------------------------*/
/* Mangle some symbols automatically.                                       */
#define stopwatch_start FC_GLOBAL_(stopwatch_start, STOPWATCH_START)
#define stopwatch_start_id FC_GLOBAL_(stopwatch_start_id, STOPWATCH_START_ID)
#define stopwatch_stop FC_GLOBAL_(stopwatch_stop, STOPWATCH_STOP)
#define stopwatch_print_all FC_GLOBAL_(stopwatch_print_all, STOPWATCH_PRINT_ALL)

#endif
