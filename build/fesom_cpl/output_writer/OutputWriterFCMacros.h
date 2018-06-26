#ifndef OutputWriterFCMacros_HEADER_INCLUDED
#define OutputWriterFCMacros_HEADER_INCLUDED

/* Mangling for Fortran global symbols without underscores. */
#define OutputWriterFCMacros_GLOBAL(name,NAME) name##_

/* Mangling for Fortran global symbols with underscores. */
#define OutputWriterFCMacros_GLOBAL_(name,NAME) name##_

/* Mangling for Fortran module symbols without underscores. */
#define OutputWriterFCMacros_MODULE(mod_name,name, mod_NAME,NAME) mod_name##_mp_##name##_

/* Mangling for Fortran module symbols with underscores. */
#define OutputWriterFCMacros_MODULE_(mod_name,name, mod_NAME,NAME) mod_name##_mp_##name##_

/*--------------------------------------------------------------------------*/
/* Mangle some symbols automatically.                                       */
#define set_output_file OutputWriterFCMacros_GLOBAL_(set_output_file, SET_OUTPUT_FILE)
#define begin_append_output OutputWriterFCMacros_GLOBAL_(begin_append_output, BEGIN_APPEND_OUTPUT)
#define end_append_output OutputWriterFCMacros_GLOBAL_(end_append_output, END_APPEND_OUTPUT)
#define set_ocean_levels_output_file OutputWriterFCMacros_GLOBAL_(set_ocean_levels_output_file, SET_OCEAN_LEVELS_OUTPUT_FILE)
#define add_ocean_levels_wetrange OutputWriterFCMacros_GLOBAL_(add_ocean_levels_wetrange, ADD_OCEAN_LEVELS_WETRANGE)
#define add_ocean_depth OutputWriterFCMacros_GLOBAL_(add_ocean_depth, ADD_OCEAN_DEPTH)

#endif
