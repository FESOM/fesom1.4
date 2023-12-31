cmake_minimum_required(VERSION 3.4)

project(fesom_ini C Fortran)

# get our source files
set(src_home ${CMAKE_CURRENT_LIST_DIR}/../fesom_cpl)
set(sources_Fortran ${src_home}/gen_partitioning.F90 ${src_home}/fesom_partition_init.F90 ${src_home}/oce_modules.F90 ${src_home}/gen_modules_config.F90)
set(sources_C ${src_home}/fort_part.c)

add_subdirectory(${src_home}/../lib/metis-4.0/Lib ${PROJECT_BINARY_DIR}/metis)

add_library(${PROJECT_NAME}_C ${sources_C})
target_compile_definitions(${PROJECT_NAME}_C PRIVATE PARMS USE_MPI REAL=double DBL HAS_BLAS FORTRAN_UNDERSCORE VOID_POINTER_SIZE_8 SGI LINUX UNDER_ MPI2)
target_link_libraries(${PROJECT_NAME}_C metis)

add_executable(${PROJECT_NAME} ${sources_Fortran})
target_compile_definitions(${PROJECT_NAME} PRIVATE PARMS use_tracer_fct calcmeans use_diag use_ice use_ice_fct cpl_str allow_calcmeans allow_diag use_sw_pene)
# CMAKE_Fortran_COMPILER_ID will also work if a wrapper is being used (e.g. mpif90 wraps ifort -> compiler id is Intel)
if(${CMAKE_Fortran_COMPILER_ID} STREQUAL  Intel )
	target_compile_options(${PROJECT_NAME} PRIVATE -r8 -i4 -fp-model precise -no-prec-div -no-prec-sqrt -fast-transcendentals -xHost -ip)
elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL  GNU )
	target_compile_options(${PROJECT_NAME} PRIVATE -fdefault-real-8 -ffree-line-length-none)
endif()
target_link_libraries(${PROJECT_NAME} ${PROJECT_NAME}_C)
set_target_properties(${PROJECT_NAME} PROPERTIES LINKER_LANGUAGE Fortran)
