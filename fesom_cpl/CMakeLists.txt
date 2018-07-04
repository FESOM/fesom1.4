cmake_minimum_required(VERSION 3.4)

project(fesom C Fortran)

set(FESOM_SHORT_WAVE_PENETRATION ON CACHE BOOL "compile fesom with short wave penetration enabled")

# get our source files
set(src_home ${CMAKE_CURRENT_LIST_DIR}) # path to src directory starting from the dir containing our CMakeLists.txt
file(GLOB sources_Fortran ${src_home}/*.F90)
if(${FESOM_RECOM})
   	target_compile_definitions(${PROJECT_NAME} PRIVATE __bgc)
	if(${RECOM_LIBRARY})
		include(${CMAKE_CURRENT_LIST_DIR}/../cmake/FindREcoM.cmake)
	else()
		file(GLOB sources_Fortran ${src_home}/*.F90 ${src_home}/../recom/*.F90)
	endif()
endif()

list(REMOVE_ITEM sources_Fortran ${src_home}/fesom_partition_init.F90)
list(REMOVE_ITEM sources_Fortran ${src_home}/ice_thermo.F90)
list(REMOVE_ITEM sources_Fortran ${src_home}/cpl_ncout.F90)
list(REMOVE_ITEM sources_Fortran ${src_home}/cpl_str2unstr.F90)
file(GLOB sources_C ${src_home}/*.c)

if(NOT ${FESOM_COUPLED})
   list(REMOVE_ITEM sources_Fortran ${src_home}/cpl_driver.F90)
endif()

add_subdirectory(output_scheduler)
add_subdirectory(output_writer)
add_subdirectory(../stopwatch ${PROJECT_BINARY_DIR}/stopwatch)

# depends on the metis library
add_subdirectory(../lib/metis-4.0/Lib ${PROJECT_BINARY_DIR}/metis)
# depends on the parms library
add_subdirectory(../lib/parms ${PROJECT_BINARY_DIR}/parms)

include(${CMAKE_CURRENT_LIST_DIR}/../cmake/FindNETCDF.cmake)

add_library(${PROJECT_NAME}_C ${sources_C})
target_compile_definitions(${PROJECT_NAME}_C PRIVATE PARMS USE_MPI REAL=double DBL HAS_BLAS FORTRAN_UNDERSCORE VOID_POINTER_SIZE_8 SGI LINUX UNDER_ MPI2)
target_link_libraries(${PROJECT_NAME}_C metis parms)

# create our binary (set its name to name of this project)
add_executable(${PROJECT_NAME} ${sources_Fortran})
target_compile_definitions(${PROJECT_NAME} PRIVATE PARMS use_tracer_fct calcmeans use_diag use_ice use_ice_fct cpl_str allow_calcmeans allow_diag)
if(${FESOM_COUPLED})
   include(${CMAKE_CURRENT_LIST_DIR}/../cmake/FindOASIS.cmake)
   target_compile_definitions(${PROJECT_NAME} PRIVATE __oasis)
endif()
if(${FESOM_SHORT_WAVE_PENETRATION})
   target_compile_definitions(${PROJECT_NAME} PRIVATE use_sw_pene)
endif()
# CMAKE_Fortran_COMPILER_ID will also work if a wrapper is being used (e.g. mpif90 wraps ifort -> compiler id is Intel)
if(${CMAKE_Fortran_COMPILER_ID} STREQUAL  Intel )
	target_compile_options(${PROJECT_NAME} PRIVATE -r8 -i4 -fp-model precise -no-prec-div -no-prec-sqrt -fast-transcendentals -xHost -ip)
elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL  GNU )
	target_compile_options(${PROJECT_NAME} PRIVATE -fdefault-real-8 -ffree-line-length-none)
endif()
target_include_directories(${PROJECT_NAME} PRIVATE ${NETCDF_Fortran_INCLUDE_DIRECTORIES} ${OASIS_Fortran_INCLUDE_DIRECTORIES})
target_include_directories(${PROJECT_NAME} PRIVATE ${RECOM_Fortran_INCLUDE_DIRECTORIES})
target_include_directories(${PROJECT_NAME} PRIVATE ${MCT_Fortran_INCLUDE_DIRECTORIES} ${MPEU_Fortran_INCLUDE_DIRECTORIES})
target_include_directories(${PROJECT_NAME} PRIVATE ${SCRIP_Fortran_INCLUDE_DIRECTORIES})
target_link_libraries(${PROJECT_NAME} ${PROJECT_NAME}_C outputscheduler outputwriter stopwatch ${NETCDF_Fortran_LIBRARIES} ${NETCDF_C_LIBRARIES} ${OASIS_Fortran_LIBRARIES} ${RECOM_Fortran_LIBRARIES})
target_link_libraries(${PROJECT_NAME} ${PROJECT_NAME}_C ${MCT_Fortran_LIBRARIES} ${MPEU_Fortran_LIBRARIES} ${SCRIP_Fortran_LIBRARIES})
set_target_properties(${PROJECT_NAME} PROPERTIES LINKER_LANGUAGE Fortran)

set(FESOM_INSTALL_FILEPATH "${CMAKE_CURRENT_LIST_DIR}/../bin/fesom" CACHE FILEPATH "file path where the FESOM binary should be put")
get_filename_component(FESOM_INSTALL_PATH ${FESOM_INSTALL_FILEPATH} DIRECTORY)
get_filename_component(FESOM_INSTALL_NAME ${FESOM_INSTALL_FILEPATH} NAME)
install(PROGRAMS ${PROJECT_BINARY_DIR}/${PROJECT_NAME} DESTINATION ${FESOM_INSTALL_PATH} RENAME ${FESOM_INSTALL_NAME})
