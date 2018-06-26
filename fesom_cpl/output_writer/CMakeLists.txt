cmake_minimum_required(VERSION 3.4)

project(outputwriter CXX C Fortran) #  the FortranCInterface requires the C language to be enabled

# get our source files
set(src_home ${CMAKE_CURRENT_LIST_DIR}) # path to src directory starting from the dir containing our CMakeLists.txt
file(GLOB sources_CXX ${src_home}/*.cpp)

include(FortranCInterface)
FortranCInterface_HEADER(OutputWriterFCMacros.h MACRO_NAMESPACE "OutputWriterFCMacros_" SYMBOLS set_output_file begin_append_output end_append_output set_ocean_levels_output_file add_ocean_levels_wetrange add_ocean_depth)

include(${CMAKE_CURRENT_LIST_DIR}/../../cmake/FindNETCDF.cmake)

add_library(${PROJECT_NAME} ${sources_CXX})
target_link_libraries(${PROJECT_NAME} ${NETCDF_CXX_LIBRARIES})
target_include_directories(${PROJECT_NAME}
	PRIVATE ${NETCDF_CXX_INCLUDE_DIRECTORIES}
	INTERFACE ${CMAKE_CURRENT_LIST_DIR}
	PUBLIC ${CMAKE_CURRENT_BINARY_DIR}
)
target_compile_options(${PROJECT_NAME} PRIVATE -std=c++11)
#target_compile_options(${PROJECT_NAME} PRIVATE -gcc-name=/opt/gcc/4.8.1/bin/gcc -gxx-name=/opt/gcc/4.8.1/bin/g++ -std=gnu++11)
#target_compile_options(${PROJECT_NAME} PRIVATE -gcc-name=/opt/gcc/4.8.1/bin/gcc -gxx-name=/opt/gcc/4.8.1/bin/g++ -std=c++11)
