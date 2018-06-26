cmake_minimum_required(VERSION 3.4)

project(stopwatch CXX C Fortran) #  the FortranCInterface requires the C language to be enabled

# get our source files
set(src_home ${CMAKE_CURRENT_LIST_DIR}) # path to src directory starting from the dir containing our CMakeLists.txt
file(GLOB sources_Fortran ${src_home}/*.f90)
file(GLOB sources_CXX ${src_home}/*.cpp)

# http://www.netlib.org/lapack/lawnspdf/lawn270.pdf
include(FortranCInterface)
FortranCInterface_HEADER(FC.h MACRO_NAMESPACE "FC_" SYMBOLS stopwatch_start stopwatch_start_id stopwatch_stop stopwatch_print_all)

add_library(${PROJECT_NAME} ${sources_CXX})
target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
