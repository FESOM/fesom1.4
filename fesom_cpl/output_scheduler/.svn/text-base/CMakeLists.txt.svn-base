cmake_minimum_required(VERSION 3.4)

project(outputscheduler CXX C Fortran) #  the FortranCInterface requires the C language to be enabled

# get our source files
set(src_home ${CMAKE_CURRENT_LIST_DIR}) # path to src directory starting from the dir containing our CMakeLists.txt
file(GLOB sources_CXX ${src_home}/*.cpp)

add_subdirectory(../../3rdparty/yaml-cpp ${PROJECT_BINARY_DIR}/yaml-cpp)

include(FortranCInterface)
FortranCInterface_HEADER(OutputSchedulerFCMacros.h MACRO_NAMESPACE "OutputSchedulerFCMacros_" SYMBOLS should_output print_schedules schedule_info print_default_schedules)

add_library(${PROJECT_NAME} ${sources_CXX})
target_link_libraries(${PROJECT_NAME} yaml-cpp)
target_include_directories(${PROJECT_NAME}
	PUBLIC ${CMAKE_CURRENT_BINARY_DIR}
	INTERFACE ${CMAKE_CURRENT_LIST_DIR}
)
