cmake_minimum_required(VERSION 3.4)

project(metis C)

# get our source files
set(src_home ${CMAKE_CURRENT_LIST_DIR}) # path to src directory starting from the dir containing our CMakeLists.txt
file(GLOB all_sources ${src_home}/*.c)

include_directories(${src_home})

# create our library (set its name to name of this project)
add_library(${PROJECT_NAME} ${all_sources})