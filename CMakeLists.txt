cmake_minimum_required(VERSION 3.4)

# set default build type cache entry (do so before project(...) is called, which would create this cache entry on its own)
if(NOT CMAKE_BUILD_TYPE)
	message(STATUS "setting default build type: Release")
	set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release RelWithDebInfo MinSizeRel.")
endif()

project(fesom)
set(TOPLEVEL_DIR ${CMAKE_CURRENT_LIST_DIR})
set(FESOM_COUPLED OFF CACHE BOOL "compile fesom standalone or with oasis support (i.e. coupled)")
set(FESOM_RECOM OFF CACHE BOOL "compile fesom with or without BGC model recom")
set(RECOM_LIBRARY OFF CACHE BOOL "use BGC model recom as library")

add_subdirectory(fesom_cpl)
