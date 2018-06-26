cmake_minimum_required(VERSION 3.4)

# custom build for yaml-cpp 0.3.0:
# https://github.com/jbeder/yaml-cpp/releases/tag/release-0.3.0

project(yaml-cpp CXX)

file(GLOB sources_CXX ${CMAKE_CURRENT_LIST_DIR}/yaml-cpp-release-0.3.0/src/*.cpp)

add_library(${PROJECT_NAME} ${sources_CXX})
target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_LIST_DIR}/yaml-cpp-release-0.3.0/include)
