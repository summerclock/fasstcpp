cmake_minimum_required(VERSION 3.16.0 FATAL_ERROR)
project(fasst VERSION 0.1.0)

# execute_process(
#   COMMAND git rev-parse --short HEAD
#   WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
#   OUTPUT_VARIABLE GIT_HASH
#   OUTPUT_STRIP_TRAILING_WHITESPACE
# )
# message("GIT HASH: ${GIT_HASH}")
# configure_file(
#   "${PROJECT_SOURCE_DIR}/version.hpp.in"
#   "${PROJECT_BINARY_DIR}/version.hpp"
# )
enable_testing()

include_directories("${CMAKE_SOURCE_DIR}" "${PROJECT_BINARY_DIR}" "/opt/include")
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Debug CACHE STRING "Build type (Debug, Release, etc.)" FORCE)
endif()

if(UNIX)
  set(CMAKE_SKIP_BUILD_RPATH FALSE)
  set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
  set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
  set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
endif()

if(MSVC)
  # export all symbols
  # set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS TRUE)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/bin")
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/bin")
  add_compile_options("/Zc:preprocessor")
  add_compile_definitions("GLOG_NO_ABBREVIATED_SEVERITIES")
  add_compile_definitions("NOMINMAX")
  add_compile_options(/Zi)
endif()

# Find required packages
find_package(Boost REQUIRED COMPONENTS system filesystem program_options)
find_package(glog REQUIRED)

# Add subdirectories
add_subdirectory(fasst)
target_compile_definitions(fasstlib
    PRIVATE BUILDING_DLL
)
add_subdirectory(apps)