cmake_minimum_required(VERSION 3.16)
project(libgeoref VERSION 0.1 DESCRIPTION "libgeoref")

enable_language(Fortran)

list(APPEND CMAKE_MODULE_PATH $ENV{EC_CMAKE_MODULE_PATH})

message(STATUS "Generating libgeoref package")
message(STATUS "Build architecture ${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}")

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
link_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib)

find_package(RMN 19.6 REQUIRED)
if (RMN_FOUND) 
   add_compile_definitions(HAVE_RMN)
   include_directories(${RMN_INCLUDE_DIR})
endif()

if(${CMAKE_SYSTEM_PROCESSOR} EQUAL x86_64)
   add_compile_options(-fPIC -DSTDC_HEADERS)
endif()

find_package(OpenMP REQUIRED) 
add_compile_options(${OpenMP_C_FLAGS})

set(CMAKE_BUILD_TYPE "Release")
include("$ENV{EC_CMAKE_MODULE_PATH}/CMake_compiler.txt")

add_compile_definitions(VERSION=\"$ENV{VERSION}-r$ENV{BUILDINFO}\" _${CMAKE_SYSTEM_NAME}_ TCL_THREADS _GNU_SOURCE)
add_subdirectory(src src)



