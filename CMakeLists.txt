cmake_minimum_required(VERSION 3.16)

#----- Append EC specific module path
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
foreach(PATH $ENV{EC_CMAKE_MODULE_PATH})
   list(APPEND CMAKE_MODULE_PATH ${PATH})
endforeach()

include(ec_init)           # Initialise compilers and ec specific functions
ec_parse_manifest()        # Parse MANIFEST file
ec_build_info()            # Generate build include file
ec_git_version()
#include(doxygen)    # Doxygen target

project(${NAME} DESCRIPTION "${DESCRIPTION}")
set(PROJECT_VERSION ${GIT_VERSION})
option(BUILD_SHARED_LIBS "Build shared libraries instead of static ones." TRUE)

set(CMAKE_INSTALL_PREFIX "" CACHE PATH "..." FORCE)

enable_language(C)
enable_language(Fortran)
enable_testing()

include(ec_compiler_presets)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp")

link_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib)
add_compile_definitions(HAVE_GPC)

#TODO --- TEST with only static
if(BUILD_SHARED_LIBS)
   find_package(RMN ${RMN_REQ_VERSION} COMPONENTS SHARED THREADED REQUIRED)
else()
   find_package(RMN ${RMN_REQ_VERSION} COMPONENTS REQUIRED)
endif()
if (RMN_FOUND) 
   add_compile_definitions(HAVE_RMN)
   include_directories(${RMN_INCLUDE_DIR})
endif()

#find_package(VGRID ${VGRID_REQ_VERSION} COMPONENTS shared OPTIONAL)
#if(VGRID_FOUND)
#  add_compile_definitions(HAVE_VGRID)
#   include_directories(${VGRID_INCLUDE_DIR})
#endif()

find_package(GDAL)
if(GDAL_FOUND)
   add_compile_definitions(HAVE_GDAL)
   include_directories(${GDAL_INCLUDE_DIR})
endif()

if(${CMAKE_SYSTEM_PROCESSOR} EQUAL x86_64)
   add_compile_options(-fPIC -DSTDC_HEADERS)
endif()

if(DEFINED ENV{OMPI})
   message(STATUS "OMPI build enabled")

   find_package(MPI REQUIRED)
   include_directories(${MPI_INCLUDE_PATH})
   add_compile_definitions(_MPI)

   find_package(OpenMP REQUIRED) 
   add_compile_options(${OpenMP_C_FLAGS})
endif()

add_compile_definitions(_${CMAKE_SYSTEM_NAME}_ _GNU_SOURCE)

add_subdirectory(lib lib)
if(NOT DEFINED ENV{OMPI})
   add_subdirectory(util util)
   add_subdirectory(test test)
endif()

configure_file(config.in ${CMAKE_BINARY_DIR}/${NAME}-config @ONLY)
install(PROGRAMS ${CMAKE_BINARY_DIR}/${NAME}-config DESTINATION bin)