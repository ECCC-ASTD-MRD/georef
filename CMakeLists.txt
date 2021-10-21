cmake_minimum_required(VERSION 3.16)

#----- Append EC specific module path
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake_rpn)
foreach(PATH $ENV{EC_CMAKE_MODULE_PATH})
   list(APPEND CMAKE_MODULE_PATH ${PATH})
endforeach()

include(ec_init)           # Initialise compilers and ec specific functions
ec_git_version()           # Get version from git state
ec_parse_manifest()        # Parse MANIFEST file

project(${NAME} DESCRIPTION "${DESCRIPTION}")
set(PROJECT_VERSION ${VERSION}${STATE})
set(LIBGEOREF_VERSION ${PROJECT_VERSION} CACHE INTERNAL "libgeoref version" FORCE) # Needed for cascaded version identification

option(BUILD_SHARED_LIBS "Build shared libraries instead of static ones." TRUE)

ec_build_info()            # Generate build include file
ec_install_prefix()        # Define install prefix  

enable_language(C)
enable_language(Fortran)
enable_testing()

#include(doxygen)    # Doxygen target
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

#TODO: might have to do this - type gdal-config >/dev/null 2>&1 && CMAKE_COMP_FLAGS="$CMAKE_COMP_FLAGS -DGDAL_ROOT=`gdal-config --prefix`" 2>/dev/null
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
#   add_subdirectory(test test)
endif()

configure_file(config.in ${CMAKE_BINARY_DIR}/${NAME}-config @ONLY)
install(PROGRAMS ${CMAKE_BINARY_DIR}/${NAME}-config DESTINATION bin)

ec_build_config()                                             # Create build configuration script

# Packaging with CPack
ec_prepare_ssm()          # Prepare ssm control file and post/pre install script

set(CPACK_GENERATOR "TGZ")
set(CPACK_PACKAGE_VENDOR "ECCC")
set(CPACK_PACKAGE_CONTACT "${MAINTAINER}")
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE")
set(CPACK_RESOURCE_FILE_README "${CMAKE_CURRENT_SOURCE_DIR}/README.org")
set(CPACK_OUTPUT_FILE_PREFIX "${CMAKE_CURRENT_BINARY_DIR}/package")
set(CPACK_PACKAGE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set(CPACK_PACKAGE_FILE_NAME "${CMAKE_INSTALL_PREFIX}.ssm")
set(CPACK_SOURCE_PACKAGE_FILE_NAME "${CMAKE_INSTALL_PREFIX}.src")
include(CPack)
