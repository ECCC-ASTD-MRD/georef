cmake_minimum_required(VERSION 3.16)
include($ENV{EC_CMAKE_MODULE_PATH}/utils.cmake)

ec_parse_version()
ec_build_include()
project(${NAME} VERSION ${VERSION} DESCRIPTION "${DESCRIPTION}")

enable_language(Fortran)

include("$ENV{EC_CMAKE_MODULE_PATH}/compiler_rules.cmake")
link_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib)

find_package(RMN 19.6 REQUIRED)
if (RMN_FOUND) 
   add_compile_definitions(HAVE_RMN)
   include_directories(${RMN_INCLUDE_DIR})
endif()
#check_version(RMN) 

find_package(OpenMP REQUIRED) 
add_compile_options(${OpenMP_C_FLAGS})

# add_link_options(-lgfortran)

add_compile_definitions(_${CMAKE_SYSTEM_NAME}_ _GNU_SOURCE)
add_subdirectory(src src)



