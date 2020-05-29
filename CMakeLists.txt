cmake_minimum_required(VERSION 3.16)

#----- Append EC specific module path
foreach(PATH $ENV{EC_CMAKE_MODULE_PATH})
   list(APPEND CMAKE_MODULE_PATH ${PATH})
endforeach()

include(ec_utils)   # Include EC specific cmake utils
ec_parse_manifest() # Parse MANIFEST file
ec_build_info()     # Generate build include file

project(${NAME} VERSION ${VERSION} DESCRIPTION "${DESCRIPTION}")

enable_language(C)
enable_language(Fortran)

include(compiler_presets)
link_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib)

set(CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS} -cpp")

find_package(RMN ${RMN_REQ_VERSION} REQUIRED)
if (RMN_FOUND) 
   add_compile_definitions(HAVE_RMN)
   include_directories(${RMN_INCLUDE_DIR})
endif()
#check_version(RMN) 

find_package(OpenMP REQUIRED) 
add_compile_options(${OpenMP_C_FLAGS})

add_compile_definitions(_${CMAKE_SYSTEM_NAME}_ _GNU_SOURCE)
add_subdirectory(src src)
add_subdirectory(test test)

include(ec_bin_config)

