cmake_minimum_required(VERSION 3.16)

# Append EC specific module path
list(APPEND CMAKE_MODULE_PATH $ENV{EC_CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake_rpn/modules)

include(ec_init)           # Initialise compilers and ec specific functions
ec_git_version()           # Get version from git state
ec_parse_manifest()        # Parse MANIFEST file

project(${NAME} DESCRIPTION "${DESCRIPTION}")
set(PROJECT_VERSION ${VERSION}${STATE})
set(${NAME}_VERSION ${PROJECT_VERSION} CACHE INTERNAL "${NAME} version" FORCE) # Needed for cascaded version identification

option(BUILD_SHARED_LIBS "Build shared libraries instead of static ones." TRUE)

ec_build_info()            # Generate build include file

enable_language(C)
enable_language(Fortran)
enable_testing()

#include(doxygen)    # Doxygen target
include(ec_compiler_presets)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp")

link_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib)

if(${CMAKE_SYSTEM_PROCESSOR} EQUAL x86_64)
   add_compile_options(-fPIC -DSTDC_HEADERS)
endif()
add_compile_definitions(HAVE_GPC)

find_package(MPI REQUIRED)
#set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
#set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
add_compile_definitions(HAVE_MPI)

find_package(rmn ${rmn_REQ_VERSION} REQUIRED)
#find_package(vgrid ${vgrid_REQ_VERSION} REQUIRED)

execute_process(
   # Use this command because the output of gdal-config may be wrong depending
   # on how it was installed, but it's location can still be relied on.
   COMMAND bash -c "cd $(dirname $(which gdal-config))/.. && pwd"
   WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
   RESULT_VARIABLE GDAL_RESULT
   OUTPUT_VARIABLE GDAL_ROOT
   ERROR_VARIABLE GDAL_ERROR
   OUTPUT_STRIP_TRAILING_WHITESPACE
   ERROR_STRIP_TRAILING_WHITESPACE
)

find_package(GDAL ${GDAL_REQ_VERSION})
if(GDAL_FOUND)
   target_compile_definitions(GDAL::GDAL INTERFACE HAVE_GDAL)
endif()

find_package(OpenMP)
if (OpenMP_FOUND) 
   add_compile_options(${OpenMP_C_FLAGS})
   add_compile_definitions(HAVE_OPENMP)
endif()

add_compile_definitions(_${CMAKE_SYSTEM_NAME}_ _GNU_SOURCE)

link_libraries("-Wl,--disable-new-dtags")
add_subdirectory(lib lib)
add_subdirectory(test/C test)

configure_file(config.in ${CMAKE_BINARY_DIR}/${NAME}-config @ONLY)
install(PROGRAMS ${CMAKE_BINARY_DIR}/${NAME}-config DESTINATION bin)

#----- Generate the config file for the project to be usable via cmake's find_package command
set(INCLUDE_INSTALL_DIR include)
set(LIB_INSTALL_DIR     lib)
set(CONFIG_INSTALL_DIR  "${LIB_INSTALL_DIR}/cmake/${PROJECT_NAME}-${PROJECT_VERSION}")

install(EXPORT ${PROJECT_NAME}-shared-targets
   NAMESPACE ${PROJECT_NAME}::
   DESTINATION ${CONFIG_INSTALL_DIR})

install(EXPORT ${PROJECT_NAME}-static-targets
   NAMESPACE ${PROJECT_NAME}::
   DESTINATION ${CONFIG_INSTALL_DIR})

include(CMakePackageConfigHelpers)
configure_package_config_file(
    "Config.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
    INSTALL_DESTINATION "${CONFIG_INSTALL_DIR}"
    PATH_VARS           INCLUDE_INSTALL_DIR LIB_INSTALL_DIR
)
write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
    COMPATIBILITY SameMajorVersion
)
install(FILES   "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
                "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
    DESTINATION "${CONFIG_INSTALL_DIR}"
)

# Packaging with CPack
ec_install_prefix(${NAME} ${PROJECT_VERSION} PACKAGE_PREFIX)  # Define package prefix  
ec_build_config()                                             # Create build configuration script
ec_prepare_ssm()                                              # Prepare ssm packaging files

set(CPACK_GENERATOR "TGZ")
set(CPACK_PACKAGE_VENDOR "ECCC")
set(CPACK_PACKAGE_CONTACT "${MAINTAINER}")
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE")
set(CPACK_RESOURCE_FILE_README "${CMAKE_CURRENT_SOURCE_DIR}/README.md")
set(CPACK_OUTPUT_FILE_PREFIX "${CMAKE_CURRENT_BINARY_DIR}/package")
set(CPACK_PACKAGE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set(CPACK_PACKAGE_FILE_NAME "${PACKAGE_PREFIX}")
set(CPACK_SOURCE_PACKAGE_FILE_NAME "${NAME}_${PROJECT_VERSION}")
include(CPack)
