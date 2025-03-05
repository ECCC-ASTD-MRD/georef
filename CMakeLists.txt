cmake_minimum_required(VERSION 3.16)

# Append EC specific module path
list(APPEND CMAKE_MODULE_PATH $ENV{EC_CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake_rpn/modules)

include(ec_init)           # Initialise compilers and ec specific functions
ec_git_version()           # Get version from git state
ec_parse_manifest()        # Parse MANIFEST file

project(${NAME} DESCRIPTION "${DESCRIPTION}")
set(PROJECT_VERSION ${VERSION}${STATE})
set(${NAME}_VERSION ${PROJECT_VERSION} CACHE INTERNAL "${NAME} version" FORCE) # Needed for cascaded version identification

include(GNUInstallDirs)
set(CMAKE_INSTALL_LIBDIR lib CACHE INTERNAL "Install dir for libs" FORCE)  # Avoid using lib64, cache when in super-project
option(BUILD_SHARED_LIBS "Build shared libraries instead of static ones." TRUE)

ec_build_info()

enable_language(C)
enable_language(Fortran)
enable_testing()

include(ec_doxygen)
include(ec_compiler_presets)

set(CMAKE_CXX_STANDARD 11)

find_package(MPI REQUIRED)

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

set(MODULE_DIR ${CMAKE_CURRENT_BINARY_DIR}/include)

find_package(OpenMP)
find_package(rmn)

link_libraries("-Wl,--disable-new-dtags")
add_subdirectory(src)
add_subdirectory(python)
add_subdirectory(test/lib test)

configure_file(config.in ${CMAKE_BINARY_DIR}/${NAME}-config @ONLY)
install(PROGRAMS ${CMAKE_BINARY_DIR}/${NAME}-config DESTINATION bin)

#----- Generate the config file for the project to be usable via cmake's find_package command
set(CONFIG_INSTALL_DIR  "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}-${PROJECT_VERSION}")

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
    PATH_VARS CMAKE_INSTALL_INCLUDEDIR CMAKE_INSTALL_LIBDIR
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
