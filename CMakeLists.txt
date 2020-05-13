cmake_minimum_required(VERSION 3.16)

#----- Append EC specific module path
foreach(PATH $ENV{EC_CMAKE_MODULE_PATH})
   list(APPEND CMAKE_MODULE_PATH ${PATH})
endforeach()

include(ec_utils)   # Include EC specific cmake utils
ec_parse_version()  # Parse VERSION file
ec_getvar()         # Get EC specific var
ec_build_info()     # Generate build include file

project(${NAME} VERSION ${VERSION} DESCRIPTION "${DESCRIPTION}")

set(CMAKE_INSTALL_PREFIX "" CACHE PATH "..." FORCE)

enable_language(C)
enable_language(Fortran)

include(compiler_rules)
link_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib)

find_package(RMN ${RMN_REQ_VERSION} REQUIRED)
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

#----- Build config information script
install(CODE "execute_process(COMMAND sed -e \"s/CMAKE_CC/\\\"${CMAKE_C_COMPILER_ID} ${CMAKE_C_COMPILER_VERSION}\\\"/\" -e \"s/CMAKE_FTN/\\\"${CMAKE_Fortran_COMPILER_ID} ${CMAKE_Fortran_COMPILER_VERSION}\\\"/\" 
   -e \"s/CMAKE_NAME/${NAME}/\" -e \"s/CMAKE_VERSION/${VERSION}${STATE}/\" -e \"s/CMAKE_RMN/${RMN_VERSION}${RMN_STATE}/\" 
   -e \"s/CMAKE_GDAL/${GDAL_VERSION}/\" -e \"s/CMAKE_VGRID/${VGRID_VERSION}${VGRID_STATE}/\"
   ../config OUTPUT_FILE ${NAME}-config)")
install(PROGRAMS ${CMAKE_BINARY_DIR}/${NAME}-config DESTINATION bin)

