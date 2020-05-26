# This modules loads compiler presets for the current platform and handles
# ECCC's computing environment differently

# Input
#  COMPILER_SUITE : Lower case name of the compiler suite (gnu, intel, ...)

# Check to see if we are at ECCC
if(DEFINED ENV{EC_ARCH})
   message(STATUS "ECCC environment variable detected; ECCC presets")
   set(EC_COMPILER_PRESET_PATH "compiler_presets/ECCC/$ENV{BASE_ARCH}/${EC_COMPILER_SUITE}.cmake")
else()
   message(STATUS "Using default presets")
   set(EC_COMPILER_PRESET_PATH "compiler_presets/default/${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}/${EC_COMPILER_SUITE}.cmake")
endif()

message(STATUS "Loading preset ${EC_COMPILER_PRESET_PATH}^")
include("${CMAKE_CURRENT_LIST_DIR}/${EC_COMPILER_PRESET_PATH}")
