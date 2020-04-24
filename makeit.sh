#!/bin/bash

export EC_CMAKE_MODULE_PATH=~/Projects/libeerUtils/modules/

if [[ -n ${INTEL_LICENSE_FILE} ]]; then
   CMAKE_COMP_FLAGS="-DCMAKE_C_COMPILER=icc  -DCMAKE_CXX_COMPILER=icpc -DCMAKE_Fortran_COMPILER=ifort"
fi
if [[ -n ${CRAYPE_VERSION} ]]; then
   CMAKE_COMP_FLAGS="-DCMAKE_SYSTEM_NAME=CrayLinuxEnvironment"
fi

\rm -f -r build; mkdir build; cd build
cmake $CMAKE_COMP_FLAGS ..
make -j