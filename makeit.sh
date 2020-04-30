#!/bin/bash

`type pgsm &>/dev/null` || {
    . r.load.dot comm/eccc/all/opt/intelcomp/intelpsxe-cluster-19.0.3.199 \
      main/opt/openmpi/openmpi-3.1.2--hpcx-2.4.0-mofed-4.6--intel-19.0.3.199 \
      main/opt/openmpi-setup/openmpi-setup-0.2 \
      rpn/libs/19.6-beta \
      rpn/utils/19.6-beta \
      rpn/vgrid/6.5.b2
}

export EC_CMAKE_MODULE_PATH=`pwd`/modules

\rm -f -r build; mkdir build; cd build
cmake $CMAKE_COMP_FLAGS ..
make -j
