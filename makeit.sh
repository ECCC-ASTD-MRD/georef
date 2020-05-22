#!/bin/bash

`type pgsm &>/dev/null` || {
    . r.load.dot comm/eccc/all/opt/intelcomp/intelpsxe-cluster-19.0.3.199 \
      main/opt/openmpi/openmpi-3.1.2--hpcx-2.4.0-mofed-4.6--intel-19.0.3.199 \
      main/opt/openmpi-setup/openmpi-setup-0.2 \
      rpn/libs/19.6-beta \
      rpn/utils/19.6-beta \
      rpn/vgrid/6.5.b2
}

#----- Parse VERSION file
while IFS= read -r line; do
   if [[ `expr index "$line" ":"` -ne 0 ]]; then
      var=$(echo ${line%%:*} | xargs)
      val=$(echo ${line#*:} | xargs)
      eval $var=\""$val"\"
   fi
done < VERSION

#----- SSM packaging stuff
if [[ -n $COMP_ARCH ]]; then
   SSM_COMP=-${COMP_ARCH}
fi
ORDENV_PLAT=${ORDENV_PLAT:-`uname -s`-`uname -m`}
   
SSM_VERSION=${VERSION}${SSM_COMP}
SSM_NAME=${NAME}_${SSM_VERSION}_${ORDENV_PLAT}

export EC_CMAKE_MODULE_PATH="`pwd`/modules;$EC_CMAKE_MODULE_PATH"
export DESTDIR=${SSM_DEV}/workspace/${SSM_NAME}

\rm -f -r build; mkdir build; cd build
cmake $CMAKE_COMP_FLAGS -DCMAKE_INSTALL_PREFIX=$DESTDIR ..
#exit
make -j
make install
