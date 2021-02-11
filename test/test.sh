#!/bin/bash

GRID_PATH=/local/drive1/afsr005/GRIDS/
GRIDS_NEW="O M ZW"
GRIDS="L PS PSN PSS ZE ZL O ZW M"
GRID=L

rm /tmp/libgeoref.fstd

for grid in $GRIDS; do
   echo "Interpolating $grid to $GRID"
   ./build/util/Interpolate -g ${GRID_PATH}/${GRID}.fstd -i ${GRID_PATH}/${grid}.fstd -o /tmp/libgeoref.fstd -n GRID -v DEBUG
done

scp /tmp/libgeoref.fstd nil000@ppp3:/home/nil000/

cd.; \rm /tmp/libgeoref.fstd; make test; scp /tmp/libgeoref.fstd nil000@ppp3:/home/nil000