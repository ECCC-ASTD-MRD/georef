#!/bin/bash

GRIDS_NEW="O M ZW"
GRIDS="L PS PSN PSS ZE ZL O ZW M"
GRID=L

rm /tmp/libgeoref.fstd

for grid in $GRIDS; do
   echo "Interpolating $grid to $GRID"
   ./build/util/Interpolate -g $ENV{ECCI_DATA_DIR}/libgeoref/grids/${GRID}.fstd -i $ENV{ECCI_DATA_DIR}/libgeoref/grids/${grid}.fstd -o /tmp/libgeoref.fstd -n GRID -v DEBUG
done

scp /tmp/libgeoref.fstd nil000@ppp3:/home/nil000/

cd.; \rm /tmp/libgeoref.fstd; make test; scp /tmp/libgeoref.fstd nil000@ppp3:/home/nil000