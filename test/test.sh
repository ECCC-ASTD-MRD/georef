#!/bin/bash

GRID_PATH=/users/dor/afsr/005/Projects/RPN/libgeoref/GRIDS
GRIDS_NEW="O M ZW"
GRIDS="L PS PSN PSS ZE ZL O ZW M"
GRID=L

rm /tmp/libgeoref.fstd

for grid in $GRIDS; do
   echo "Interpolating $grid to $GRID"
   ./build/util/Interpolate -g ${GRID_PATH}/${GRID}.fstd -i ${GRID_PATH}/${grid}.fstd -o /tmp/libgeoref.fstd -n GRID -v DEBUG
done

scp /tmp/libgeoref.fstd nil000@ppp3:/home/nil000/

# Test latlon
#./build/util/Interpolate -g  data/GRIDS/in/L.fstd -i $CMCGRIDF/prog/glbeta/2020091412_000 -o /tmp/libgeoeer.fst -n TT 
#./build/util/Interpolate -g  data/GRIDS/in/L.fstd -i data/GRIDS/out/L_sinus -o /tmp/libgeoeer.fst -n XX
#./build/util/Interpolate -g  data/GRIDS/in/L.fstd -i data/GRIDS/O.fstd -o /tmp/libgeoeer.fst -n TM 
