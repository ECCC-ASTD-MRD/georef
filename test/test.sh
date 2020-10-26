#!/bin/bash

# Test latlon
\rm /tmp/toto.fst; 
#./build/util/Interpolate -g  data/GRIDS/in/L.fstd -i $CMCGRIDF/prog/glbeta/2020091412_000 -o /tmp/libgeoeer.fst -n TT 
./build/util/Interpolate -g  data/GRIDS/in/L.fstd -i data/GRIDS/out/L_sinus -o /tmp/libgeoeer.fst -n XX
./build/util/Interpolate -g  data/GRIDS/in/L.fstd -i data/GRIDS/O.fstd -o /tmp/libgeoeer.fst -n TM 
scp /tmp/libgeoeer.fst nil000@ppp3:/home/nil000/