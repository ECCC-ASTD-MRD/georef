#!/bin/bash

# Test latlon
\rm /tmp/toto.fst; 
#./build/util/Interpolate -g  data/GRIDS/in/L.fstd -i $CMCGRIDF/prog/glbeta/2020091412_000 -o /tmp/toto.fst -n TT 
./build/util/Interpolate -g  data/GRIDS/in/L.fstd -i data/GRIDS/out/L_sinus -o /tmp/toto.fst -n XX
./build/util/Interpolate -g  data/GRIDS/in/L.fstd -i data/GRIDS/O.fstd -o /tmp/toto.fst -n TM 
scp /tmp/toto.fst nil000@ppp3:/home/nil000/