#!/usr/bin/env python3

# . ssmuse-sh -x comm/eccc/all/opt/intelcomp/intelpsxe-cluster-19.0.3.199 \
# -x rpn/libs/19.5 \
# -x rpn/MIG/ENV/rpnpy/2.1.0

import rpnpy.librmn.all as rmn
import numpy as np

def gen_data():
    'Generate random data'

    (lat0, lon0, dlat, dlon) = (0.,180.,1.,0.5)
    (ni, nj) = (5, 3)
    params = rmn.defGrid_L(ni, nj, lat0, lon0, dlat, dlon)
    data = np.array(np.random.rand(ni, nj), order='F')

    return data, params
    
def write_fst(data, params, filename):
    'Write data to RPN Standard file'
    
    funitOut = rmn.fstopenall(filename, rmn.FST_RW)
    rmn.fstecr(funitOut, data, params)
    rmn.fstcloseall(funitOut)

def main():
    data, params = gen_data()
    write_fst(data, params, 'newfile.fst')

if __name__ == "__main__":
    main()
