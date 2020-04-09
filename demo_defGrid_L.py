#!/usr/bin/env python3

# FIXME Conda requis car Cartopy pas installé sur gpscc2, référence :
# https://gitlab.science.gc.ca/hpc/support/issues/5
# export PATH=~kro001/miniconda3/bin:"$PATH"
# . ssmuse-sh -x comm/eccc/all/opt/intelcomp/intelpsxe-cluster-19.0.3.199 \
# -x rpn/libs/19.5 \
# -x rpn/MIG/ENV/rpnpy/2.1.0

'Démonstration de la fonction defGrid_L.'

import cartopy.crs as ccrs
import numpy as np
import matplotlib.pyplot as plt
import rpnpy.librmn.all as rmn

def gen_data():
    'Generate random data'

    (nlat, nlon) = (9, 6)
    (lat0, lon0, dlat, dlon) = (-80, 0, 20, 30)
    params = rmn.defGrid_L(nlon, nlat, lat0, lon0, dlat, dlon)

    rng = np.random.default_rng()
    data = np.array(rng.random((nlon, nlat)), order='F')

    return data, params

def write_fst(data, params, filename):
    'Write data to RPN Standard file'

    funit_out = rmn.fstopenall(filename, rmn.FST_RW)
    rmn.fstecr(funit_out, data, params)
    rmn.fstcloseall(funit_out)

def plot_grid(params):
    'Display grid on map'

    axes = plt.axes(projection=ccrs.PlateCarree())
    axes.set_global()
    axes.coastlines()

    lats = np.arange(params['lat0'], params['lat0'] + params['nj']*params['dlat'], params['dlat'])
    lons = np.arange(params['lon0'], params['lon0'] + params['ni']*params['dlon'], params['dlon'])
    lats, lons = np.meshgrid(lats, lons)
    plt.plot(lons, lats, linestyle='None', color='black', marker='.')

    plt.savefig('l_grid.png')

def main():
    'Demonstrate all functions by calling them in correct order'

    data, params = gen_data()
    write_fst(data, params, 'l_grid.fst')
    plot_grid(params)

if __name__ == "__main__":
    main()
