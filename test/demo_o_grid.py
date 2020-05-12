#!/usr/bin/env python3

'Démonstration pour grille de type Orca'

import os
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np
import rpnpy.librmn.all as rmn
from mpl_toolkits.mplot3d import Axes3D

def plot_grid(lat, lon):
    'Display grid on map'

    axes = plt.axes(projection=ccrs.PlateCarree())
    axes.set_global()
    axes.coastlines()

    # Downsample from 0.25 to 10°
    STEP = 40

    plt.plot(lon[::STEP, ::STEP], lat[::STEP, ::STEP], linestyle='None', color='black', marker='.')
    plt.savefig('o_grid.png')

def plot_data(lat, lon, data):
    'Plot data on map'

    fig = plt.figure()
    axes = fig.gca(projection='3d')

    lon = lon.flatten()
    lat = lat.flatten()

    axes.plot(lon, lat, data.flatten(), linestyle='None', marker='.')

    plt.savefig('o_data.png')

def main():
    # Read lat-lon points
    #funit = rmn.fstopenall(os.path.join('..', 'GRIDS', 'O.fstd'))
    funit = rmn.fstopenall(os.path.join('..', '..', 'O.fstd'))
    lat = rmn.fstlir(funit, nomvar='^^', dtype=np.float32)['d']
    lon = rmn.fstlir(funit, nomvar='>>', dtype=np.float32)['d']
    rmn.fstcloseall(funit)

    plot_grid(lat, lon)

    STEP = 35
    lon = lon[::STEP, ::STEP]
    lat = lat[::STEP, ::STEP]
    data = np.sin(np.pi*lon/180)*np.sin(np.pi*lat/90)

    plot_data(lat, lon, data)

if __name__ == "__main__":
    main()
