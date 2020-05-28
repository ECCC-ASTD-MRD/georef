#!/usr/bin/env python3

'Démonstration de la fonction defGrid_L.'

import os
from pathlib import Path
import cartopy.crs as ccrs
import numpy as np
import matplotlib.pyplot as plt
import rpnpy.librmn.all as rmn
import stage_2020

def gen_data():
    'Generate data and parameters'

    (nlat, nlon) = (9, 6)
    (lat0, lon0, dlat, dlon) = (-80, 0, 20, 30)
    params = rmn.defGrid_L(nlon, nlat, lat0, lon0, dlat, dlon)

    lalo = rmn.gdll(params['id'])
    data = np.sin(np.pi*lalo['lon']/180)*np.sin(np.pi*lalo['lat']/90)

    return data, params

def plot_grid(params):
    'Display grid on map'

    axes = plt.axes(projection=ccrs.PlateCarree())
    axes.set_global()
    axes.coastlines()

    lats = np.arange(params['lat0'], params['lat0'] \
                     + params['nj']*params['dlat'], params['dlat'])
    lons = np.arange(params['lon0'], params['lon0'] \
                     + params['ni']*params['dlon'], params['dlon'])
    lats, lons = np.meshgrid(lats, lons)
    plt.plot(lons, lats, linestyle='None', color='black', marker='.')

    plt.savefig(os.path.join('out', 'l_grid.png'))

def main():
    'Call all functions in order'

    Path('out').mkdir(exist_ok=True)
    data, params = gen_data()
    stage_2020.write_fst(data, params, os.path.join('out', 'l_grid.fst'))
    plot_grid(params)

if __name__ == "__main__":
    main()
