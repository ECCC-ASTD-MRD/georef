#!/usr/bin/env python3

'Démonstration pour grille de type A'

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np
import rpnpy.librmn.all as rmn
from mpl_toolkits.mplot3d import Axes3D

def gen_params(nlon, nlat, hemisphere):
    'Generate grid parameters'

    params = rmn.FST_RDE_META_DEFAULT

    params['grtyp'] = 'A'
    params['ni'] = nlon
    params['nj'] = nlat
    hemisphere_xg1 = {'global': 0, 'north': 1, 'south': 2}

    params['ig1'], params['ig2'], params['ig3'], params['ig4'] = \
        rmn.cxgaig(params['grtyp'], hemisphere_xg1[hemisphere])

    return params

def plot_grid(params):
    'Display grid on map'

    gid = rmn.ezqkdef(params)
    lalo = rmn.gdll(gid)

    axes = plt.axes(projection=ccrs.PlateCarree())
    axes.coastlines()
    axes.gridlines(xlocs=(-180, 0, 180), ylocs=(-90, 0, 90))
    axes.set_global()

    lalo['lon'] = np.mod((lalo['lon'] + 180), 360) - 180

    plt.plot(lalo['lon'], lalo['lat'],
             linestyle='None', color='red', marker='.')

    plt.savefig('a_grid.png')

def plot_data(params):
    'Plot data on map'

    gid = rmn.ezqkdef(params)
    lalo = rmn.gdll(gid)
    lalo['lon'] = np.mod((lalo['lon'] + 180), 360) - 180

    fig = plt.figure()
    ax = fig.gca(projection='3d')

    lon = lalo['lon'].flatten()
    lat = lalo['lat'].flatten()
    z = np.sin(np.pi*lon/180)*np.sin(np.pi*lat/90)

    ax.plot(lon, lat, z, linestyle='None', marker='.')

    plt.savefig('a_data.png')

def main():
    'Call all functions in order'

    nlon = 360//45
    nlat = 180//30
    hemisphere = 'global'
    params = gen_params(nlon, nlat, hemisphere)
    plot_grid(params)

    # Fine mesh required for 3-D surface
    nlon = 360//5
    nlat = 180//5
    params = gen_params(nlon, nlat, hemisphere)
    plot_data(params)

if __name__ == "__main__":
    main()
