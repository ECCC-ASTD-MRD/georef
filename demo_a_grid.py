#!/usr/bin/env python3

'Démonstration pour grille de type A'

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np
import rpnpy.librmn.all as rmn

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

    plt.plot(lalo['lon'], lalo['lat'], linestyle='None', color='red', marker='.')

    plt.savefig('a_grid.svg')

def main():
    'Call all functions in order'

    nlon = 360//45
    nlat = 180//30
    hemisphere = 'global'
    params = gen_params(nlon, nlat, hemisphere)
    plot_grid(params)

if __name__ == "__main__":
    main()
