#!/usr/bin/env python3

'Démonstration pour grille de type A'

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import rpnpy.librmn.all as rmn
import libgeoref.interp as libgeoref
import stage_2020

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

def plot_data(params, data):
    'Plot data on map'

    gid = rmn.ezqkdef(params)
    lalo = rmn.gdll(gid)

    fig = plt.figure()
    axes = fig.gca(projection='3d')

    lon = lalo['lon'].flatten()
    lat = lalo['lat'].flatten()

    axes.plot(lon, lat, data.flatten(), linestyle='None', marker='.')

    plt.savefig('a_data.png')

def error(params, data):
    'Calculate error with respect to analytical truth'

    #gid = libgeoref.wrapper_ezqkdef(params)
    gid = rmn.ezqkdef(params)

    nlat, nlon = (9, 6)
    lat0, lon0, dlat, dlon = (-80, 0, 20, 30)
    #out_gid = libgeoref.defGrid_L(nlon, nlat, lat0, lon0, dlat, dlon)
    out_gid = rmn.defGrid_L(nlon, nlat, lat0, lon0, dlat, dlon)
    out_lalo = rmn.gdll(out_gid)

    out_data = libgeoref.wrapper_ezsint(out_gid, gid, data)
    #out_data = rmn.ezsint(out_gid, gid, data)
    true_data = np.sin(np.pi*out_lalo['lon']/180)\
        *np.sin(np.pi*out_lalo['lat']/90)
    return np.linalg.norm(out_data - true_data)

def main():
    'Call all functions in order'

    nlon = 360//45
    nlat = 180//30
    hemisphere = 'global'
    params = gen_params(nlon, nlat, hemisphere)
    #plot_grid(params)

    # Fine mesh required to recognize analytic 3-D surface
    nlon = 360//5
    nlat = 180//5
    params = gen_params(nlon, nlat, hemisphere)

    gid = rmn.ezqkdef(params)
    lalo = rmn.gdll(gid)
    data = np.sin(np.pi*lalo['lon']/180)*np.sin(np.pi*lalo['lat']/90)

    #plot_data(params, data)

    #stage_2020.write_fst(data, params, 'a_grid.fst')
    #print(error(params, data))
    error(params, data)

if __name__ == "__main__":
    main()
