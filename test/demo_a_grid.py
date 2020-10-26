#!/usr/bin/env python3

'Démonstration pour grille de type A'

import argparse
import os
from pathlib import Path
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

    params['nomvar'] = 'TT'
    params['grtyp'] = 'A'
    params['ni'] = nlon
    params['nj'] = nlat
    hemisphere_xg1 = {'global': 0, 'north': 1, 'south': 2}

    params['EZ_IG1'], params['ig2'], params['ig3'], params['ig4'] = \
        rmn.cxgaig(params['grtyp'], hemisphere_xg1[hemisphere])

    return params

def plot_grid(params):
    'Display grid on map'

    gid = rmn.ezqkdef(params)
    lalo = rmn.gdll(gid)

    axes = plt.axes(projection=ccrs.PlateCarree())
    axes.coastlines()
    axes.set_global()

    lalo['lon'] = np.mod((lalo['lon'] + 180), 360) - 180

    plt.plot(lalo['lon'], lalo['lat'],
             linestyle='None', color='red', marker='.')

    plt.savefig(os.path.join('out', 'a_grid.png'))

def plot_data(params, data):
    'Plot data on map'

    gid = rmn.ezqkdef(params)
    lalo = rmn.gdll(gid)

    fig = plt.figure()
    axes = fig.gca(projection='3d')

    lon = lalo['lon'].flatten()
    lat = lalo['lat'].flatten()

    axes.plot(lon, lat, data.flatten(), linestyle='None', marker='.')

    plt.savefig(os.path.join('out', 'a_data.png'))

def error(params, data, georef=False):
    """Calculate error with respect to analytical truth.

    Parameters
    ----------
    params : dict
        Produce with gen_data
    data : array_like
        Data on source grid
    georef : bool, optional
        Use georef for interpolation. librmn used otherwise.

    Returns
    -------
    int
        norm of difference between interpolated and true data

    """

    nlat, nlon = (9, 6)
    lat0, lon0, dlat, dlon = (-80, 0, 20, 30)

    if georef:
        gid = libgeoref.GeoRef_RPNCreate(params)
        out_gid = libgeoref.defGrid_L(nlon, nlat, lat0, lon0, dlat, dlon)
        out_data = libgeoref.GeoRef_Interp(out_gid, gid, data)
        out_lalo = libgeoref.gdll(out_gid)
    else:
        gid = rmn.ezqkdef(params)
        out_gid = rmn.defGrid_L(nlon, nlat, lat0, lon0, dlat, dlon)
        out_data = rmn.ezsint(out_gid, gid, data)
        out_lalo = rmn.gdll(out_gid)

    true_data = np.sin(np.pi*out_lalo['lon']/180)\
        *np.sin(np.pi*out_lalo['lat']/90)
    return np.linalg.norm(out_data - true_data)

def main(georef=False):
    """Call all functions in order.

    Parameter
    ---------
    georef : bool, optional
        Use georef for interpolation. librmn used otherwise.

    """

    Path('out').mkdir(exist_ok=True)

    nlon = 360//45
    nlat = 180//30
    hemisphere = 'global'
    params = gen_params(nlon, nlat, hemisphere)
    plot_grid(params)

    gid = rmn.ezqkdef(params)
    lalo = rmn.gdll(gid)
    data = np.sin(np.pi*lalo['lon']/180)*np.sin(np.pi*lalo['lat']/90)

    plot_data(params, data)

    stage_2020.write_fst(data, params, os.path.join('out', 'a_grid.fst'))
    print('Interpolation error: {:g}'.format(error(params, data, georef)))

if __name__ == "__main__":

    PARSER = argparse.ArgumentParser()
    PARSER.add_argument('-g', '--georef', action='store_true',
                        help='use libgeoref')
    ARGS = PARSER.parse_args()

    if 'DISPLAY' not in os.environ:
        plt.switch_backend('agg')

    main(ARGS.georef)
