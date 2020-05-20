#!/usr/bin/env python3

'Démonstration pour grille de type Orca'

import os
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import rpnpy.librmn.all as rmn

def plot_grid(lat, lon):
    'Display grid on map'

    axes = plt.axes(projection=ccrs.PlateCarree())
    axes.set_global()
    axes.coastlines()

    # Downsample from 0.25 to 10°
    step = 40

    plt.plot(lon[::step, ::step], lat[::step, ::step], linestyle='None', color='black', marker='.')
    plt.savefig('o_grid.png')

def plot_data(lat, lon, data):
    'Plot data on map'

    fig = plt.figure()
    axes = fig.gca(projection='3d')

    lon = lon.flatten()
    lat = lat.flatten()

    axes.plot(lon, lat, data.flatten(), linestyle='None', marker='.')

    plt.savefig('o_data.png')

def write_fst(lat_record, lon_record, data):
    'Write data to RPN Standard file'

    funit_out = rmn.fstopenall('O_sinus.fst', rmn.FST_RW)
    rmn.fstecr(funit_out, lat_record)
    rmn.fstecr(funit_out, lon_record)
    params = rmn.FST_RDE_META_DEFAULT
    params['nomvar'] = 'XX'
    params['grtyp'] = 'O'
    params['ni'] = data.shape[0]
    params['nj'] = data.shape[1]
    params['ig1'] = lat_record['ip1']
    params['ig2'] = lat_record['ip2']
    params['ig3'] = lat_record['ip3']
    rmn.fstecr(funit_out, data, params)
    rmn.fstcloseall(funit_out)

def error_cstintrp():
    'Calculate error with respect to analytical truth'

    funit = rmn.fstopenall(os.path.join('GRIDS', 'L_sinus'))
    out_data = rmn.fstlir(funit, nomvar='XX', typvar='P@')['d']
    rmn.fstcloseall(funit)

    nlat, nlon = (9, 6)
    lat0, lon0, dlat, dlon = (-80, 0, 20, 30)
    out_gid = rmn.defGrid_L(nlon, nlat, lat0, lon0, dlat, dlon)
    out_lalo = rmn.gdll(out_gid)

    true_data = np.sin(np.pi*out_lalo['lon']/180)\
        *np.sin(np.pi*out_lalo['lat']/90)
    difference = out_data - true_data
    np.savetxt('dif.txt', difference)
    return np.linalg.norm(difference)

def main():
    'Call all functions in order'

    # Read lat-lon points
    funit = rmn.fstopenall(os.path.join('..', '..', 'GRIDS', 'O.fstd'))
    lat_record = rmn.fstlir(funit, nomvar='^^', dtype=np.float32)
    lon_record = rmn.fstlir(funit, nomvar='>>', dtype=np.float32)
    lat = lat_record['d']
    lon = lon_record['d']
    rmn.fstcloseall(funit)

    plot_grid(lat, lon)

    step = 35
    data = np.sin(np.pi*lon/180)*np.sin(np.pi*lat/90)
    plot_data(lat[::step, ::step], lon[::step, ::step], data[::step, ::step])

    write_fst(lat_record, lon_record, data)

    print(error_cstintrp())

if __name__ == "__main__":
    main()
