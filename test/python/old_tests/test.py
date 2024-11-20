#!/usr/bin/python3

from ctypes import *
cdll.LoadLibrary('build/src/libgeoref.so.0.1.0')


def main():
    'Call all functions in order'

    nlon = 360//45
    nlat = 180//30
    hemisphere = 'global'
    #params = gen_params(nlon, nlat, hemisphere)
    #plot_grid(params)
    params = {}
    gid = libgeoref.ezqkdef(params)

if __name__ == "__main__":
    main()