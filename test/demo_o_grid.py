#!/usr/bin/env python3

'Démonstration pour grille de type Orca'

import os
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np
import rpnpy.librmn.all as rmn

funit = rmn.fstopenall(os.path.join('..', 'GRIDS', 'O.fstd'))
lat = rmn.fstlir(funit, nomvar='^^', dtype=np.float32)['d']
lon = rmn.fstlir(funit, nomvar='>>', dtype=np.float32)['d']
rmn.fstcloseall(funit)

axes = plt.axes(projection=ccrs.PlateCarree())
axes.set_global()
axes.coastlines()

# Downsample from 0.25 to 10°
STEP = 40

plt.plot(lon[::STEP, ::STEP], lat[::STEP, ::STEP], linestyle='None', color='black', marker='.')
plt.savefig('o_grid.png')
