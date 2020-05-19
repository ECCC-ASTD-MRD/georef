#!/usr/bin/env python3

'Comparaison entre différentes interpolations de la grille ORCA vers la grille L'

import os
import numpy as np
import rpnpy.librmn.all as rmn

def main():
    'Call all functions in order'

    # Read TM field
    funit = rmn.fstopenall(os.path.join('..', '..', 'GRIDS', 'out', 'out.spi'))
    tm_spi = rmn.fstlir(funit, nomvar='TM', typvar='P@')['d']
    tm_mask_spi = rmn.fstlir(funit, nomvar='TM', typvar='@@')['d']
    rmn.fstcloseall(funit)

    funit = rmn.fstopenall(os.path.join('..', '..', 'GRIDS', 'out', 'out.csintrp'))
    tm_cs = rmn.fstlir(funit, nomvar='TM', typvar='P@')['d']
    tm_mask_cs = rmn.fstlir(funit, nomvar='TM', typvar='@@')['d']
    rmn.fstcloseall(funit)

    tm_spi = np.ma.array(tm_spi, mask=np.logical_not(tm_mask_spi))
    tm_spi = tm_spi.filled(fill_value=0)
    tm_cs = np.ma.array(tm_cs, mask=np.logical_not(tm_mask_cs))
    tm_cs = tm_cs.filled(fill_value=0)

    difference = tm_spi - tm_cs
    error = np.linalg.norm(difference)
    np.savetxt('dif.txt', difference)
    print(error)

if __name__ == "__main__":
    main()
