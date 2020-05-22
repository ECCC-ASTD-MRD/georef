#!/usr/bin/env python3

'Comparaison entre différentes interpolations de la grille ORCA vers la grille L'

import os
from pathlib import Path
import numpy as np
import rpnpy.librmn.all as rmn

def main():
    'Call all functions in order'

    # Read XX field
    funit = rmn.fstopenall(os.path.join('GRIDS', 'out.spi'))
    xx_spi = rmn.fstlir(funit, nomvar='XX', typvar='P')['d']
    #xx_mask_spi = rmn.fstlir(funit, nomvar='XX', typvar='@@')['d']
    rmn.fstcloseall(funit)

    funit = rmn.fstopenall(os.path.join('GRIDS', 'L_sinus'))
    xx_cs = rmn.fstlir(funit, nomvar='XX', typvar='P@')['d']
    xx_mask_cs = rmn.fstlir(funit, nomvar='XX', typvar='@@')['d']
    rmn.fstcloseall(funit)

    funit = rmn.fstopenall(os.path.join('GRIDS', 'L_sinus_bilin'))
    xx_cs_bilin = rmn.fstlir(funit, nomvar='XX', typvar='P@')['d']
    xx_mask_cs_bilin = rmn.fstlir(funit, nomvar='XX', typvar='@@')['d']
    rmn.fstcloseall(funit)

    #xx_spi = np.ma.array(xx_spi, mask=np.logical_not(xx_mask_spi))
    #xx_spi = xx_spi.filled(fill_value=0)
    xx_cs = np.ma.array(xx_cs, mask=np.logical_not(xx_mask_cs))
    xx_cs = xx_cs.filled(fill_value=0)
    xx_cs_bilin = np.ma.array(xx_cs_bilin, mask=np.logical_not(xx_mask_cs_bilin))
    xx_cs_bilin = xx_cs_bilin.filled(fill_value=0)

    Path("out").mkdir(parents=True, exist_ok=True)

    difference = xx_spi - xx_cs
    error = np.linalg.norm(difference)
    np.savetxt(os.path.join('out', 'xx_cs.txt'), xx_cs)
    np.savetxt(os.path.join('out', 'dif.txt'), difference)
    print(error)

    difference = xx_spi - xx_cs_bilin
    error = np.linalg.norm(difference)
    np.savetxt(os.path.join('out', 'dif_bilin.txt'), difference)
    print(error)

if __name__ == "__main__":
    main()
