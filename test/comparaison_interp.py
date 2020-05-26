#!/usr/bin/env python3

'Comparaison entre différentes interpolations de la grille ORCA vers la grille L'

import os
import numpy as np
import rpnpy.librmn.all as rmn

def main():
    'Call all functions in order'

    # Read XX field
    funit = rmn.fstopenall(os.path.join('GRIDS', 'out.spi'))
    xx_spi = rmn.fstlir(funit, nomvar='XX', typvar='P')['d']
    rmn.fstcloseall(funit)

    funit = rmn.fstopenall(os.path.join('GRIDS', 'out.csintrp.avg'))
    xx_cs = rmn.fstlir(funit, nomvar='XX', typvar='P@')['d']
    xx_mask_cs = rmn.fstlir(funit, nomvar='XX', typvar='@@')['d']
    rmn.fstcloseall(funit)

    funit = rmn.fstopenall(os.path.join('GRIDS', 'out.csintrp'))
    xx_cs_bilin = rmn.fstlir(funit, nomvar='XX', typvar='P@')['d']
    xx_mask_cs_bilin = rmn.fstlir(funit, nomvar='XX', typvar='@@')['d']
    rmn.fstcloseall(funit)

    xx_cs = np.ma.array(xx_cs, mask=np.logical_not(xx_mask_cs))
    xx_cs = xx_cs.filled(fill_value=0)
    xx_cs_bilin = np.ma.array(xx_cs_bilin, mask=np.logical_not(xx_mask_cs_bilin))
    xx_cs_bilin = xx_cs_bilin.filled(fill_value=0)

    difference = xx_spi - xx_cs
    print(np.linalg.norm(difference))

    difference = xx_spi - xx_cs_bilin
    print(np.linalg.norm(difference))

if __name__ == "__main__":
    main()
