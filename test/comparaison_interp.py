#!/usr/bin/env python3

'Comparaison entre différentes interpolations de la grille ORCA vers la grille L'

import os
import numpy as np
import rpnpy.librmn.all as rmn

def gen_true_data():
    'Generate true sinus data for L grid'

    nlat, nlon = (9, 6)
    lat0, lon0, dlat, dlon = (-80, 0, 20, 30)
    out_gid = rmn.defGrid_L(nlon, nlat, lat0, lon0, dlat, dlon)
    out_lalo = rmn.gdll(out_gid)

    true_data = np.sin(np.pi*out_lalo['lon']/180)\
        *np.sin(np.pi*out_lalo['lat']/90)
    
    return true_data

def error(data, data_compared_with):
    'Calculate error'

    difference = data - data_compared_with
    return np.linalg.norm(difference)

def main():
    'Call all functions in order'

    # Read XX field
    funit = rmn.fstopenall(os.path.join('out', 'out.spi'))
    xx_spi = rmn.fstlir(funit, nomvar='XX', typvar='P')['d']
    rmn.fstcloseall(funit)

    funit = rmn.fstopenall(os.path.join('out', 'out.csintrp.avg'))
    xx_cs = rmn.fstlir(funit, nomvar='XX', typvar='P@')['d']
    rmn.fstcloseall(funit)

    funit = rmn.fstopenall(os.path.join('out', 'out.csintrp'))
    xx_cs_bilin = rmn.fstlir(funit, nomvar='XX', typvar='P@')['d']
    rmn.fstcloseall(funit)

    true_data = gen_true_data()

    # Compare with respect to analytical truth
    print('*** Comparaison avec les vraies valeurs de sinus pour la grille L ***')
    print('cstintrp: {}'.format(error(xx_cs, true_data)))
    print('cstintrp_bilin: {}'.format(error(xx_cs_bilin, true_data)))
    print('SPI: {}'.format(error(xx_spi, true_data)))

    # Compare SPI and cstintrp
    print('*** Comparaison entre SPI et cstintrp ***')
    print('SPI et cstintrp: {}'.format(error(xx_spi, xx_cs))) 
    print('SPI et cstintrp_bilin: {}'.format(error(xx_spi, xx_cs_bilin)))

if __name__ == "__main__":
    main()
