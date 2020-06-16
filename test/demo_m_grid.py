#!/usr/bin/env python3

'Démonstration pour grille de type M'

import os
import pathlib
import shutil
import numpy as np
import rpnpy.librmn.all as rmn

def main():
    'Ajoute champs sinus à fichier Standard RPN'

    pathlib.Path('out').mkdir(exist_ok=True)
    shutil.copyfile(os.path.join('in', 'mesh.std'), os.path.join('out', 'mesh.std'))
    funit = rmn.fstopenall(os.path.join('out', 'mesh.std'), rmn.FST_RW)
    lat = rmn.fstlir(funit, nomvar='^^')['d']
    lon = rmn.fstlir(funit, nomvar='>>')['d'][0]
    hb = rmn.fstlir(funit, nomvar='hb')
    x = hb
    x['nomvar'] = 'x'
    x['etiket'] = ''
    x['d'] = np.sin(np.pi*lon/180)*np.sin(np.pi*lat/90)
    rmn.fstecr(funit, x)
    rmn.fstcloseall(funit)

if __name__ == "__main__":

    main()
