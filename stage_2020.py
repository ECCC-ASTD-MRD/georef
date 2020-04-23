'Fonctions commune à multiples exemples'

import rpnpy.librmn.all as rmn

def write_fst(data, params, filename):
    'Write data to RPN Standard file'

    funit_out = rmn.fstopenall(filename, rmn.FST_RW)
    rmn.fstecr(funit_out, data, params)
    rmn.fstcloseall(funit_out)
