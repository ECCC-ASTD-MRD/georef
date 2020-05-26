#!/usr/bin/env python3

'Wrapper pour la fonction c_ezsint de libgeoref'

import numpy  as _np
import rpnpy.librmn.interp as inter
from . import proto as rp

_isftn = lambda x, t: x.dtype == t and x.flags['F_CONTIGUOUS']
_ftn   = lambda x, t: x if _isftnf32(x) else _np.asfortranarray(x, dtype=t)
_isftnf32 = lambda x: _isftn(x, _np.float32)
_ftnf32   = lambda x: _ftn(x, _np.float32)
_ftnOrEmpty = lambda x, s, t: \
    _np.empty(s, dtype=t, order='F') if x is None else _ftn(x, t)

def wrapper_ezsint(gdidout, gdidin, zin, zout=None):
    gdidout = inter._getCheckArg(int, gdidout, gdidout, 'id')
    gdidin  = inter._getCheckArg(int, gdidin, gdidin, 'id')
    zin     = inter._getCheckArg(_np.ndarray, zin, zin, 'd')
    zout    = inter._getCheckArg(None, zout, zout, 'd')
    gridsetid = inter.ezdefset(gdidout, gdidin)
    gridParams = inter.ezgxprm(gdidin)
    zin  = _ftnf32(zin)
    if zin.shape != gridParams['shape']:
        raise TypeError("Provided zin array have inconsistent " +
                        "shape compared to the input grid")
    dshape = inter.ezgprm(gdidout)['shape']
    zout = _ftnOrEmpty(zout, dshape, zin.dtype)
    if not (isinstance(zout, _np.ndarray) and zout.shape == dshape):
        raise TypeError("Wrong type,shape for zout: {0}, {1}"\
                        .format(type(zout), repr(dshape)))
    #print(rp.c_ezsint)
    istat = rp.c_ezsint(zout, zin)
    if istat >= 0:
        return zout
    else:
        print('Erreur')
    #raise inter.EzscintError()