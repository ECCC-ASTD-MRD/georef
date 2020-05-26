#!/usr/bin/env python3

'Wrapper pour la fonction c_ezsint de libgeoref'

import ctypes as ct
import numpy  as _np
import rpnpy.librmn.base as rb
import rpnpy.librmn.interp as inter
from rpnpy.librmn  import const as rc
from rpnpy import C_WCHAR2CHAR as _C_WCHAR2CHAR
from rpnpy import C_CHAR2WCHAR as _C_CHAR2WCHAR
from rpnpy import C_MKSTR as _C_MKSTR
from rpnpy import integer_types as _integer_types
from . import proto as rp

_isftn = lambda x, t: x.dtype == t and x.flags['F_CONTIGUOUS']
_ftn = lambda x, t: x if _isftnf32(x) else _np.asfortranarray(x, dtype=t)
_isftnf32 = lambda x: _isftn(x, _np.float32)
_ftnf32 = lambda x: _ftn(x, _np.float32)
_ftnOrEmpty = lambda x, s, t: \
    _np.empty(s, dtype=t, order='F') if x is None else _ftn(x, t)

def wrapper_ezdefset(gdidout, gdidin):
    gdidout = inter._getCheckArg(int, gdidout, gdidout, 'id')
    gdidin = inter._getCheckArg(int, gdidin, gdidin, 'id')
    istat = rp.c_ezdefset(gdidout, gdidin)
    if istat < 0:
        raise inter.EzscintError()
    return istat

def wrapper_ezgprm(gdid, doSubGrid=False):
    gdid = inter._getCheckArg(int, gdid, gdid, 'id')
    (cni, cnj) = (ct.c_int(), ct.c_int())
    (cgrtyp, cig1, cig2, cig3, cig4) = (_C_MKSTR(' '*rc.FST_GRTYP_LEN),
                                        ct.c_int(), ct.c_int(), ct.c_int(),
                                        ct.c_int())
    istat = rp.c_ezgprm(gdid, cgrtyp, cni, cnj, cig1, cig2, cig3, cig4)
    if istat < 0:
        raise inter.EzscintError()
    params = {
        'id'    : gdid,
        'shape' : (max(1, cni.value), max(1, cnj.value)),
        'ni'    : cni.value,
        'nj'    : cnj.value,
        'grtyp' : _C_CHAR2WCHAR(cgrtyp.value),
        'ig1'   : cig1.value,
        'ig2'   : cig2.value,
        'ig3'   : cig3.value,
        'ig4'   : cig4.value
            }
    if doSubGrid:
        params['nsubgrids'] = inter.ezget_nsubgrids(gdid)
        params['subgridid'] = inter.ezget_subgridids(gdid)
        params['subgrid'] = []
        if params['nsubgrids'] > 0:
            for gid2 in params['subgridid']:
                params['subgrid'].append(wrapper_ezgprm(gid2))
    return params

def wrapper_ezgxprm(gdid, doSubGrid=False):
    gdid = inter._getCheckArg(int, gdid, gdid, 'id')
    (cni, cnj) = (ct.c_int(), ct.c_int())
    cgrtyp = _C_MKSTR(' '*rc.FST_GRTYP_LEN)
    (cig1, cig2, cig3, cig4) = (ct.c_int(), ct.c_int(),
                                ct.c_int(), ct.c_int())
    cgrref = _C_MKSTR(' '*rc.FST_GRTYP_LEN)
    (cig1ref, cig2ref, cig3ref, cig4ref) = (ct.c_int(), ct.c_int(),
                                            ct.c_int(), ct.c_int())
    istat = rp.c_ezgxprm(gdid, cni, cnj, cgrtyp, cig1, cig2, cig3, cig4,
                         cgrref, cig1ref, cig2ref, cig3ref, cig4ref)
    if istat < 0:
        raise inter.EzscintError()
    params = {
        'id'    : gdid,
        'shape' : (max(1, cni.value), max(1, cnj.value)),
        'ni'    : cni.value,
        'nj'    : cnj.value,
        'grtyp' : _C_CHAR2WCHAR(cgrtyp.value),
        'ig1'   : cig1.value,
        'ig2'   : cig2.value,
        'ig3'   : cig3.value,
        'ig4'   : cig4.value,
        'grref' : _C_CHAR2WCHAR(cgrref.value),
        'ig1ref'   : cig1ref.value,
        'ig2ref'   : cig2ref.value,
        'ig3ref'   : cig3ref.value,
        'ig4ref'   : cig4ref.value
            }
    if doSubGrid:
        params['nsubgrids'] = inter.ezget_nsubgrids(gdid)
        params['subgridid'] = inter.ezget_subgridids(gdid)
        params['subgrid'] = []
        if params['nsubgrids'] > 0:
            for gid2 in params['subgridid']:
                params['subgrid'].append(wrapper_ezgxprm(gid2))
    return params

def wrapper_ezqkdef(ni, nj=None, grtyp=None, ig1=None, ig2=None, ig3=None, ig4=None,
                    iunit=0):
    if isinstance(ni, dict):
        gridParams = ni
        try:
            (ni, nj) = gridParams['shape']
        except:
            (ni, nj) = (None, None)
        try:
            if not ni:
                ni = gridParams['ni']
            if not nj:
                nj = gridParams['nj']
            grtyp = gridParams['grtyp']
            ig1 = gridParams['ig1']
            ig2 = gridParams['ig2']
            ig3 = gridParams['ig3']
            ig4 = gridParams['ig4']
        except:
            raise TypeError('ezgdef_fmem: provided incomplete grid description')
        try:
            iunit = gridParams['iunit']
        except:
            iunit = 0
    if (type(ni), type(nj), type(grtyp), type(ig1), type(ig2), type(ig3),
            type(ig4), type(iunit)) != (int, int, str, int, int, int, int, int):
        raise TypeError('ezqkdef: wrong input data type')
    if grtyp.strip() in ('', 'X'):
        raise inter.EzscintError('ezqkdef: Grid type {0} Not supported'.format(grtyp))
    if iunit <= 0 and grtyp.strip() in ('Z', '#', 'Y', 'U'):
        raise inter.EzscintError('ezqkdef: A valid opened file unit ({0}) is needed for Grid type {1}'.format(iunit, grtyp))
    gdid = rp.c_ezqkdef(ni, nj, _C_WCHAR2CHAR(grtyp), ig1, ig2, ig3, ig4, iunit)
    if gdid >= 0:
        return gdid
    raise inter.EzscintError()

def wrapper_ezsint(gdidout, gdidin, zin, zout=None):
    gdidout = inter._getCheckArg(int, gdidout, gdidout, 'id')
    gdidin  = inter._getCheckArg(int, gdidin, gdidin, 'id')
    zin     = inter._getCheckArg(_np.ndarray, zin, zin, 'd')
    zout    = inter._getCheckArg(None, zout, zout, 'd')
    gridsetid = wrapper_ezdefset(gdidout, gdidin)
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
    istat = rp.c_ezsint(zout, zin)
    if istat >= 0:
        return zout
    raise inter.EzscintError()

def defGrid_L(ni, nj=None, lat0=None, lon0=None, dlat=None, dlon=None,
              setGridId=True):
    params = {
        'ni'   : ni,
        'nj'   : nj,
        'lat0' : lat0,
        'lon0' : lon0,
        'dlat' : dlat,
        'dlon' : dlon
        }
    if isinstance(ni, dict):
        params.update(ni)
        try:
            setGridId = ni['setGridId']
        except:
            pass
    for k in ('ni', 'nj'):
        v = params[k]
        if not isinstance(v, _integer_types):
            raise TypeError('defGrid_L: wrong input data type for ' +
                            '{0}, expecting int, Got ({1})'.format(k, type(v)))
        if v <= 0:
            raise ValueError('defGrid_L: grid dims must be >= 0, got {0}={1}'.format(k, v))
    for k in ('lat0', 'lon0', 'dlat', 'dlon'):
        v = params[k]
        if isinstance(v, _integer_types):
            v = float(v)
        if not isinstance(v, float):
            raise TypeError('defGrid_L: wrong input data type for ' +
                            '{0}, expecting float, Got ({1})'.format(k, type(v)))
        params[k] = v
    params['grtyp'] = 'L'
    ig1234 = rb.cxgaig(params['grtyp'], params['lat0'], params['lon0'],
                        params['dlat'], params['dlon'])
    params['ig1'] = ig1234[0]
    params['ig2'] = ig1234[1]
    params['ig3'] = ig1234[2]
    params['ig4'] = ig1234[3]
    params['id'] = wrapper_ezqkdef(params) if setGridId else -1
    params['shape'] = (params['ni'], params['nj'])
    return params
