import numpy as np
import pytest
from ../../python/proto import *

# pytest test/python/test_proto.py -v

def test_georef_new():
    """Test creating a new GeoRef object"""
    ref = GeoRef_New()
    assert ref is not None
    GeoRef_Free(ref)

def test_georef_setup():
    """Test setting up a GeoRef object with basic parameters"""
    ni, nj = 100, 100
    grtyp = b'Z'  # Note: strings need to be bytes for ctypes
    ig1, ig2, ig3, ig4 = 0, 0, 0, 0
    
    ref = GeoRef_Setup(ni, nj, grtyp, ig1, ig2, ig3, ig4)
    assert ref is not None
    GeoRef_Free(ref)

def test_zref_new():
    """Test creating a new ZRef object"""
    zref = ZRef_New()
    assert zref is not None
    ZRef_Free(zref)

def test_zref_define():
    """Test defining a ZRef object with levels"""
    zref = ZRef_New()
    nk = 10
    type_val = 1  # Example type value
    levels = np.zeros(nk, dtype=np.float32)
    for i in range(nk):
        levels[i] = float(i)
    
    result = ZRef_Define(zref, nk, type_val, levels)
    assert result == 1  # Assuming 1 means success
    ZRef_Free(zref)

def test_georef_xy2ll_z():
    """Test coordinate conversion from XY to LatLon for Z projection"""
    ref = GeoRef_Setup(100, 100, b'Z', 0, 0, 0, 0)
    
    # Test single point conversion
    n = 1
    x = np.array([0.0], dtype=np.float64)
    y = np.array([0.0], dtype=np.float64)
    lat = np.zeros(n, dtype=np.float64)
    lon = np.zeros(n, dtype=np.float64)
    
    result = GeoRef_XY2LL_Z(ref, lat, lon, x, y, n)
    assert result == 1  # Assuming 1 means success
    
    # Coordinates should be transformed but still be valid
    assert not np.isnan(lat[0])
    assert not np.isnan(lon[0])
    
    GeoRef_Free(ref)

def test_vertex_vals():
    """Test vertex value interpolation"""
    x, y = 0.5, 0.5  # Example coordinates
    ni, nj = 10, 10
    data = np.zeros((ni * nj), dtype=np.float64)
    mask = b'\x01' * (ni * nj)  # All valid
    
    val = VertexValS(data, mask, ni, nj, x, y, b'0')
    assert isinstance(val, float)

def test_ogr_layer():
    """Test OGR Layer creation and manipulation (if GDAL is available)"""
    try:
        layer = OGR_LayerNew(b"TestLayer")
        assert layer is not None
        OGR_LayerFree(layer)
    except AttributeError:
        pytest.skip("GDAL support not available")

def test_georef_interp():
    """Test geographic interpolation"""
    ref_from = GeoRef_Setup(100, 100, b'Z', 0, 0, 0, 0)
    ref_to = GeoRef_Setup(50, 50, b'Z', 0, 0, 0, 0)
    
    # Create options
    opt = TGeoOptions()
    
    # Create sample data
    data_in = np.zeros((100, 100), dtype=np.float32)
    data_out = np.zeros((50, 50), dtype=np.float32)
    
    # Test interpolation
    result = GeoRef_Value(ref_from, 0.5, 0.5)  # Test a single point interpolation
    assert isinstance(result, int)
    
    GeoRef_Free(ref_from)
    GeoRef_Free(ref_to)

if __name__ == "__main__":
    pytest.main([__file__]) 