"""Unit tests for the GeoOptions ctypes structure.

Focuses on the pointer-backed numpy-array properties (table, lutDef, ancilliary)
that follow the fst_record pattern: the C pointer field (underscore-prefixed) is
kept in sync with a numpy array held in __slots__, and the public property exposes
the array.  The tests verify that the C pointer actually points at the numpy data.

The tests are self-contained and need no external data files.

Run directly:
    python TestGeoOptions.py
or with pytest:
    pytest TestGeoOptions.py
"""

import ctypes

import numpy as np

from georef import GeoOptions


# --- helpers -----------------------------------------------------------------

def read_c_double_array(address, size):
    """Read a C double array at the given address and return it as a numpy array."""
    return np.asarray((ctypes.c_double * size).from_address(address))


def read_lut_rows(address, nrows, ncols):
    """Read a C double** (array of row pointers) and return a list of numpy rows."""
    row_ptrs = (ctypes.POINTER(ctypes.c_double) * nrows).from_address(address)
    rows = []
    for i in range(nrows):
        addr = ctypes.cast(row_ptrs[i], ctypes.c_void_p).value
        rows.append(read_c_double_array(addr, ncols))
    return rows


# --- tests -------------------------------------------------------------------

def test_defaults_and_kwargs():
    """Constructor starts from C defaults and keyword arguments override fields."""
    o = GeoOptions(PolarCorrect=False, Interp=2, NoData=0.0)
    assert o.PolarCorrect == 0
    assert o.Interp == 2
    assert o.NoData == 0.0
    # Default interpolation is IR_CUBIC (3)
    assert GeoOptions().Interp == 3


def test_table_property():
    """The table property wraps a 1D numpy array and points the C field at it."""
    o = GeoOptions()
    o.table = np.array([1.0, 2.0, 3.0])
    assert o.table is not None
    assert np.allclose(o.table, [1.0, 2.0, 3.0])
    # The C pointer must reference the same data
    assert np.allclose(read_c_double_array(o._table, 3), [1.0, 2.0, 3.0])


def test_ancilliary_property():
    """The ancilliary property wraps a 1D numpy array and points the C field at it."""
    o = GeoOptions()
    o.ancilliary = np.array([9.5, 8.5])
    assert o.ancilliary is not None
    assert np.allclose(o.ancilliary, [9.5, 8.5])
    assert np.allclose(read_c_double_array(o._ancilliary, 2), [9.5, 8.5])


def test_lutdef_property():
    """The lutDef property wraps a 2D numpy array as a C double** and sets lutSize/lutDim."""
    o = GeoOptions()
    lut = np.array([[0, 10, 20], [1, 11, 21], [2, 12, 22]], dtype=np.float64)
    o.lutDef = lut
    assert o.lutDef is not None
    assert np.allclose(o.lutDef, lut)
    assert o.lutSize == 3
    assert o.lutDim == 3
    # The C double** must resolve to the same rows
    rows = read_lut_rows(o._lutDef, 3, 3)
    for i in range(3):
        assert np.allclose(rows[i], lut[i]), f"row {i} mismatch"


def test_deletion():
    """Deleting a pointer property clears both the array and the C field."""
    o = GeoOptions()
    o.table = np.array([1.0, 2.0, 3.0])
    o.lutDef = np.array([[0, 1], [2, 3]], dtype=np.float64)

    del o.table
    assert o.table is None
    assert o._table is None

    del o.lutDef
    assert o.lutDef is None
    assert o.lutSize == 0
    assert o.lutDim == 0


def test_private_kwarg_rejected():
    """Passing an underscore-prefixed field as a keyword raises ValueError."""
    try:
        GeoOptions(_table=1)
    except ValueError:
        return
    raise AssertionError("expected ValueError for underscore-prefixed keyword")


def test_validation():
    """Non-array values and wrong dimensionality are rejected."""
    o = GeoOptions()
    try:
        o.table = [1.0, 2.0, 3.0]
    except TypeError:
        pass
    else:
        raise AssertionError("expected TypeError for non-array table")

    try:
        o.lutDef = np.array([1.0, 2.0, 3.0])
    except ValueError:
        pass
    else:
        raise AssertionError("expected ValueError for 1D lutDef")


def main():
    tests = [
        test_defaults_and_kwargs,
        test_table_property,
        test_ancilliary_property,
        test_lutdef_property,
        test_deletion,
        test_private_kwarg_rejected,
        test_validation,
    ]
    for t in tests:
        t()
        print(f"[PASS] {t.__name__}")
    print(f"\nAll {len(tests)} GeoOptions tests passed.")


if __name__ == "__main__":
    main()
