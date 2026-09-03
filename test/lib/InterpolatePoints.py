"""Unit tests for GeoRef.interpolate.

Covers both the grid-to-grid path (GeoRef_Interp) and the lat/lon point-sampling
path (GeoRef_LLVal) that interpolate() dispatches to when the target is given as a
(lat, lon) tuple of numpy arrays.

The tests are self-contained: they build an in-memory global lat/lon 'A' grid and a
smooth analytic field, so no external FST data files (or ECCI_DATA_DIR) are needed.

Run directly:
    python InterpolatePoints.py
or with pytest:
    pytest InterpolatePoints.py
"""

import numpy as np
from georef import GeoOptions, GeoRef

# --- helpers -----------------------------------------------------------------


def make_grid(ni=36, nj=18):
    """Create an in-memory global lat/lon 'A' grid (ig1=GRID_GLOBAL=0)."""
    return GeoRef(ni, nj, "A", 0, 0, 0, 0)


def analytic(lat, lon):
    """Smooth analytic field used as ground truth (lat/lon in degrees)."""
    lat_r = np.radians(lat)
    lon_r = np.radians(lon)
    return np.sin(lat_r) * np.cos(lon_r) + 0.5 * np.cos(2.0 * lat_r)


def make_field(ref):
    """Sample the analytic field on the grid nodes (float32, Fortran order)."""
    lat, lon = ref.getll()
    return np.asfortranarray(analytic(lat, lon).astype(np.float32))


# --- tests -------------------------------------------------------------------


def test_point_sampling_at_nodes():
    """interpolate(field, (lat, lon)) at grid nodes recovers the field values."""
    ref = make_grid()
    field = make_field(ref)
    lat, lon = ref.getll()

    result = ref.interpolate(field, (lat, lon))

    assert result.shape == (lat.size,), f"expected 1-D output of size {lat.size}, got {result.shape}"
    assert result.dtype == np.float32
    assert np.all(np.isfinite(result)), "interpolated values contain non-finite entries"

    expected = analytic(lat.reshape(-1), lon.reshape(-1))
    max_err = np.max(np.abs(result - expected))
    assert max_err < 1e-3, f"point sampling at nodes deviates too much (max err {max_err:.3e})"


def test_point_sampling_arbitrary_points():
    """interpolate(field, (lat, lon)) at arbitrary interior points matches the analytic field."""
    ref = make_grid()
    field = make_field(ref)

    # Interior points (avoid the poles / dateline for a fair linear-interpolation check)
    lat = np.array([10.0, -20.0, 30.0, 5.0, -45.0], dtype=np.float64)
    lon = np.array([45.0, 120.0, 200.0, 300.0, 90.0], dtype=np.float64)

    result = ref.interpolate(field, (lat, lon))

    assert result.shape == (lat.size,)
    expected = analytic(lat, lon)
    # Linear interpolation of a smooth field: allow a few percent error
    max_err = np.max(np.abs(result - expected))
    assert max_err < 0.05, f"point sampling at arbitrary points deviates too much (max err {max_err:.3e})"


def test_grid_to_grid():
    """The original grid-to-grid path still works and stays close to the analytic field."""
    src = make_grid(ni=36, nj=18)
    field = make_field(src)
    tgt = make_grid(ni=18, nj=9)

    options = GeoOptions(
        Interp=2,  # IR_LINEAR
        PolarCorrect=False,
    )

    result = src.interpolate(field, tgt, options)

    assert result.shape == tgt.shape, f"expected {tgt.shape}, got {result.shape}"
    assert result.dtype == np.float32
    assert np.all(np.isfinite(result)), "interpolated values contain non-finite entries"

    lat, lon = tgt.getll()
    expected = analytic(lat, lon)
    max_err = np.max(np.abs(result - expected))
    assert max_err < 0.05, f"grid-to-grid interpolation deviates too much (max err {max_err:.3e})"


def test_invalid_target_type():
    """A target that is neither a GeoRef nor a (lat, lon) tuple raises TypeError."""
    ref = make_grid()
    field = make_field(ref)
    try:
        ref.interpolate(field, "not a valid target")
    except TypeError:
        return
    raise AssertionError("expected TypeError for invalid target_ref")


def test_mismatched_lat_lon_shapes():
    """A (lat, lon) tuple with mismatched shapes raises ValueError."""
    ref = make_grid()
    field = make_field(ref)
    lat = np.array([10.0, 20.0], dtype=np.float64)
    lon = np.array([10.0, 20.0, 30.0], dtype=np.float64)
    try:
        ref.interpolate(field, (lat, lon))
    except ValueError:
        return
    raise AssertionError("expected ValueError for mismatched lat/lon shapes")


def test_non_array_points():
    """A (lat, lon) tuple containing non-array elements raises TypeError."""
    ref = make_grid()
    field = make_field(ref)
    try:
        ref.interpolate(field, (10.0, 20.0))
    except TypeError:
        return
    raise AssertionError("expected TypeError for non-array lat/lon")


def main():
    tests = [
        test_point_sampling_at_nodes,
        test_point_sampling_arbitrary_points,
        test_grid_to_grid,
        test_invalid_target_type,
        test_mismatched_lat_lon_shapes,
        test_non_array_points,
    ]
    for t in tests:
        t()
        print(f"[PASS] {t.__name__}")
    print(f"\nAll {len(tests)} interpolate tests passed.")


if __name__ == "__main__":
    main()
