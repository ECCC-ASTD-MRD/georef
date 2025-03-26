#ifndef GEOREF_TYPE_C_H
#define GEOREF_TYPE_C_H


//! \defgroup cubed_sphere Cubed-sphere Grid
//! The cubed-sphere panelization is as follows:
//! ```
//!        +----+
//!        | 4  |
//!   +----+----+----+----+
//!   | 3  | 0  | 1  | 2  |
//!   +----+----+----+----+
//!        | 5  |
//!        +----+
//! ```
//! The center of panel 0 has longitude `Lon0`, latitude `Lat0` and its central
//! meridian has a clockwise rotation of `Yaw0` from the vertical (all in radians)
//!
//! Each panel is divided into `NumElem` x `NumElem` elements, and each element contains
//! `Degree` x `Degree` grid points. The position of points within an element are determined by
//! Gauss-Legendre quadrature. Along one side of a panel, there are `NumElem` x `Degree` points.
//!
//! There are several coordinate systems used when referring to a cubed-sphere grid:
//! 
//! UV
//! These are local coordinates, identical on each panel, in the [-pi/4, pi/4] range.
//! They represent angles in radians. In these coordinates, every element has the exact
//! same size; however, the points within an element do not have a uniform spacing.
//!
//! XY
//! These are local coordinates, but are different for each panel. On the X axis, they
//! are in the range [-0.5, `NumAxisPoints` - 0.5], on the Y axis, [-0.5, 6 * `NumAxisPoints` - 0.5]
//! Each panel has a different slice of the Y range:
//!     Panel 0: [-0.5, `NumAxisPoints` - 0.5]
//!     Panel 1: [1 * `NumAxisPoints` - 0.5, 2 * `NumAxisPoints` - 0.5]
//!     Panel 2: [2 * `NumAxisPoints` - 0.5, 3 * `NumAxisPoints` - 0.5]
//!       ...
//! In XY coordinates, all grid points are equidistant (along each axis).
//! There is a piecewise linear mapping between XY and UV coordinates: Each X (or Y) segment between
//! contiguous grid points has a linear correspondance to a segment on the U (or V) axis.
//! When (X, Y) coordinates have integer value, they correspond to the index of the corresponding point in data arrays.
//!
//! Cartesian (3D)
//! These are mostly used globally, but they are also used locally when transforming between UV and global coordinates.
//! The origin is at the center of the sphere, with the Y axis pointing at the north pole, the Z axis at the
//! intersection of the equator and meridian 0, and the X axis at the intersection of the equator and meridian 90 (pi/2
//! radians):
//!```
//!      y
//!      ^
//!      |
//!      |----> x
//!     /
//!    z
//!```
//! Cartesian coordinates are needed to apply the rotation needed to map the cube on the sphere (according to
//! `Lon0`, `Lat0` and `Yaw0`), when converting from XY to Lat-Lon, and the sphere to the cube, when converting from
//! Lat-Lon to XY. Note that several cartesian points can map to the same lat-lon point.
//!
//! Latitude-longitude (Lat-Lon)
//! Used for displaying data, and as universal coordinate system to convert between various grids.
//!
//! When converting from Lat-lon to XY coordinates, the following transformations occur:
//!    lat-lon -> Cartesian (global) -> rotation -> cartesian (local) -> UV -> XY
//! and the (almost) reverse when converting from XY to lat-lon:
//!    XY -> UV -> cartesian (local) -> cartesian (global) -> rotation -> lat-lon

//! A rotation matrix that can be applied to a 3D (cartesian) point
//! \ingroup cubed_sphere
typedef struct {
    double mat[3][3];
} RotationParam;

//! Set of parameters that describe a cubed-sphere grid
//! \ingroup cubed_sphere
typedef struct {
    double Lon0; //!< Longitude of central point of panel 0 (radians)
    double Lat0; //!< Latitude of central point of panel 0 (radians)
    double Yaw0; //!< Clockwise rotation of central meridian of panel 0 with respect to true north (at equator) (radians)
    int32_t NumElem;        //!< How many elements there are along the side of a panel
    int32_t Degree;         //!< Degree of the discretization (number of grid points along one side of a single element)
    int32_t NumAxisPoints;  //!< How many grid points there are along the side of a panel (= `NumElem` * `Degree`)
    int32_t NumPanelPoints; //!< How many grid points there are in a panel (= `NumAxisPoints`^2)
    RotationParam LocalToGlobal; //!< Rotation matrix to transform local (panel) cartesian coordinates into global
    RotationParam GlobalToLocal; //!< Rotation matrix to transform global cartesian coordinates into local ones
} CubedSphereParams;

//! Decode an angle stored as an integer (IG, 24 bits) into its real value, in the range [-pi, pi[
static inline double decode_cs_angle(const int32_t angle24) {
    const double interval = M_2PI / 0x1000000;
    return ((angle24 & 0xffffff) - 0x800000) * interval;
}

#endif 
