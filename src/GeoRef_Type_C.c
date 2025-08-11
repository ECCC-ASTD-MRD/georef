#include <App.h>

#include "GeoRef.h"
#include "georef/GeoRef_Type_C.h"

static const double AXIS_MIN = -M_PI4;
static const double AXIS_MAX = M_PI4;
static const double AXIS_RANGE = M_PI2;

//!> Coordinates of Gauss-Legendre quadrature points of degree 1-8, in the interval [-1, 1]
//!> These positions describe angular coordinates
static const double QUAD_POINTS[8][8] = {
    { 0.0, NAN, NAN, NAN,  NAN, NAN, NAN, NAN},
    {-0.5773502691896257,  0.5773502691896257,  NAN, NAN, NAN, NAN,  NAN, NAN},
    {-0.7745966692414834,  0.0               ,  0.7745966692414834,  NAN, NAN, NAN, NAN, NAN},
    {-0.8611363115940526, -0.3399810435848563,  0.3399810435848563,  0.8611363115940526, NAN, NAN, NAN, NAN},
    {-0.9061798459386640, -0.5384693101056831,  0.0               ,  0.5384693101056831, 0.9061798459386640, NAN, NAN, NAN},
    {-0.9324695142031521, -0.6612093864662645, -0.2386191860831969,  0.2386191860831969, 0.6612093864662645, 0.9324695142031521, NAN, NAN},
    {-0.9491079123427585, -0.7415311855993945, -0.4058451513773972,  0.0               , 0.4058451513773972, 0.7415311855993945, 0.9491079123427585, NAN},
    {-0.9602898564975363, -0.7966664774136267, -0.5255324099163290, -0.1834346424956498, 0.1834346424956498, 0.5255324099163290, 0.7966664774136267, 0.9602898564975363}
};

typedef enum {
    PANEL_NONE = -1,
    PANEL0 = 0,
    PANEL1 = 1,
    PANEL2 = 2,
    PANEL3 = 3,
    PANEL4 = 4,
    PANEL5 = 5,
} PanelID;

//! Describes a grid point, identified by 2D integer coordinates + which panel it lies on.
typedef struct {
    int32_t x, y;
    PanelID p;
} CSPoint;

//! Position of a point on a cubed sphere grid in UV coordinates
typedef struct {
    double u, v;
    PanelID p;
} CoordUV;
typedef CoordUV CoordXY; //! Alias

//! Position of a point in 3D space in cartesian coordinates. It does not have to be on the surface of the Earth.
typedef struct {
    double x, y, z;
} Coord3D;

//! Position of a point in 2D space, usually in XY coordinates, but it could also be cartesian.
typedef struct {
    double x, y;
} Coord2D;


// static inline uint32_t encode_cs_angle(const double angle) {
//     // # Keep in [-pi/2, pi/2[ range
//     // while f64 >= (math.pi / 2):
//     //     print(f"Adjusting from {f64} to {f64 - math.pi}")
//     //     f64 -= math.pi
//     // while f64 < (-math.pi / 2):
//     //     print(f"Adjusting from {f64} to {f64 + math.pi}")
//     //     f64 += math.pi

//     // # Scale, then shift, then truncate to 24 bits
//     // return (round(f64 / _INTERVAL) + 0x800000) & 0xFFFFFF
//     double small_angle = angle;
//     while (small_angle >= M_PI2) 
// }

//! Decode the IG4 parameter into its 2 components
static inline void decode_cs_ig4(
    const int32_t ig4,      //!< Value of the IG4 parameters
    int32_t* num_elements,  //!< [out] Number of elements along the side of a panel
    int32_t* num_solpts     //!< [out] Number of solution points per element side (= order of discretization)
) {
    *num_elements = ig4 >> 7;
    *num_solpts = ig4 & 0x7f;
}

//! Encode the IG4 parameter from its 2 components
static inline int32_t encode_cs_ig4(
    const int32_t num_elements, //!< Number of elements along the side of a panel
    const int32_t num_solpts    //!< Number of solution points per element side (= order of discretization)
) {
    //! \return The corresponding IG4
    return ((num_elements & 0x1ffff) << 7) | (num_solpts & 0x3f);
}

//! Create a rotation matrix that composes three rotations in cartesian coordinates (1 around each axis), 
//! in this order: Z, X, Y
//!```
//!      y
//!      ^
//!      |
//!      |----> x
//!     /
//!    z
//!```
//! The origin is at the center of the Earth, the Y axis points to the north pole, the Z axis to
//! point (0, 0) in lon-lat, the X axis to point (pi/2, 0)
static inline RotationParam make_cs_rotation(
    const double lambda,    //!< Angle of counterclockwise rotation in radians around Y (from center to north pole)
    const double phi,       //!< Angle of clockwise rotation in radians around X (from center to lon/lat (pi/2, 0))
    const double alpha      //!< Angle of clockwise rotation in radians around Z (from center to lon/lat (0, 0))
) {
    RotationParam rot;
    rot.mat[0][0] =  cos(lambda) * cos(alpha) + sin(lambda) * sin(phi) * sin(alpha);
    rot.mat[0][1] =  cos(lambda) * sin(alpha) - sin(lambda) * sin(phi) * cos(alpha);
    rot.mat[0][2] =  sin(lambda) * cos(phi);
    rot.mat[1][0] = -cos(phi)    * sin(alpha);
    rot.mat[1][1] =  cos(phi)    * cos(alpha);
    rot.mat[1][2] =  sin(phi);
    rot.mat[2][0] = -sin(lambda) * cos(alpha) + cos(lambda) * sin(phi) * sin(alpha);
    rot.mat[2][1] = -sin(lambda) * sin(alpha) - cos(lambda) * sin(phi) * cos(alpha);
    rot.mat[2][2] =  cos(lambda) * cos(phi);

    return rot;
}

//! Invert the given rotation matrix (just a transpose)
static inline RotationParam invert_cs_rotation(const RotationParam rot) {
    RotationParam inverse;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            inverse.mat[i][j] = rot.mat[j][i];
        }
    }
    return inverse;
}

//! Apply a rotation to a 3D point (in cartesian coordinates)
static inline Coord3D apply_rotation(
    const RotationParam rot, //!< The rotation we want to apply
    const Coord3D pt         //!< The point we want to rotate
) {
    //! \return Position of the point in 3D space after being rotated
    Coord3D result;
    result.x = rot.mat[0][0] * pt.x + rot.mat[0][1] * pt.y + rot.mat[0][2] * pt.z;
    result.y = rot.mat[1][0] * pt.x + rot.mat[1][1] * pt.y + rot.mat[1][2] * pt.z;
    result.z = rot.mat[2][0] * pt.x + rot.mat[2][1] * pt.y + rot.mat[2][2] * pt.z;
    return result;
}

//! Convert a position from longitude-latitude to cartesian coordinates.
//! When looking at (equator, meridian 0) from the sky, the x axis is towards the right, 
//! the y axis towards the top, and the z axis comes out of the screen towards us
//!```
//!      y
//!      ^
//!      |
//!      |----> x
//!     /
//!    z
//!```
static inline Coord3D ll_to_cart(const double lon, const double lat) {
    return (Coord3D){
        .x = sin(lon) * cos(lat),
        .y = sin(lat),
        .z = cos(lon) * cos(lat),
    };
}

//! Retrieve the closest grid point located to the lower left of the given XY coordinates
static inline CSPoint lower_left(const CoordXY Coord) {
    return (CSPoint){
        .x = (int32_t)floor(Coord.u),
        .y = (int32_t)floor(Coord.v),
        .p = Coord.p
    };
}

//! Create a CSPoint from the given coordinates
static inline CSPoint make_xy(const int32_t x, const int32_t y, const PanelID p) {
    return (CSPoint){.x = x, .y = y, .p = p};
}

//! Determine the "local" point that corresponds to the given point, which may lie outside the border
//! of the panel used to describe it. If the point is outside the given panel, we determine in which panel
//! it actually is, and its local coordinates within that panel.
static inline CSPoint local_xy(
    const CSPoint Coord,        //!< The initial point, which may or may not have the correct panel
    const int32_t PanelWidth    //!< Number of grid points per panel side
) {
    const CSPoint error = (CSPoint) {.x = -1, .y = -1, .p = PANEL_NONE};

    if (Coord.x <= -PanelWidth || Coord.x > 2*PanelWidth || Coord.y <= -PanelWidth || Coord.y >= 2*PanelWidth ||
        Coord.p == PANEL_NONE) {
        Lib_Log(APP_LIBGEOREF, APP_ERROR, "%s: Point (%d, %d)[%d] is too far out! (Panel width %d points)\n",
                __func__, Coord.x, Coord.y, Coord.p, PanelWidth);
        return error;
    }

    CSPoint result = Coord;

    for (int i = 0; i < 2; i++) {
        // Do this twice, for x and y (they might switch during an iteration, so we have to check both, both times)
        if (result.x < 0) {
            switch (result.p)
            {
            case PANEL0: result = make_xy(result.x + PanelWidth, result.y, PANEL3); break;
            case PANEL1: result = make_xy(result.x + PanelWidth, result.y, PANEL0); break;
            case PANEL2: result = make_xy(result.x + PanelWidth, result.y, PANEL1); break;
            case PANEL3: result = make_xy(result.x + PanelWidth, result.y, PANEL2); break;
            case PANEL4: result = make_xy(PanelWidth - 1 - result.y, result.x + PanelWidth, PANEL3); break;
            case PANEL5: result = make_xy(result.y, -result.x - 1, PANEL3); break;
            default: break; // Should not be possible
            }
        }
        else if (result.x >= PanelWidth) {
            switch (result.p)
            {
            case PANEL0: result = make_xy(result.x - PanelWidth, result.y, PANEL1); break;
            case PANEL1: result = make_xy(result.x - PanelWidth, result.y, PANEL2); break;
            case PANEL2: result = make_xy(result.x - PanelWidth, result.y, PANEL3); break;
            case PANEL3: result = make_xy(result.x - PanelWidth, result.y, PANEL0); break;
            case PANEL4: result = make_xy(result.y, 2*PanelWidth - result.x - 1, PANEL1); break;
            case PANEL5: result = make_xy(PanelWidth - 1 - result.y, result.x - PanelWidth, PANEL1); break;
            default: break; // Should not be possible
            }
        }
        else if (result.y < 0) {
            switch (result.p)
            {
            case PANEL0: result = make_xy(result.x, result.y + PanelWidth, PANEL5); break;
            case PANEL1: result = make_xy(result.y + PanelWidth, PanelWidth - result.x - 1, PANEL5); break;
            case PANEL2: result = make_xy(PanelWidth - 1 - result.x, -(1 + result.y), PANEL5); break;
            case PANEL3: result = make_xy(-result.y - 1, result.x, PANEL5); break;
            case PANEL4: result = make_xy(result.x, result.y + PanelWidth, PANEL0); break;
            case PANEL5: result = make_xy(PanelWidth - 1 - result.x, -(1 + result.y), PANEL2); break;
            default: break; // Should not be possible
            }
        }
        else if (result.y >= PanelWidth) {
            switch (result.p)
            {
            case PANEL0: result = make_xy(result.x, result.y - PanelWidth, PANEL4); break;
            case PANEL1: result = make_xy(2*PanelWidth - result.y - 1, result.x, PANEL4); break;
            case PANEL2: result = make_xy(PanelWidth - 1 - result.x, 2*PanelWidth - result.y - 1, PANEL4); break;
            case PANEL3: result = make_xy(result.y - PanelWidth, PanelWidth - 1 - result.x, PANEL4); break;
            case PANEL4: result = make_xy(PanelWidth - 1 - result.x, 2*PanelWidth - result.y - 1, PANEL2); break;
            case PANEL5: result = make_xy(result.x, result.y - PanelWidth, PANEL0); break;
            default: break; // Should not be possible
            }
        }
        else {
            break; // Get out of the loop
        }
    }

    return result;
}

//! Compute angular distance between given point and lon-lat (0, 0), in radians
static inline double angular_dist_0(const double lon, const double lat) {
    const double f1 = cos(lon) * sin(lat);
    const double f2 = sin(lon);
    return fabs(atan2(
        sqrt(f1*f1 + f2*f2),
        cos(lon) * cos(lat)
    ));
}

//! Compute euclidian distance between two (lon, lat) points
static inline double ll_dist(const double lon_a, const double lat_a, const double lon_b, const double lat_b) {
    const Coord3D pa = ll_to_cart(lon_a, lat_a);
    const Coord3D pb = ll_to_cart(lon_b, lat_b);
    const double x = pb.x - pa.x;
    const double y = pb.y - pa.y;
    const double z = pb.z - pa.z;

    // Lib_Log(APP_LIBGEOREF, APP_VERBATIM, "dist b/w (%5.2f, %5.2f) and (%5.2f, %5.2f): "
    //     "pts (%5.2f, %5.2f, %5.2f),  (%5.2f, %5.2f, %5.2f) "
    //     "diff (%5.2f, %5.2f, %5.2f)\n",
    //     lon_a, lat_a, lon_b, lat_b, pa.x, pa.y, pa.z, pb.x, pb.y, pb.z, x, y, z
    // );

    const double dist = sqrt(x*x + y*y + z*z);
    return dist;
}

//! Compute cartesian distance between two 2D points
static inline double cart_dist(const double x1, const double y1, const double x2, const double y2) {
    const double x = x2 - x1;
    const double y = y2 - y1;
    return sqrt(x*x + y*y);
}

//! Compute the bilinear interpolation of the 4 given values. The values are assumed to be arranged in
//! counterclockwise order. `weight_1` is the weight given to points `a` and `d` in each of [a, b] and [d, c] pairs;
//! `weight_2` is the weight given to points `a` and `b` in each of [a, d] and [b, c] pairs.
//! Each weight should be in the range [0, 1] and its opposite "1 - weight" is applied to the other member of
//! each pair.
static inline double bilinear_interp(
    const double a,         //!< 1st point "lower left"
    const double b,         //!< 2nd point "lower right"
    const double c,         //!< 3rd point "upper right"
    const double d,         //!< 4th point "upper left"
    const double weight_1,  //!< Weight applied to a in [a, b] and d in [d, c]
    const double weight_2   //!< Weigth applied to a in [a, d] and b in [b, c]
) {
    const double mid_1 = a * weight_1 + b * (1.0 - weight_1);
    const double mid_2 = d * weight_1 + c * (1.0 - weight_1);
    const double result = mid_1 * weight_2 + mid_2 * (1.0 - weight_2);
    return result;
}

//! 32-bit version of bilinear_interp()
//! \sa bilinear_interp
static inline float bilinear_interp_32(
    const float a,
    const float b,
    const float c,
    const float d,
    const float alpha_1,
    const float alpha_2
) {
    const float mid_1 = a * alpha_1 + b * (1.0 - alpha_1);
    const float mid_2 = d * alpha_1 + c * (1.0 - alpha_1);
    const float result = mid_1 * alpha_2 + mid_2 * (1.0 - alpha_2);
    return result;
}

//! Convert given 3D cartesian coordinates to lon-lat coordinates
static inline Coord2D cart_to_ll(const double x, const double y, const double z) {
    const double xz = sqrt(x*x + z*z);
    return (Coord2D){
        .x = atan2(x, z),
        .y = atan2(y, xz),
    };
}

//! Convert given 3D cartesian point from describing a position on panel 0 to a position
//! on the target panel. Even if the resulting position can be considered global, it is
//! *before* any rotation is applied to the grid; this means the resulting position might
//! need to be rotated before being used for other calculations.
static inline Coord3D rotate_to_panel(const Coord3D pt, const PanelID target_panel) {
    switch (target_panel)
    {
    case PANEL0:
        return pt;
    case PANEL1:
        return (Coord3D){.x = pt.z, .y = pt.y, .z = -pt.x};
    case PANEL2:
        return (Coord3D){.x = -pt.x, .y = pt.y, .z = -pt.z};
    case PANEL3:
        return (Coord3D){.x = -pt.z, .y = pt.y, .z = pt.x};
    case PANEL4:
        return (Coord3D){.x = pt.x, .y = pt.z, .z = -pt.y};
    case PANEL5:
        return (Coord3D){.x = pt.x, .y = -pt.z, .z = pt.y};
    default:
        Lib_Log(APP_LIBGEOREF, APP_ERROR, "%s: Invalid panel ID %d\n", __func__, target_panel);
        return (Coord3D){.x = 0.0, .y = 0.0, .z = 0.0};
    }
}

//! Convert given 3D cartesian point from describing a position on `origin_panel` to a position
//! on panel 0. Even if the initial position could be considered global, it must be a position
//! on a non-rotated grid; this means a truly global position must first be "un-rotated" before
//! calling this function.
static inline Coord3D rotate_from_panel(const Coord3D pt, const PanelID origin_panel) {
    switch (origin_panel)
    {
    case PANEL0:
        return pt;
    case PANEL1:
        return (Coord3D){.x = -pt.z, .y = pt.y, .z = pt.x};
    case PANEL2:
        return (Coord3D){.x = -pt.x, .y = pt.y, .z = -pt.z};
    case PANEL3:
        return (Coord3D){.x = pt.z, .y = pt.y, .z = -pt.x};
    case PANEL4:
        return (Coord3D){.x = pt.x, .y = -pt.z, .z = pt.y};
    case PANEL5:
        return (Coord3D){.x = pt.x, .y = pt.z, .z = -pt.y};
    default:
        Lib_Log(APP_LIBGEOREF, APP_ERROR, "%s: Invalid panel ID %d\n", __func__, origin_panel);
        return (Coord3D){.x = 0.0, .y = 0.0, .z = 0.0};
    }
}

//! Find the panel to which the given 3D cartesian point belongs. To call this function, the position must
//! first be un-rotated from the global position.
//! For a cube of side length 2, centered at the origin
//!   panel 0 lies on the `z = 1` plane, with x and y in the [-1, 1] range
//!   panel 1 lies on the `x = 1` plane, with y and z in the [-1, 1] range
//!   panel 2 lies on the `z = -1` plane, with x and y in the [-1, 1] range
//!   panel 3 lies on the `x = -1` plane, with y and z in the [-1, 1] range
//!   panel 4 lies on the `y = 1` plane, with x and z in the [-1, 1] range
//!   panel 5 lies on the `y = -1` plane, with x and z in the [-1, 1] range
//!
//! More generally, there are 6 zones extending from the origin that are delimited by the 6 planes
//! `x = y`, `x = -y`, `x = z`, `x = -z`, `y = z`, `y = -z`. This function determines in which of these zones
//! a point falls.
static inline PanelID find_panel(const Coord3D pt) {

    if (pt.x >= pt.y) {
        // Panel 0, 1, 2, 5

        if (pt.x >= pt.z) {
            // Panel 1, 2, 5
            if (pt.x >= -pt.z) {
                // Panel 1, 5
                if (pt.x >= -pt.y) return PANEL1;
                else return PANEL5;
            }
            else {
                // Panel 2, 5
                if (pt.y >= pt.z) return PANEL2;
                else return PANEL5;
            }
        }
        else {
            // Panel 0, 5
            if (pt.y >= -pt.z) return PANEL0;
            else return PANEL5;
        }
    }
    else {
        // Panel 0, 2, 3, 4
        if (pt.x >= pt.z) {
            // Panel 2, 4
            if (pt.y >= -pt.z) return PANEL4;
            else return PANEL2;
        }
        else {
            // Panel 0, 3, 4
            if (pt.x >= -pt.z) {
                // Panel 0, 4
                if (pt.y > pt.z) return PANEL4;
                else return PANEL0;
            }
            else {
                // Panel 3, 4
                if (pt.x >= -pt.y) return PANEL4;
                else return PANEL3;
            }
        }
    }
}

//! Convert a point from UV to 3D cartesian coordinates (un-rotated grid)
static inline Coord3D uv_to_cart(const CoordUV pt) {
    Coord3D tmp = (Coord3D){
        .x = tan(pt.u),
        .y = tan(pt.v),
        .z = 1.0,
    };
    return rotate_to_panel(tmp, pt.p);
}

//! Convert a point from 3D cartesian (un-rotated) to UV coordinates
static inline CoordUV cart_to_uv(const Coord3D pt) {
    const PanelID p = find_panel(pt);
    const Coord3D tmp = rotate_from_panel(pt, p);
    return (CoordUV) {
        .u = atan2(tmp.x, tmp.z),
        .v = atan2(tmp.y, tmp.z),
        .p = p,
    };
}

//! Convert a point from (1D) X (grid-point numbering) to U (angular) coordinates
static inline double x_to_u(
    const double x,             //!< The coordinate we want to convert. Must be in range [-0.5, `num_points` - 0.5]
    const double* axis_points,  //!< The list of U (angular) coordinates of each axis of a grid panel
    const int32_t num_points    //!< How many grid points there are along each axis of a grid panel
) {
    // If outside the range of points, we extrapolate (linearly)
    if (x < 0.0) {
        return (axis_points[0] * (0.5+x) - AXIS_MIN * x) / 0.5;
    }
    else if (x > (num_points - 1)) {
        const double alpha = x - (num_points - 1);
        return (axis_points[num_points - 1] * (0.5-alpha) + AXIS_MAX * alpha) / 0.5;
    }

    // Within the range of points, we interpolate (linearly)
    double x_int;
    const double alpha = modf(x, &x_int);
    const int32_t low = (int32_t)x_int;
    const int32_t high = low + 1;
    return axis_points[low] * (1-alpha) + axis_points[high] * alpha;
}

//! Convert a point from local XY to local UV coordinates
//! The `y` value of the points indicates in which panel it falls:
//!   Panel 0: [-0.5, `num_points` - 0.5[
//!   Panel 1: ]1 * `num_points` - 0.5, 2 * `num_points` - 0.5[
//!   Panel 2: ]2 * `num_points` - 0.5, 3 * `num_points` - 0.5[
//!   Panel 3: ]3 * `num_points` - 0.5, 4 * `num_points` - 0.5[
//!   Panel 4: ]4 * `num_points` - 0.5, 5 * `num_points` - 0.5[
//!   Panel 5: ]5 * `num_points` - 0.5, 6 * `num_points` - 0.5[
//!   Panel 0 again: [6 * `num_points` - 0.5]
//! When y is exactly between two panels (except 5 and 0), we simply can't determine on which panel
//! it goes, so that coordinate is not valid.
static inline CoordUV xy_to_uv(
    const double x,     //!< X coord of the point
    const double y,     //!< Y coord of the point
    const double* axis_points,  //!< U (angular) coord of grid points along each axis of a panel
    const int32_t num_points    //!< Number of grid points along each axis of a panel
) {
    const CoordUV error = (CoordUV){.u = 0.0, .v = 0.0, .p = PANEL_NONE};
    if (x < -0.5 || x > num_points - 0.5 || y < -0.5 || y > 6*num_points - 0.5) {
        Lib_Log(APP_LIBGEOREF, APP_ERROR, "%s: (%f, %f) out of range ([%f, %f], [%f, %f])\n",
                __func__, x, y, -0.5, num_points - 0.5, -0.5, 6 * num_points - 0.5);
        return error;
    }

    // We accept when panel == 6 (when y is exactly 6 * `num_points` - 0.5), because the top of
    // panel 5 wraps around to the bottom of panel 0
    const int panel = (int)((y + 0.5) / num_points);
    const double local_y = y - (panel * num_points);

    // if (local_y <= -0.5 && panel > 0 && panel < 6) {
    //     Lib_Log(APP_LIBGEOREF, APP_ERROR, "%s: Coordinate y %f lies on a discontinuity within the valid Y range.\n",
    //         __func__, y);
    //     return error;
    // }

    CoordUV c;
    c.p = panel % 6; // accept panel == 6
    c.u = x_to_u(x, axis_points, num_points);
    c.v = x_to_u(local_y, axis_points, num_points);
    return c;
}

//! Convert from angular to pointwise coordinates *in a reference element*
//! Angular positions are in the [-1, 1] range
//! Pointwise positions are in the [-0.5, `degree` - 0.5] range
static inline double refu_to_refx(const double u, const int32_t degree) {
    const double* points = QUAD_POINTS[degree - 1];
    // If less than first point, interpolate between -0.5 and 1st point
    if (u <= points[0]) {
        return 0.5 * (u - points[0]) / (points[0] + 1.0);
    }

    // Find between which pair of points the position falls
    for (int i = 1; i < degree; i++) {
        if (u <= points[i]) {
            return (u - points[i]) / (points[i] - points[i-1]) + i;
        }
    }

    // The position is beyond the last point: interpolate between last point and `degree` - 0.5
    return degree - 1 + 0.5 * (u - points[degree - 1]) / (1.0 - points[degree-1]);
}

static inline Coord2D uv_to_xy(
    const CoordUV pt,
    const int32_t num_elem,
    const int32_t degree
) {
    const double elem_size = AXIS_RANGE / num_elem;
    Coord2D result;

    const double u_offset = pt.u - AXIS_MIN;
    const int32_t elem_x = u_offset / elem_size;
    const double local_u = u_offset - (elem_x * elem_size);
    const double ref_u = local_u / elem_size * 2.0 - 1.0;
    result.x = elem_x * degree + refu_to_refx(ref_u, degree);
    // Lib_Log(APP_LIBGEOREF, APP_WARNING,
    //     "%s: Point %.3f (offset %.3f) is in element %d, remainder %g, reference element u %g, ref_x %g\n",
    //     __func__, pt.u, u_offset, elem_x, local_u, ref_u, refu_to_refx(ref_u, degree));

    const double v_offset = pt.v - AXIS_MIN;
    const int32_t elem_y = v_offset / elem_size;
    const double local_v = v_offset - (elem_y * elem_size);
    const double ref_v = local_v / elem_size * 2.0 - 1.0;
    result.y = elem_y * degree + refu_to_refx(ref_v, degree);
    result.y += num_elem * degree * pt.p;
    
    return result;
}

//! Retrieve the index of point (X, Y) in the array of all points
static inline size_t xy_to_index(
    const CSPoint Coord,        //! The point whose index we want
    const int32_t PanelWidth    //! Number of points along each axis of a panel
) {
    // Convert to correct panel (if it goes beyond the boundaries of its given panel)
    const CSPoint local = local_xy(Coord, PanelWidth);
    if (local.p == PANEL_NONE) return (size_t)-1; // if error
    return (local.y + PanelWidth * (int32_t)local.p) * PanelWidth + local.x;
}

/*----------------------------------------------------------------------------
 * @brief  Transforms XY grid coordinates to LatLon for a C (cube-sphere) grid
 *    @param[in]  Ref     Georeference pointer
 *    @param[out] Lat     Latitude array
 *    @param[out] Lon     Longitude array
 *    @param[in]  X       X array
 *    @param[in]  Y       Y array
 *    @param[in]  Nb      Number of coordinates
 * 
 * The following transformations take place:
 *    XY -> UV -> cartesian (local) -> cartesian (global) -> rotation -> lat-lon
 *
 *    @return             Number of points that were converted
 */
int32_t GeoRef_XY2LL_C(TGeoRef *Ref, double *Lat, double *Lon, double *X, double *Y, int32_t Nb) {

    int32_t num_points = 0;
    const int32_t num_axis_points = Ref->CGrid->NumAxisPoints;
    const RotationParam rot = Ref->CGrid->LocalToGlobal;
    for (int i = 0; i < Nb; i++) {
        const CoordUV pt_face = xy_to_uv(X[i], Y[i], Ref->AX, num_axis_points);
        if (pt_face.p == PANEL_NONE) continue;
        const Coord3D pt_cube = uv_to_cart(pt_face);
        const Coord3D pt_global = apply_rotation(rot, pt_cube);
        const Coord2D ll = cart_to_ll(pt_global.x, pt_global.y, pt_global.z);
        Lon[i] = ll.x;
        Lat[i] = ll.y;
        num_points++;
    }

    return num_points;
}

/*----------------------------------------------------------------------------
 * @brief  Transforms LatLon coordinates to XY for a C (cube-sphere) grid
 *    @param[in]  Ref     Georeference pointer
 *    @param[out] X       X array
 *    @param[out] Y       Y array
 *    @param[in]  Lat     Latitude array
 *    @param[in]  Lon     Longitude array
 *    @param[in]  Nb      Number of coordinates
 *
 * The following transformations take place:
 *    lat-lon -> Cartesian (global) -> rotation -> cartesian (local) -> UV -> XY
 *
 *    @return             Number of converted points
 */
int32_t GeoRef_LL2XY_C(TGeoRef *Ref, double *X, double *Y, double *Lat, double *Lon, int32_t Nb) {
    int32_t num_points = 0;

    const RotationParam rot = Ref->CGrid->GlobalToLocal;
    for (int i = 0; i < Nb; i++) {
        if (Lat[i] > M_PI2 || Lat[i] < -M_PI2) continue; // Invalid latitude
        const Coord3D pt_global = ll_to_cart(Lon[i], Lat[i]);
        const Coord3D pt_cube = apply_rotation(rot, pt_global);
        const CoordUV pt_face = cart_to_uv(pt_cube);
        const Coord2D xy = uv_to_xy(pt_face, Ref->CGrid->NumElem, Ref->CGrid->Degree);
        X[i] = xy.x;
        Y[i] = xy.y;
        num_points++;
    }

    return num_points;
}


static inline double rel_diff(const double a, const double b) {
    double diff = a - b;
    if (fabs(a) > 1e-7) diff /= a;
    return fabs(diff);
}

static inline double rel_diff_rad(const double a, const double b) {
    double a2 = a;
    if (a2 < 0.0) a2 += M_2PI;
    double b2 = b;
    if (b2 < 0.0) b2 += M_2PI;
    double diff = a2 - b2;
    if (fabs(a2) > 1e-7) diff /= a2;
    return fabs(diff);
}

static int check_diffs(const double* ref_lon, const double* ref_lat, const double* comp_lon, const double* comp_lat,
                        const int num_points, const char* label, const int verbose) {
    int num_errors = 0;
    double total_error = 0.0;
    for (int i = 0; i < num_points; i++) {
        const double diff = ll_dist(ref_lon[i], ref_lat[i], comp_lon[i], comp_lat[i]);
        total_error += diff;
        if (diff > 3e-7) {
            num_errors++;
            if (verbose) {
                Lib_Log(APP_LIBGEOREF, APP_ERROR, "Error[%4d] %8.2e (%6.3f, %6.3f) vs (%6.3f, %6.3f)\n",
                        i, diff, ref_lon[i], ref_lat[i], comp_lon[i], comp_lat[i]);
            }
        }
    }
    const double error_pc = num_errors * 100.0 / num_points;

    const int log_lvl = num_errors > 0 ? APP_ERROR : APP_INFO;
    Lib_Log(APP_LIBGEOREF, log_lvl, "%s: [%s] Num errors %3d (%4.1f%%) avg dist %10.2e\n",
            __func__, label, num_errors, error_pc, total_error / num_points);
    return num_errors;
}

//! Create an array with the coordinates of every panel point along the U (or V) axis. These coordinates
//! are the same on either axis. They vary between -pi/4 and pi/4. Elements have the same size in these
//! coordinates, but the points within an element have varying spacings.
static double* make_axis(
    const int32_t num_elem, //!< Number of elements along a panel side
    const int32_t degree    //!< Degree of discretization (= number of points in an element side)
) {
    //! \return A pointer to the array, with `num_elem` * `degree` entries. Will  need to be freed by the caller.
    double* axis = (double*) malloc(num_elem * degree * sizeof(double));

    const double elem_size = M_PI2 / num_elem;

    // UV coordinates of the first element
    // reference points are in the range [-1, 1], we need to shift+scale them to the first element
    double elem_base[degree];
    for (int i = 0; i < degree; i++) {
        elem_base[i] = (QUAD_POINTS[degree - 1][i] + 1.0) / 2.0 * elem_size;
    }

    // Translate to other elements
    for (int i = 0; i < num_elem; i++) {
        const double elem_start = -M_PI4 + elem_size * i;
        for (int j = 0; j < degree; j++) {
            axis[i*degree + j] = elem_start + elem_base[j];
        }
    }

    return axis;
}

//! Initialize data structures used to facilitate computations on a cubed-sphere grid
TGeoRef *GeoRef_DefineC(TGeoRef *Ref) {
    
    Ref->CGrid = (CubedSphereParams*)malloc(sizeof(CubedSphereParams));
    CubedSphereParams* param = Ref->CGrid;

    param->Lon0 = decode_cs_angle(Ref->RPNHead.ig1);
    param->Lat0 = decode_cs_angle(Ref->RPNHead.ig2);
    param->Yaw0 = decode_cs_angle(Ref->RPNHead.ig3);
    decode_cs_ig4(Ref->RPNHead.ig4, &(param->NumElem), &(param->Degree));

    Lib_Log(APP_LIBGEOREF, APP_WARNING, "%s: ig1-3 = %d (%g), %d (%g), %d (%g), num elem %d, degree %d\n",
            __func__, Ref->RPNHead.ig1, param->Lon0, Ref->RPNHead.ig2, param->Lat0,
            Ref->RPNHead.ig3, param->Yaw0, param->NumElem, param->Degree);

    param->LocalToGlobal = make_cs_rotation(param->Lon0, param->Lat0, param->Yaw0);
    param->GlobalToLocal = invert_cs_rotation(param->LocalToGlobal);

    param->NumAxisPoints = param->NumElem * param->Degree;
    param->NumPanelPoints = param->NumAxisPoints * param->NumAxisPoints;

    if (Ref->AX) free(Ref->AX);
    if (Ref->AY) { free(Ref->AY); Ref->AY = NULL; } // Unused
    if (Ref->Lon) free(Ref->Lon);
    if (Ref->Lat) free(Ref->Lat);

    Ref->Lon = (double*) malloc(param->NumPanelPoints * 6 * sizeof(double));
    Ref->Lat = (double*) malloc(param->NumPanelPoints * 6 * sizeof(double));
    Ref->AX  = make_axis(param->NumElem, param->Degree);

    Ref->RPNHead.grtyp[0] = 'X';
    Ref->RPNHead.grtyp[1] = '\0';
    Ref->GRTYP[0] = 'C';
    Ref->GRTYP[1] = '\0';

    for (int i_panel = 0; i_panel < 6; i_panel++) {
        for (int j = 0; j < param->NumAxisPoints; j++) {
            for (int i = 0; i < param->NumAxisPoints; i++) {
                const Coord3D pt_cube = uv_to_cart((CoordUV){.u = Ref->AX[i], .v = Ref->AX[j], .p = i_panel});
                const Coord3D pt_global = apply_rotation(param->LocalToGlobal, pt_cube);
                const Coord2D ll = cart_to_ll(pt_global.x, pt_global.y, pt_global.z);
                const int index = param->NumPanelPoints * i_panel + j * param->NumAxisPoints + i;
                Ref->Lon[index] = ll.x;
                Ref->Lat[index] = ll.y;
            }
        }
    }

    return Ref;
}

//! Check whether the given value is valid (finite, not NaN and not `no_data`)
static inline int is_valid_data_32(const float z, const float no_data) {
    return isfinite(z) && z != no_data;
}

//! Check whether the given value is valid (finite, not NaN and not `no_data`)
static inline int is_valid_data_64(const double z, const double no_data) {
    return isfinite(z) && z != no_data;
}


//! Compute the weights and grid point indices needed to extrapolate to the given positions
int32_t ComputeLinearInterpIndicesC(
    TGeoRef *Ref,               //!< Grid on which we are interpolating
    const double* X,            //!< X position of the points to where we are interpolating
    const double* Y,            //!< Y position of the points to where we are interpolating
    const int32_t NumPoints,    //!< How many positions to interpolate
    float indices[][6]          //!< [out] The weights and indices of interpolation grid points
) {
    const int32_t numPanelPoints = Ref->CGrid->NumPanelPoints;
    const int32_t numAxisPoints = Ref->CGrid->NumAxisPoints;
    int32_t num_interp = 0;
    for (int i = 0; i < NumPoints; i++) {
        const int panel = (int)((Y[i] + 0.5) / numAxisPoints);
        const double local_y = Y[i] - (panel * numAxisPoints);
        const CoordXY local_xy_pt = (CoordXY){.u = X[i], .v = local_y, .p = panel%6};
        const CSPoint ll = lower_left(local_xy_pt);
        const CSPoint lr = (CSPoint){.x = ll.x + 1, .y = ll.y, .p = ll.p};
        const CSPoint ul = (CSPoint){.x = ll.x, .y = ll.y + 1, .p = ll.p};
        const CSPoint ur = (CSPoint){.x = ll.x + 1, .y = ll.y + 1, .p = ll.p};

        indices[i][0] = X[i] - ll.x; // interpolation factor along X
        indices[i][1] = local_y - ll.y; // interpolation factor along Y
        indices[i][2] = (float)xy_to_index(ll, numAxisPoints);
        indices[i][3] = (float)xy_to_index(lr, numAxisPoints);
        indices[i][4] = (float)xy_to_index(ul, numAxisPoints);
        indices[i][5] = (float)xy_to_index(ur, numAxisPoints);

        if (indices[i][2] > numPanelPoints * 6 || indices[i][3] > numPanelPoints * 6 ||
            indices[i][4] > numPanelPoints * 6 || indices[i][5] > numPanelPoints * 6) {
            Lib_Log(APP_LIBGEOREF, APP_ERROR, "%s: bad index (%d, %d, %d, %d) at point %d. Max index %d\n", __func__,
                (int)indices[i][2], (int)indices[i][3], (int)indices[i][4], (int)indices[i][5], i, numPanelPoints * 6);
            Lib_Log(APP_LIBGEOREF, APP_VERBATIM,
                "The 4 points are (%d, %d)[%d], (%d, %d)[%d], (%d, %d)[%d], (%d, %d)[%d]\n",
                ll.x, ll.y, ll.p, lr.x, lr.y, lr.p, ul.x, ul.y, ul.p, ur.x, ur.y, ur.p);
            continue;
        }

        num_interp++;
    }

    return num_interp;
}

void ApplyLinearInterpC(
    const float indices[][6],   //!< Weights and grid point indices for interpolating
    const int32_t NumPoints,    //!< How many points to interpolate
    const double NoData,        //!< Value that indicates absence of data
    const double* z_in,         //!< The input field from which to interpolate
    double* z_out               //!< [out] The interpolated field
) {
    for (int i = 0; i < NumPoints; i++) {
        const float a1 = indices[i][0];
        const float a2 = indices[i][1];
        const double v1 = z_in[(int)indices[i][2]];
        const double v2 = z_in[(int)indices[i][3]];
        const double v3 = z_in[(int)indices[i][4]];
        const double v4 = z_in[(int)indices[i][5]];
        if (is_valid_data_64(v1, NoData) && is_valid_data_64(v2, NoData) && is_valid_data_64(v3, NoData) &&
            is_valid_data_64(v4, NoData)) {
            z_out[i] = bilinear_interp(v1, v2, v4, v3, 1.0 - a1, 1.0 - a2);
        }
        else {
            z_out[i] = NoData;
        }
    }
}

void ApplyLinearInterpC_32(
    const float indices[][6],   //!< Weights and grid point indices for interpolating
    const int32_t NumPoints,    //!< How many points to interpolate
    const float NoData,         //!< Value that indicates absence of data
    const float* z_in,          //!< The input field from which to interpolate
    float* z_out                //!< [out] The interpolated field
) {
    for (int i = 0; i < NumPoints; i++) {
        const float a1 = indices[i][0];
        const float a2 = indices[i][1];
        const float v1 = z_in[(int)indices[i][2]];
        const float v2 = z_in[(int)indices[i][3]];
        const float v3 = z_in[(int)indices[i][4]];
        const float v4 = z_in[(int)indices[i][5]];
        if (is_valid_data_32(v1, NoData) && is_valid_data_32(v2, NoData) && is_valid_data_32(v3, NoData) &&
            is_valid_data_32(v4, NoData)) {
            z_out[i] = bilinear_interp_32(v1, v2, v4, v3, a1, a2);
        }
        else {
            z_out[i] = NoData;
        }
    }
}

int test_cubed_sphere(void) {
    const int num_elem = 20;
    const int num_solpts = 5;
    TGeoRef* ref = GeoRef_New();
    // ref->RPNHead.ig1 = 0x700000;
    // ref->RPNHead.ig2 = 0x900000;
    // ref->RPNHead.ig3 = 0x800400;
    ref->RPNHead.ig1 = 0x910000;
    ref->RPNHead.ig2 = 0x904000;
    ref->RPNHead.ig3 = 0x7f0000;
    // ref->RPNHead.ig1 = 0x800000;
    // ref->RPNHead.ig2 = 0x800000;
    // ref->RPNHead.ig3 = 0x800000;
    ref->RPNHead.ig4 = encode_cs_ig4(num_elem, num_solpts);
    if (GeoRef_DefineC(ref) != ref) {
        App_Log(APP_ERROR, "Error while creating grid!\n");
        return -1;
    }

    const CubedSphereParams* param = ref->CGrid;

    {
        // Compute array of UV coordinates for every grid point
        double* full_u = (double*) malloc(param->NumPanelPoints * 6 * sizeof(double));
        double* full_v = (double*) malloc(param->NumPanelPoints * 6 * sizeof(double));
        for (int i = 0; i < 6; i++) {
            for (int j = 0; j < param->NumAxisPoints; j++) {
                memcpy(full_u + i * param->NumPanelPoints + j * param->NumAxisPoints, ref->AX,
                       param->NumAxisPoints * sizeof(double));
                for (int k = 0; k < param->NumAxisPoints; k++) {
                    full_v[i * param->NumPanelPoints + j * param->NumAxisPoints + k] = ref->AX[j];
                }
            }
        }

        // Inverse operation (just to check if we get back the correct UV)
        double* tmp_u = (double*) malloc(param->NumPanelPoints * 6 * sizeof(double));
        double* tmp_v = (double*) malloc(param->NumPanelPoints * 6 * sizeof(double));
        for (int i = 0; i < param->NumPanelPoints * 6; i++) {
            const Coord3D pt_global = ll_to_cart(ref->Lon[i], ref->Lat[i]);
            const Coord3D pt_cube = apply_rotation(param->GlobalToLocal, pt_global);
            const CoordUV pt_face = cart_to_uv(pt_cube);
            tmp_u[i] = pt_face.u;
            tmp_v[i] = pt_face.v;
        }

        if (check_diffs(full_u, full_v, tmp_u, tmp_v, param->NumPanelPoints * 6, "Going back", 0) > 0) {
            return -1;
        }

        free(tmp_u);
        free(tmp_v);

        int escape = 0;
        size_t num_errors = 0;
        // Testing xy/uv conversion
        for (int i_panel = 0; i_panel < 6 && !escape; i_panel++) {
            const size_t offset_p = i_panel * param->NumPanelPoints;
            for (int j = 0; j < param->NumAxisPoints && !escape; j++) {
                const size_t offset = j * param->NumAxisPoints + offset_p;
                for (int i = 0; i < param->NumAxisPoints && !escape; i++) {
                    const Coord2D xy = uv_to_xy(
                        (CoordUV){.u = full_u[offset + i], .v = full_v[offset + i], .p = i_panel},
                        num_elem, num_solpts);

                    {
                        const double diff = cart_dist(xy.x, xy.y, i, j + param->NumAxisPoints * i_panel);

                        if (diff > 1e-15 * ref->CGrid->NumAxisPoints) {
                            Lib_Log(APP_LIBGEOREF, APP_ERROR,
                                "%s: Difference! Expected (%2g, %2g), got (%7.2g, %7.2g), diff %.2e\n",
                                __func__, (float)i, (float)j + param->NumAxisPoints * i_panel,
                                xy.x, xy.y, diff);
                            num_errors++;
                        }
                    }

                    {
                        const CoordUV uv = xy_to_uv(xy.x, xy.y, ref->AX, param->NumAxisPoints);
                        const double diff = cart_dist(uv.u, uv.v, full_u[offset + i], full_v[offset + i]);
                        if (diff > 1e-15 || uv.p != i_panel) {
                            Lib_Log(APP_LIBGEOREF, APP_ERROR,
                                "%s: (%g, %g) xy->uv error! Expected (%.3f, %.3f), got (%.3f, %.3f) - %d\n",
                                __func__, xy.x, xy.y, full_u[offset + i], full_v[offset + i], uv.u, uv.v, uv.p);
                            num_errors++;
                        }
                    }

                }
            }
        }
        if (num_errors > 0) {
            Lib_Log(APP_LIBGEOREF, APP_ERROR, "%s: %zu errors in XY-UV conversion\n", __func__, num_errors);
            return -1;
        }
    }

    const int num_pts = num_elem * num_solpts;
    const int num_total = num_pts * num_pts * 6;

    double* grid_x = (double*) malloc(num_total * sizeof(double));
    double* grid_y = (double*) malloc(num_total * sizeof(double));

    for (int j = 0; j < num_pts * 6; j++) {
        for (int i = 0; i < num_pts; i++) {
            const size_t index = j * num_pts + i;
            grid_x[index] = i;
            grid_y[index] = j;
        }
    }

    Lib_Log(APP_LIBGEOREF, APP_INFO, "%s: Checking XY2LL\n", __func__);
    {
        double* lon = (double*) malloc(num_total * sizeof(double));
        double* lat = (double*) malloc(num_total * sizeof(double));

        const int num_converted = GeoRef_XY2LL_C(ref, lat, lon, grid_x, grid_y, num_total);
        if (num_converted != num_total) {
            Lib_Log(APP_LIBGEOREF, APP_ERROR, "%s: We have a problem (XY2LL)\n", __func__);
            return -1;
        }

        size_t num_errors = 0;
        for (int i = 0; i < num_total; i++) {
            const double diff = ll_dist(lon[i], lat[i], ref->Lon[i], ref->Lat[i]);
            if (diff > 1e-16) {
                App_Log(APP_ERROR, "AAAAHhhhh not the same lat-lon (%.2e)\n", diff);
                num_errors++;
            }
        }
        
        const int log_lvl = num_errors > 0 ? APP_ERROR : APP_INFO;
        Lib_Log(APP_LIBGEOREF, log_lvl, "%s: Num XY2LL errors = %zu\n", __func__, num_errors);
        if (num_errors > 0) {
            return -1;
        }
    }
    Lib_Log(APP_LIBGEOREF, APP_INFO, "Checking LL2XY\n");
    {
        double* grid_back_x = (double*) malloc(num_total * sizeof(double));
        double* grid_back_y = (double*) malloc(num_total * sizeof(double));
        const int num_converted = GeoRef_LL2XY_C(ref, grid_back_x, grid_back_y, ref->Lat, ref->Lon, num_total);

        if (num_converted != num_total) {
            Lib_Log(APP_LIBGEOREF, APP_ERROR, "Didn't convert the same number of points!\n");
            return -1;
        }

        size_t num_errors = 0;
        for (int i = 0; i < num_total; i++) {
            const double diff = cart_dist(grid_x[i], grid_y[i], grid_back_x[i], grid_back_y[i]);
            if (diff > 1e-15 * ref->CGrid->NumAxisPoints) {
                App_Log(APP_ERROR, "AAAAHhhhh not the same grid coord (%.2e)\n", diff);
                num_errors++;
            }
        }

        const int log_lvl = num_errors > 0 ? APP_ERROR : APP_INFO;
        Lib_Log(APP_LIBGEOREF, log_lvl, "%s: Num LL2XY errors = %zu\n", __func__, num_errors);
        if (num_errors > 0) {
            return -1;
        }
    }


    Lib_Log(APP_LIBGEOREF, APP_INFO, "%s: Checking LL interp\n", __func__);
    {
        // Create a field with values defined at grid points
        double* grid_point_field = (double*) malloc(param->NumPanelPoints * 6 * sizeof(double));
        for (int i = 0; i < param->NumPanelPoints * 6; i++) {
            grid_point_field[i] = angular_dist_0(ref->Lon[i], ref->Lat[i]);
        }

        // Write grid and field to file
        {
            fst_file* f = fst24_open("C.fst", "R/W");
            if (f == NULL) {
                Lib_Log(APP_LIBGEOREF, APP_ERROR, "%s: Unable to open file to store test grid\n", __func__);
                return -1;
            }
            GeoRef_WriteFST(ref, "^>", ref->RPNHead.ig1, ref->RPNHead.ig2, ref->RPNHead.ig3, ref->RPNHead.ig4, f);
            fst_record rec = default_fst_record;
            
            rec.data = grid_point_field;
            rec.ni   = param->NumAxisPoints;
            rec.nj   = param->NumAxisPoints * 6;
            rec.nk   = 1;
            rec.dateo = 0;
            rec.deet  = 0;
            rec.npas  = 0;
            rec.ip1  = ref->RPNHead.ig1;
            rec.ip2  = ref->RPNHead.ig2;
            rec.ip3  = ref->RPNHead.ig3;
            // strncpy(rec.typvar,"X", FST_TYPVAR_LEN);
            strncpy(rec.nomvar, "DIST", FST_NOMVAR_LEN);
            strncpy(rec.grtyp, "C", FST_GTYP_LEN);
            strncpy(rec.etiket, "TEST_FIELD", FST_ETIKET_LEN);
            rec.ig1   = 0;
            rec.ig2   = 0;
            rec.ig3   = 0;
            rec.ig4   = 0;
            rec.data_type = FST_TYPE_REAL;
            rec.data_bits = 64;
            rec.pack_bits = 32;
            if (fst24_write(f, &rec, FST_SKIP) <= 0) {
                Lib_Log(APP_LIBGEOREF, APP_ERROR, "%s: Could not write test field (fst24_write failed)\n",
                    __func__);
                return -1;
            }
            fst24_close(f);
        }

        TApp_Timer t1 = NULL_TIMER;
        TApp_Timer t2 = NULL_TIMER;
        TApp_Timer t3 = NULL_TIMER;

        const size_t num_samples = 1980;
        double* x = (double*) malloc(num_samples * num_samples * sizeof(double));
        double* y = (double*) malloc(num_samples * num_samples * sizeof(double));
        double* lon = (double*) malloc(num_samples * num_samples * sizeof(double));
        double* lat = (double*) malloc(num_samples * num_samples * sizeof(double));
        double* field_interp = (double*) malloc(num_samples * num_samples * sizeof(double));

        // const double corner_x = 0.0;
        // const double corner_y = 0.0;
        // const double corner_x2 = 1.0;
        // const double corner_y2 = 1.0;
        const double corner_x = -0.5;
        const double corner_y = -0.5;
        const double corner_x2 = ref->CGrid->NumAxisPoints - 0.5;
        const double corner_y2 = ref->CGrid->NumAxisPoints * 6 - 0.5;

        for (int j = 0; j < num_samples; j++) {
            for (int i = 0; i < num_samples; i++) {
                const size_t index = j * num_samples + i;
                x[index] = corner_x + i * (corner_x2 - corner_x) / (num_samples - 1);
                y[index] = corner_y + j * (corner_y2 - corner_y) / (num_samples - 1);
                // App_Log(APP_WARNING, "%zu %.3f %.3f\n", index, x[index], y[index]);
            }
        }

        // Compute lat-lon coords for given XY points
        App_TimerStart(&t1);
        if (GeoRef_XY2LL_C(ref, lat, lon, x, y, num_samples * num_samples) != num_samples * num_samples) {
            Lib_Log(APP_LIBGEOREF, APP_ERROR, "%s: Could not convert everything (XY2LL)\n", __func__);
            return -1;
        }
        App_TimerStop(&t1);

        // Interpolate field to these points
        App_TimerStart(&t3);
        float (*indices)[6] = malloc(num_samples * num_samples * 6 * sizeof(float));
        ComputeLinearInterpIndicesC(ref, x, y, num_samples * num_samples, indices);
        ApplyLinearInterpC(indices, num_samples * num_samples, 10e38, grid_point_field, field_interp);
        App_TimerStop(&t3);

        
        // We will skip error checking at the min/max poles, since interpolation error will be high there.
        double pole_x[2], pole_y[2];
        double pole_lon[] = {0.0, M_PI};
        double pole_lat[] = {0.0, 0.0};
        GeoRef_LL2XY_C(ref, pole_x, pole_y, pole_lat, pole_lon, 2);
        {
            const int panel = (int)((pole_y[0] + 0.5) / ref->CGrid->NumAxisPoints);
            double local_y = pole_y[0] - (panel * ref->CGrid->NumAxisPoints);
            if (local_y < 0.0   || local_y > ref->CGrid->NumAxisPoints - 1 ||
                pole_x[0] < 0.0 || pole_x[0] > ref->CGrid->NumAxisPoints) {
                Lib_Log(APP_LIBGEOREF, APP_ERROR,
                    "%s: Min/max poles are at the edge of a panel. Change your grid parameters\n", __func__);
                return -1;
            }
        }
        const double min_corner_x0 = (int)pole_x[0];
        const double min_corner_x1 = min_corner_x0 + 1;
        const double max_corner_x0 = (int)pole_x[1];
        const double max_corner_x1 = max_corner_x0 + 1;

        const double min_corner_y0 = (int)pole_y[0];
        const double min_corner_y1 = min_corner_y0 + 1;
        const double max_corner_y0 = (int)pole_y[1];
        const double max_corner_y1 = max_corner_y0 + 1;

        Lib_Log(APP_LIBGEOREF, APP_WARNING, "%s: Skipping around (%f, %f) and (%f, %f)\n",
            __func__, pole_x[0], pole_y[0], pole_x[1], pole_y[1]);
        Lib_Log(APP_LIBGEOREF, APP_WARNING, "%s: Corners %f %f %f %f %f %f %f %f\n", __func__,
            min_corner_x0, min_corner_x1, min_corner_y0, min_corner_y1,
            max_corner_x0, max_corner_x1, max_corner_y0, max_corner_y1
        );

        double max_error2 = 0.0;
        double total_error2 = 0.0;
        int num_large_errors = 0;
        int num_skipped = 0;
        int num_printed = 0;
        for (int i = 0; i < num_samples * num_samples; i++) {
            // const double dist = ll_dist(lon[i], lat[i], lon_interp[i], lat_interp[i]);
            // total_error1 += dist;
            // if (dist > max_error1) max_error1 = dist;
            if ((x[i] > min_corner_x0 && x[i] < min_corner_x1 && y[i] > min_corner_y0 && y[i] < min_corner_y1) ||
                (x[i] > max_corner_x0 && x[i] < max_corner_x1 && y[i] > max_corner_y0 && y[i] < max_corner_y1))
            {
                // Skip this point, it's near a pole
                num_skipped++;
                continue;
            }

            const double expected = angular_dist_0(lon[i], lat[i]);
            const double dist2 = fabs(field_interp[i] - expected);
            total_error2 += dist2;
            if (dist2 > 2e-3) {
                if (num_printed < 5 && (expected > 0.056 && expected < 3.09)) {
                    Lib_Log(APP_LIBGEOREF, APP_ERROR,
                        "%s: [%6d] Large error %f at (%6.3f, %7.3f) -> (%6.3f, %6.3f) "
                        "Got %.3f, expected %.3f\n",
                        __func__, i, dist2, x[i], y[i], lon[i], lat[i], field_interp[i], expected);
                    const int ind[] = {(int)indices[i][2], (int)indices[i][3], (int)indices[i][4], (int)indices[i][5]};
                    Lib_Log(APP_LIBGEOREF, APP_VERBATIM,
                        "Weights %f %f, indices %d, %d, %d, %d\n",
                        indices[i][0], indices[i][1], ind[0], ind[1], ind[2], ind[3]);
                    Lib_Log(APP_LIBGEOREF, APP_VERBATIM,
                        "Points\n"
                        "  (%6.3f, %6.3f) -> %.3f\n"
                        "  (%6.3f, %6.3f) -> %.3f\n"
                        "  (%6.3f, %6.3f) -> %.3f\n"
                        "  (%6.3f, %6.3f) -> %.3f\n",
                        ref->Lon[ind[0]], ref->Lat[ind[0]], grid_point_field[ind[0]],
                        ref->Lon[ind[1]], ref->Lat[ind[1]], grid_point_field[ind[1]],
                        ref->Lon[ind[2]], ref->Lat[ind[2]], grid_point_field[ind[2]],
                        ref->Lon[ind[3]], ref->Lat[ind[3]], grid_point_field[ind[3]]
                    );
                    num_printed++;
                }
                num_large_errors++;
            }

            if (dist2 > max_error2) max_error2 = dist2;
        }

        const double method_accurate = App_TimerTotalTime_ms(&t1);
        // const double method_interp = App_TimerTotalTime_ms(&t2);
        const double method_func = App_TimerTotalTime_ms(&t3);
        // const double error_avg1 = total_error1 / (num_samples * num_samples);
        const double error_avg2 = total_error2 / (num_samples * num_samples - num_skipped);
        Lib_Log(APP_LIBGEOREF, APP_INFO,
            "Avg error %.2e, max %.2e (%d samples, %d large errors, %.1f%%, %d skipped). "
            "Exact method %.2f ms, function method %.2f\n",
            error_avg2, max_error2, num_samples * num_samples, num_large_errors, num_skipped,
            num_large_errors * 100.0 / (num_samples * num_samples - num_skipped),
            method_accurate, method_func);
        if (error_avg2 > 3e-5 || max_error2 > 5e-3) {
            Lib_Log(APP_LIBGEOREF, APP_ERROR, "%s: Error is too large (avg %.2e, max %.2e)\n",
                __func__, error_avg2, max_error2);
            return -1;
        }
    }

    return 0;
}
