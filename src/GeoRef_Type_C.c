#include <App.h>

#include "GeoRef.h"
#include "georef/GeoRef_Type_C.h"

static const double AXIS_MIN = -M_PI4;
static const double AXIS_MAX = M_PI4;
static const double AXIS_RANGE = M_PI2;

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

typedef struct {
    int32_t u, v;
    PanelID p;
} CSPoint;

typedef struct {
    double u, v;
    PanelID p;
} UVCoord;

typedef struct {
    double x, y, z;
} Coord3D;

typedef struct {
    double x, y;
} Coord2D;

typedef struct {
} RotationSet;

// static RotationSet cs_rotation;
// static double* axis_rad = NULL;
// static double* lats = NULL;
// static double* lons = NULL;
// static int num_axis_points = 0;
// static int num_panel_points = 0;
// static int num_elem = 0;
// static int degree = 0;

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

static inline double decode_cs_angle(const int32_t angle24) {
    const double interval = M_2PI / 0x1000000;
    return ((angle24 & 0xffffff) - 0x800000) * interval;
}

static inline void decode_cs_ig4(const int32_t ig4, int32_t* num_elements, int32_t* num_solpts) {
    *num_elements = ig4 >> 7;
    *num_solpts = ig4 & 0x7f;
}

static inline int32_t encode_cs_ig4(const int32_t num_elements, const int32_t num_solpts) {
    return ((num_elements & 0x1ffff) << 7) | (num_solpts & 0x3f);
}

//! Create a rotation matrix that composes three rotations in cartesian coordinates (1 around each axis), 
//! in this order: Z, X, Y
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

//! Apply a rotation to a point (in cartesian coordinates)
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

//! Retrieve the closest point located to the lower left of the given coordinates
static inline CSPoint lower_left(const UVCoord Coord) {
    return (CSPoint){.u = (int32_t)floor(Coord.u), .v = (int32_t)floor(Coord.v), .p = Coord.p};
}

static inline CSPoint make_uv(const int32_t u, const int32_t v, const PanelID p) {
    return (CSPoint){.u = u, .v = v, .p = p};
}

static inline CSPoint local_uv(const CSPoint Coord, const int32_t PanelSize) {
    const CSPoint error = (CSPoint) {.u = -1, .v = -1, .p = PANEL_NONE};

    if (Coord.u <= -PanelSize || Coord.u > 2*PanelSize || Coord.v <= -PanelSize || Coord.v >= 2*PanelSize) {
        Lib_Log(APP_LIBGEOREF, APP_ERROR, "%s: Point (%d, %d) is too far out! (Panel size %d points)\n",
                __func__, Coord.u, Coord.v, PanelSize);
        return error;
    }

    switch (Coord.p) {
    case PANEL0:
        if (Coord.u < 0)          return make_uv(Coord.u + PanelSize, Coord.v, PANEL3);
        if (Coord.u >= PanelSize) return make_uv(Coord.u - PanelSize, Coord.v, PANEL1);
        if (Coord.v < 0)          return make_uv(Coord.u, Coord.v + PanelSize, PANEL4);
        if (Coord.v >= PanelSize) return make_uv(Coord.u, Coord.v - PanelSize, PANEL5);
        break;
    case PANEL1:
        if (Coord.u < 0)          return make_uv(Coord.u + PanelSize, Coord.v, PANEL0);
        if (Coord.u >= PanelSize) return make_uv(Coord.u - PanelSize, Coord.v, PANEL2);
        if (Coord.v < 0)          return make_uv(Coord.v + PanelSize, PanelSize - Coord.u - 1, PANEL5);
        if (Coord.v >= PanelSize) return make_uv(Coord.v - PanelSize, Coord.u, PANEL4);
        break;
    case PANEL2:
        if (Coord.u < 0)          return make_uv(Coord.u + PanelSize, Coord.v, PANEL1);
        if (Coord.u >= PanelSize) return make_uv(Coord.u - PanelSize, Coord.v, PANEL3);
        if (Coord.v < 0)          return make_uv(PanelSize - 1 - Coord.u, 1 + Coord.v, PANEL5);
        if (Coord.v >= PanelSize) return make_uv(PanelSize - 1 - Coord.u, 2*PanelSize - Coord.v - 1, PANEL4);
        break;
    case PANEL3:
        if (Coord.u < 0)          return make_uv(Coord.u + PanelSize, Coord.v, PANEL2);
        if (Coord.u >= PanelSize) return make_uv(Coord.u - PanelSize, Coord.v, PANEL0);
        if (Coord.v < 0)          return make_uv(-Coord.v - 1, Coord.u, PANEL5);
        if (Coord.v >= PanelSize) return make_uv(Coord.v - PanelSize, PanelSize - 1 - Coord.u, PANEL4);
        break;
    case PANEL4:
        if (Coord.u < 0)          return make_uv(PanelSize - 1 - Coord.v, Coord.u + PanelSize, PANEL3);
        if (Coord.u >= PanelSize) return make_uv(Coord.v, 2*PanelSize - Coord.u - 1, PANEL1);
        if (Coord.v < 0)          return make_uv(Coord.u, Coord.v + PanelSize, PANEL0);
        if (Coord.v >= PanelSize) return make_uv(Coord.u, Coord.v - PanelSize, PANEL5);
        break;
    case PANEL5:
        if (Coord.u < 0)          return make_uv(Coord.v, -Coord.u - 1, PANEL3);
        if (Coord.u >= PanelSize) return make_uv(PanelSize - 1 - Coord.v, Coord.u - PanelSize, PANEL1);
        if (Coord.v < 0)          return make_uv(Coord.u, Coord.v + PanelSize, PANEL4);
        if (Coord.v >= PanelSize) return make_uv(Coord.u, Coord.v - PanelSize, PANEL0);
        break;
    default:
        Lib_Log(APP_LIBGEOREF, APP_ERROR, "%s: Invalid panel! (%d)\n", __func__, Coord.p);
        return error;
    }

    return Coord;
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
static inline double bilinear_interp(const double a, const double b, const double c, const double d, const double alpha_1,
                             const double alpha_2
) {
    const double mid_1 = a * alpha_1 + b * (1.0 - alpha_1);
    const double mid_2 = d * alpha_1 + c * (1.0 - alpha_1);
    const double result = mid_1 * alpha_2 + mid_2 * (1.0 - alpha_2);
    return result;
}

static inline Coord2D cart_to_ll(const double x, const double y, const double z) {
    const double xz = sqrt(x*x + z*z);
    return (Coord2D){
        .x = atan2(x, z),
        .y = atan2(y, xz),
    };
}

//! Assumes from panel 0
static inline Coord3D rotate_to_panel(const Coord3D pt, const PanelID target_panel) {
    // Lib_Log(APP_LIBGEOREF, APP_VERBATIM, "Rotating (%f, %f, %f) to panel %d\n",
    //         pt.x, pt.y, pt.z, target_panel);
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

//! To panel 0
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

static inline Coord3D uv_to_cart(const UVCoord pt) {
    Coord3D tmp = (Coord3D){
        .x = tan(pt.u),
        .y = tan(pt.v),
        .z = 1.0,
    };
    return rotate_to_panel(tmp, pt.p);
}

static inline UVCoord cart_to_uv(const Coord3D pt) {
    const PanelID p = find_panel(pt);
    const Coord3D tmp = rotate_from_panel(pt, p);
    return (UVCoord) {
        .u = atan2(tmp.x, tmp.z),
        .v = atan2(tmp.y, tmp.z),
        .p = p,
    };
}

static inline double x_to_u(const double x, const double* axis_points, const int32_t num_points) {
    if (x < 0.0) {
        return (axis_points[0] * (0.5+x) - AXIS_MIN * x) / 0.5;
    }
    else if (x > (num_points - 1)) {
        const double alpha = x - (num_points - 1);
        return (axis_points[num_points - 1] * (0.5-alpha) + AXIS_MAX * alpha) / 0.5;
    }

    double x_int;
    const double alpha = modf(x, &x_int);
    const int32_t low = (int32_t)x_int;
    const int32_t high = low + 1;
    return axis_points[low] * (1-alpha) + axis_points[high] * alpha;
}

static inline UVCoord xy_to_uv(const double x, const double y, const double* axis_points, const int32_t num_points) {
    const UVCoord error = (UVCoord){.u = 0.0, .v = 0.0, .p = PANEL_NONE};
    if (x < -0.5 || x > num_points - 0.5 || y < -0.5 || y > 6*num_points - 0.5) {
        Lib_Log(APP_LIBGEOREF, APP_ERROR, "%s: (%f, %f) out of range ([%f, %f], [%f, %f])\n",
                __func__, x, y, -0.5, num_points - 0.5, -0.5, 6 * num_points - 0.5);
        return error;
    }

    const int panel = (int)((y + 0.5) / num_points);
    const double local_y = y - (panel * num_points);
    // Lib_Log(APP_LIBGEOREF, APP_VERBATIM, "%s: panel %d, y = %.3f (local %.3f)\n", __func__, panel, y, local_y);

    UVCoord c;
    c.p = panel;
    c.u = x_to_u(x, axis_points, num_points);
    c.v = x_to_u(local_y, axis_points, num_points);
    return c;
}

static inline double refu_to_refx(const double u, const int32_t degree) {
    const double* points = QUAD_POINTS[degree - 1];
    if (u <= points[0]) {
        // Lib_Log(APP_LIBGEOREF, APP_VERBATIM, "%s: %.3f is between edge and point 0\n", __func__, u);
        return 0.5 * (u - points[0]) / (points[0] + 1.0);
    }
    for (int i = 1; i < degree; i++) {
        if (u <= points[i]) {
            // Lib_Log(APP_LIBGEOREF, APP_VERBATIM, "%s: %.3f is between points %d and %d\n", __func__, u,
            //     i-1, i);
            return (u - points[i]) / (points[i] - points[i-1]) + i;
        }
    }

    // Lib_Log(APP_LIBGEOREF, APP_VERBATIM, "%s: %.3f is between point %d and edge\n", __func__, u, degree - 1);
    return degree - 1 + 0.5 * (u - points[degree - 1]) / (1.0 - points[degree-1]);
}

static inline Coord2D uv_to_xy(
    const UVCoord pt,
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

//! Points are ordered counter-clockwise
static inline double interp4(const double* values, const size_t indices[4], const double weights[2]) {
    Lib_Log(APP_LIBGEOREF, APP_WARNING, "%s: indices = %6lu %6lu %6lu %6lu, weights = %f, %f\n",
            __func__, indices[0], indices[1], indices[2], indices[3], weights[0], weights[1]);
    const double v1 = values[indices[0]] * weights[0] + values[indices[1]] * (1.0 - weights[0]);
    const double v2 = values[indices[3]] * weights[0] + values[indices[2]] * (1.0 - weights[0]);
    return v1 * weights[1] + v2 * (1.0 - weights[1]);
}

//! Retrieve the index of point (X, Y) in the array of all points
static inline size_t uv_to_index(const CSPoint Coord, const int32_t PanelSize) {
    const CSPoint local = local_uv(Coord, PanelSize);
    if (local.p == PANEL_NONE) return (size_t)-1;
    return (local.v + PanelSize * (int32_t)local.p) * PanelSize + local.u;
}

static inline size_t point_index(const int32_t x, const int32_t y, const int32_t panel, const int32_t num_axis_points) {
    const size_t panel_offset = num_axis_points * num_axis_points * panel;
    const size_t row_offset = y * num_axis_points;
    return panel_offset + row_offset + x;
}

/*----------------------------------------------------------------------------
 * @brief  Transforms XY grid coordinates to LatLon for a C (cube-sphere) grid
 *    @param[in]  Ref     Georeference pointer
 *    @param[out] Lat     Latitude array
 *    @param[out] Lon     Longitude array
 *    @param[in]  X       X array
 *    @param[in]  Y       Y array
 *    @param[in]  Nb      Number of coordinates

 *    @return             Error code (0=ok)
*/
int32_t GeoRef_XY2LL_C(TGeoRef *Ref, double *Lat, double *Lon, double *X, double *Y, int32_t Nb) {

    int32_t num_points = 0;
    const int32_t num_axis_points = Ref->CGrid->NumAxisPoints;
    const RotationParam rot = Ref->CGrid->LocalToGlobal;
    for (int i = 0; i < Nb; i++) {
        const UVCoord pt_face = xy_to_uv(X[i], Y[i], Ref->AX, num_axis_points);
        if (pt_face.p == PANEL_NONE) continue;
        const Coord3D pt_cube = uv_to_cart(pt_face);
        const Coord3D pt_global = apply_rotation(rot, pt_cube);
        const Coord2D ll = cart_to_ll(pt_global.x, pt_global.y, pt_global.z);
        Lon[i] = ll.x;
        Lat[i] = ll.y;
        num_points++;
    }

    // {
    //     char buffer[2048 * 10];
    //     char* ptr = buffer;
    //     for (int i = 0; i < Nb; i++) {
    //         if (i % 2 == 0) ptr += sprintf(ptr, "\n");
    //         ptr += sprintf(ptr, "(%8.3f, %8.3f) -> (%8.3f, %8.3f)  |  ", X[i], Y[i], Lon[i], Lat[i]);
    //     }
    //     Lib_Log(APP_LIBGEOREF, APP_ALWAYS, "%s: Values = %s\n", __func__, buffer);
    // }

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

 *    @return             Error code (0=ok)
*/
int32_t GeoRef_LL2XY_C(TGeoRef *Ref, double *X, double *Y, double *Lat, double *Lon, int32_t Nb) {
    int32_t num_points = 0;

    const RotationParam rot = Ref->CGrid->GlobalToLocal;
    for (int i = 0; i < Nb; i++) {
        const Coord3D pt_global = ll_to_cart(Lon[i], Lat[i]);
        const Coord3D pt_cube = apply_rotation(rot, pt_global);
        const UVCoord pt_face = cart_to_uv(pt_cube);
        const Coord2D xy = uv_to_xy(pt_face, Ref->CGrid->NumElem, Ref->CGrid->Degree);
        X[i] = xy.x;
        Y[i] = xy.y;
        num_points++;
    }

    // {
    //     char buffer[2048 * 10];
    //     char* ptr = buffer;
    //     for (int i = 0; i < Nb && i < 16; i++) {
    //         if (i % 4 == 0) ptr += sprintf(ptr, "\n");
    //         ptr += sprintf(ptr, "(%8.3f, %8.3f) ", Lat[i], Lon[i]);
    //     }
    //     Lib_Log(APP_LIBGEOREF, APP_ALWAYS, "%s: Values = %s\n", __func__, buffer);
    // }
    
    return num_points;
}


static inline double r2d(const double r) { return r * 180 / M_PI; }
static inline double d2r(const double d) { return d / 180 * M_PI; }

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

    Lib_Log(APP_LIBGEOREF, APP_ERROR, "%s: [%s] Num errors %3d (%4.1f%%) avg dist %10.2e\n",
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

    for (int i_panel = 0; i_panel < 6; i_panel++) {
        for (int j = 0; j < param->NumAxisPoints; j++) {
            for (int i = 0; i < param->NumAxisPoints; i++) {
                const Coord3D pt_cube = uv_to_cart((UVCoord){.u = Ref->AX[i], .v = Ref->AX[j], .p = i_panel});
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

int test_cubed_sphere(void) {
    const int num_elem = 2;
    const int num_solpts = 3;
    TGeoRef* ref = GeoRef_New();
    // ref->RPNHead.ig1 = 0x700000;
    // ref->RPNHead.ig2 = 0x900000;
    // ref->RPNHead.ig3 = 0x800400;
    ref->RPNHead.ig1 = 0x810000;
    ref->RPNHead.ig2 = 0x804000;
    ref->RPNHead.ig3 = 0x8f0000;
    ref->RPNHead.ig4 = encode_cs_ig4(num_elem, num_solpts);
    if (GeoRef_DefineC(ref) != ref) {
        App_Log(APP_ERROR, "Error while creating grid!\n");
        return -1;
    }

    const CubedSphereParams* param = ref->CGrid;

    {
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

        // Inverse operation (just to check)
        double* tmp_u = (double*) malloc(param->NumPanelPoints * 6 * sizeof(double));
        double* tmp_v = (double*) malloc(param->NumPanelPoints * 6 * sizeof(double));
        for (int i = 0; i < param->NumPanelPoints * 6; i++) {
            const Coord3D pt_global = ll_to_cart(ref->Lon[i], ref->Lat[i]);
            const Coord3D pt_cube = apply_rotation(param->GlobalToLocal, pt_global);
            const UVCoord pt_face = cart_to_uv(pt_cube);
            tmp_u[i] = pt_face.u;
            tmp_v[i] = pt_face.v;
        }

        {
            char buffer[1024 * 10];
            char* ptr = buffer;
            const int CHUNK = 6;
            for (int i = 0; i < param->NumPanelPoints * 6; i += CHUNK) {
                for (int j = i; j < i + CHUNK; j++) {
                    ptr += sprintf(ptr, "%6.3f ", tmp_u[j]);
                }
                ptr += sprintf(ptr, "  ");
                for (int j = i; j < i + CHUNK; j++) {
                    ptr += sprintf(ptr, "%6.3f ", full_u[j]);
                }
                ptr += sprintf(ptr, "\n");
            }
            Lib_Log(APP_LIBGEOREF, APP_WARNING, "%s: U: (tmp vs full)\n%s\n", __func__, buffer);

            ptr = buffer;
            for (int i = 0; i < param->NumPanelPoints * 6; i += CHUNK) {
                for (int j = i; j < i + CHUNK; j++) {
                    ptr += sprintf(ptr, "%6.3f ", tmp_v[j]);
                }
                ptr += sprintf(ptr, "  ");
                for (int j = i; j < i + CHUNK; j++) {
                    ptr += sprintf(ptr, "%6.3f ", full_v[j]);
                }
                ptr += sprintf(ptr, "\n");
            }
            Lib_Log(APP_LIBGEOREF, APP_WARNING, "%s: V: (tmp vs full)\n%s\n", __func__, buffer);
        }

        if (check_diffs(full_u, full_v, tmp_u, tmp_v, param->NumPanelPoints * 6, "Going back", 0) > 0) {
            return -1;
        }

        free(tmp_u);
        free(tmp_v);

        int escape = 0;
        // Testing xy/uv conversion
        for (int i_panel = 0; i_panel < 6 && !escape; i_panel++) {
            const size_t offset_p = i_panel * param->NumPanelPoints;
            for (int j = 0; j < param->NumAxisPoints && !escape; j++) {
                const size_t offset = j * param->NumAxisPoints + offset_p;
                for (int i = 0; i < param->NumAxisPoints && !escape; i++) {
                    const Coord2D xy = uv_to_xy(
                        (UVCoord){.u = full_u[offset + i], .v = full_v[offset + i], .p = i_panel},
                        num_elem, num_solpts);

                    // Lib_Log(APP_LIBGEOREF, APP_VERBATIM, "(%g, %g)\n", xy.x, xy.y);
                    {
                        const double diff_x = rel_diff(xy.x, i);
                        const double diff_y = rel_diff(xy.y, (j + param->NumAxisPoints * i_panel));

                        if (diff_x > 1e-14 || diff_y > 1e-14) {
                            Lib_Log(APP_LIBGEOREF, APP_ERROR,
                                "%s: Difference! Expected (%2g, %2g), got (%7.2g, %7.2g), diff (%.2e, %.2e)\n",
                                __func__, (float)i, (float)j + param->NumAxisPoints * i_panel,
                                xy.x, xy.y, diff_x, diff_y);
                            escape = 1;
                        }
                    }

                    {
                        const UVCoord uv = xy_to_uv(xy.x, xy.y, ref->AX, param->NumAxisPoints);
                        const double diff_u = rel_diff(uv.u, full_u[offset + i]);
                        const double diff_v = rel_diff(uv.v, full_v[offset + i]);
                        if (diff_u > 1e-14 || diff_v > 1e-14 || uv.p != i_panel) {
                            Lib_Log(APP_LIBGEOREF, APP_ERROR,
                                "%s: (%g, %g) xy->uv error! Expected (%.3f, %.3f), got (%.3f, %.3f) - %d\n",
                                __func__, xy.x, xy.y, full_u[offset + i], full_v[offset + i], uv.u, uv.v, uv.p);
                                escape = 1;
                        }
                    }

                }
            }
        }
        if (escape) return -1;

        escape = 0;
        for (double y = -0.5; (y < 6 * param->NumAxisPoints - 0.5) && !escape; y += 0.12345) {
            for (double x = -0.5; (x < param->NumAxisPoints - 0.5) && !escape; x += 0.12456) {
                // Compute expected uv values
                // We make use of the fact that u is constant along the vertical axis and v is constant along the
                // horizontal axis. Otherwise we would need to do a bilinear interpolation rather than 2 linear ones.
                const int32_t x_im = floor(x);
                const int32_t x_ip = x_im + 1;

                const int panel = (int)((y + 0.5) / param->NumAxisPoints);
                const double y_local = y - panel * param->NumAxisPoints;
                const int32_t y_im = floor(y_local);
                const int32_t y_ip = y_im + 1;

                const double alpha_x = x_im >= 0 ? x_ip < param->NumAxisPoints ? 
                                x - x_im :
                                (x - param->NumAxisPoints + 1) * 2.0:
                                (x + 0.5) * 2.0;
                const double alpha_y = y_im >= 0 ? y_ip < param->NumAxisPoints ?
                                y_local - y_im :
                                (y_local - param->NumAxisPoints + 1) * 2.0 :
                                (y_local + 0.5) * 2.0;

                const size_t i_um = point_index(x_im, y_im >= 0 ? y_im : y_ip, panel, param->NumAxisPoints);
                const size_t i_up = point_index(x_ip, y_im >= 0 ? y_im : y_ip, panel, param->NumAxisPoints);
                const size_t i_vm = point_index(x_im >= 0 ? x_im : x_ip, y_im, panel, param->NumAxisPoints);
                const size_t i_vp = point_index(x_im >= 0 ? x_im : x_ip, y_ip, panel, param->NumAxisPoints);

                const double u_m = x_im < 0 ? AXIS_MIN : full_u[i_um];
                const double u_p = x_ip >= param->NumAxisPoints ? AXIS_MAX : full_u[i_up];
                const double v_m = y_im < 0 ? AXIS_MIN : full_v[i_vm];
                const double v_p = y_ip >= param->NumAxisPoints ? AXIS_MAX : full_v[i_vp];

                const double target_u = u_m * (1.0 - alpha_x) + u_p * alpha_x;
                const double target_v = v_m * (1.0 - alpha_y) + v_p * alpha_y;

                // Lib_Log(APP_LIBGEOREF, APP_VERBATIM,
                //     "%s: (%.3g, %.3g) -> ([%2d, %2d], [%2d, %2d]), uv ([%5.2f, %5.2f], [%5.2f, %5.2f]) "
                //     "indices (%3llu, %3llu, %3llu, %3llu), alphas (%.2f, %.2f)\n",
                //     __func__, x, y, x_im, x_ip, y_im, y_ip, u_m, u_p, v_m, v_p, i_um, i_up, i_vm, i_vp, alpha_x, alpha_y);

                // Do the conversion (back and forth)
                const UVCoord uv = xy_to_uv(x, y, ref->AX, param->NumAxisPoints);
                const Coord2D xy = uv_to_xy(uv, num_elem, num_solpts);

                const double diff_x = fabs(xy.x - x);
                const double diff_y = fabs(xy.y - y);
                if (diff_x > 2e-14 || diff_y > 2e-14) {
                    Lib_Log(APP_LIBGEOREF, APP_ERROR,
                        "%s: Difference! Expected (%7.2g, %7.2g), got (%7.2g, %7.2g), diff (%.2e, %.2e)\n",
                        __func__, x, y, xy.x, xy.y, diff_x, diff_y);
                    escape = 1;
                }

                const double diff_u = fabs(uv.u - target_u);
                const double diff_v = fabs(uv.v - target_v);
                if (diff_u > 1e-14 || diff_v > 1e-14 || uv.p != panel) {
                    Lib_Log(APP_LIBGEOREF, APP_ERROR,
                        "%s: (%g, %g) xy->uv error! Expected (%.3f, %.3f), got (%.3f, %.3f) - %d\n",
                        __func__, x, y, target_u, target_v, uv.u, uv.v, uv.p);
                    escape = 1;
                }
            }
        }

        if (escape) return -1;
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

    App_Log(APP_INFO, "Checking XY2LL\n");
    {
        double* lon = (double*) malloc(num_total * sizeof(double));
        double* lat = (double*) malloc(num_total * sizeof(double));

        const int num_converted = GeoRef_XY2LL_C(ref, lat, lon, grid_x, grid_y, num_total);
        if (num_converted != num_total) {
            App_Log(APP_ERROR, "We have a problem\n");
            return -1;
        }

        size_t num_errors = 0;
        for (int i = 0; i < num_total; i++) {
            const double diff1 = fabs(lon[i] - ref->Lon[i]);
            const double diff2 = fabs(lat[i] - ref->Lat[i]);
            if (diff1 > 1e-16 || diff2 > 1e-16) {
                App_Log(APP_ERROR, "AAAAHhhhh not the same (%e, %e)\n", diff1, diff2);
                num_errors++;
            }
        }
        App_Log(APP_INFO, "Num errors = %zu\n", num_errors);
        if (num_errors > 0) {
            App_Log(APP_ERROR, "Num errors = %zu\n", num_errors);
            return -1;
        }
    }
    App_Log(APP_INFO, "Checking LL2XY\n");
    {
        double* grid_back_x = (double*) malloc(num_total * sizeof(double));
        double* grid_back_y = (double*) malloc(num_total * sizeof(double));
        const int num_converted = GeoRef_LL2XY_C(ref, grid_back_x, grid_back_y, ref->Lat, ref->Lon, num_total);

        if (num_converted != num_total) {
            App_Log(APP_ERROR, "Didn't convert the same number of points!\n");
            return -1;
        }

        size_t num_errors = 0;
        for (int i = 0; i < num_total; i++) {
            const double diff1 = fabs(grid_back_x[i] - grid_x[i]);
            const double diff2 = fabs(grid_back_y[i] - grid_y[i]);
            if (diff1 > 3e-14 || diff2 > 3e-14) {
                App_Log(APP_ERROR, "AAAAHhhhh not the same grid coord (%e, %e)\n", diff1, diff2);
                num_errors++;
            }
        }

        App_Log(APP_INFO, "Num errors = %zu\n", num_errors);
        if (num_errors > 0) {
            App_Log(APP_ERROR, "Num errors = %zu\n", num_errors);
            return -1;
        }
    }


    App_Log(APP_INFO, "Checking LL interp\n");
    {
        TApp_Timer t1 = NULL_TIMER;
        TApp_Timer t2 = NULL_TIMER;
        // const int base0 = 0;
        // const int base1 = num_pts;
        size_t num_errors = 0;
        double max_error = 0.0;

        const size_t num_samples = 3000;
        double* x = (double*) malloc(num_samples * num_samples * sizeof(double));
        double* y = (double*) malloc(num_samples * num_samples * sizeof(double));
        double* lon = (double*) malloc(num_samples * num_samples * sizeof(double));
        double* lat = (double*) malloc(num_samples * num_samples * sizeof(double));
        double* lon_interp = (double*) malloc(num_samples * num_samples * sizeof(double));
        double* lat_interp = (double*) malloc(num_samples * num_samples * sizeof(double));

        const double corner_x = 0.0;
        const double corner_y = 0.0;

        for (int j = 0; j < num_samples; j++) {
            for (int i = 0; i < num_samples; i++) {
                const size_t index = j * num_samples + i;
                x[index] = corner_x + i * 1.0 / num_samples;
                y[index] = corner_y + j * 1.0 / num_samples;
                // App_Log(APP_WARNING, "%zu %.3f %.3f\n", index, x[index], y[index]);
            }
        }

        App_TimerStart(&t1);
        if (GeoRef_XY2LL_C(ref, lat, lon, x, y, num_samples * num_samples) != num_samples * num_samples) {
            App_Log(APP_ERROR, "Could not convert everythin\n");
            return -1;
        }
        App_TimerStop(&t1);

        App_TimerStart(&t2);
        const int num_panel_points = num_pts * num_pts;
        for (int i = 0; i < num_samples * num_samples; i++) {
            const int panel = (int)((y[i] + 0.5) / num_panel_points);
            const double local_y = y[i] - (panel * num_panel_points);

            // const CSPoint ll = lower_left((UVCoord));
            // App_Log(APP_VERBATIM, "x = %.3f, local_y = %.3f, lon/lat = (%.3f, %3f)\n",
            //         x[i], local_y,
            //         ref->AX[(int)x[i] + (int)local_y * num_pts],
            //         ref->AY[(int)x[i] + (int)local_y * num_pts]
            //     );

            const int base0 = (int)local_y * num_pts + (int)x[i];
            const int base1 = base0 + num_pts;
            const double alpha_x = 1.0 - (x[i] - (int)x[i]);
            const double alpha_y = 1.0 - (local_y - (int)local_y);
            
            const double x0 = ref->Lon[base0];
            const double x1 = ref->Lon[base0 + 1];
            const double x2 = ref->Lon[base1 + 1];
            const double x3 = ref->Lon[base1];
            const double y0 = ref->Lat[base0];
            const double y1 = ref->Lat[base0 + 1];
            const double y2 = ref->Lat[base1 + 1];
            const double y3 = ref->Lat[base1];

            lon_interp[i] = bilinear_interp(x0, x1, x2, x3, alpha_x, alpha_y);
            lat_interp[i] = bilinear_interp(y0, y1, y2, y3, alpha_x, alpha_y);

            // App_Log(APP_WARNING, "%d %.3f %.3f\n", i, lon_interp[i], lat_interp[i]);
        }
        App_TimerStop(&t2);

        double total_error = 0.0;
        // double total_size_sq = 0.0;

        double area = 0.0;
        {
            const int base0 = (int)corner_y * num_pts + (int)corner_x;
            const int base1 = base0 + num_pts;
            const double x1 = ref->Lon[base0];
            const double x2 = ref->Lon[base0 + 1];
            const double x3 = ref->Lon[base1 + 1];
            const double x4 = ref->Lon[base1];
            const double y1 = ref->Lat[base0];
            const double y2 = ref->Lat[base0 + 1];
            const double y3 = ref->Lat[base1 + 1];
            const double y4 = ref->Lat[base1];

            area = fabs(x1*y2 - x2*y1 + x2*y3 - x3*y2 + x3*y4 - x4*y3 + x4*y1 - x1*y4) / 2.0;
            App_Log(APP_VERBATIM, "area = %.4f, [(%.3f, %.3f), (%.3f, %.3f), (%.3f, %.3f), (%.3f, %.3f)]\n",
                    area, x1, y1, x2, y2, x3, y3, x4, y4);
        }

        for (int i = 0; i < num_samples * num_samples; i++) {
            // const double diff1 = fabs((lon_interp[i] - lon[i]) / lon[i]);
            // const double diff2 = fabs((lat_interp[i] - lat[i]) / lat[i]);
            // const double dist1 = (lon_interp[i] - lon[i]);
            // const double dist2 = (lat_interp[i] - lat[i]);
            const double dist = ll_dist(lon[i], lat[i], lon_interp[i], lat_interp[i]);
            total_error += dist;
            // total_size_sq += lon[i]*lon[i] + lat[i]*lat[i];

            if (dist > max_error) max_error = dist;
            const double threshold = 1e-2;
            if (dist > threshold) {
                // App_Log(APP_WARNING, "%d (%.3f %.3f) (%.3f %.3f)\n", i, lon[i], lat[i], lon_interp[i], lat_interp[i]);
                // App_Log(APP_ERROR, "Interpolation is different at (%.3f, %.3f) (error %.2e, %.2e)\n",
                //         x[i], y[i], dist/area, dist);
                num_errors++;
            }
        }

        const double method_accurate = App_TimerTotalTime_ms(&t1);
        const double method_interp = App_TimerTotalTime_ms(&t2);
        const double error_norm = total_error / (num_samples * num_samples);
        App_Log(APP_INFO, "%zu errors (%.1f%%), rel %.2e, max %.2e. Exact method %.2f ms, interpolation %.2f ms, %.1f slowdown\n",
                num_errors, num_errors * 100.0 / (num_samples * num_samples), error_norm, max_error,
                method_accurate, method_interp, method_accurate / method_interp);
    }

    return 0;
}
