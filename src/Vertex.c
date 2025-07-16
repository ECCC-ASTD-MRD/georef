#include "georef/Vertex.h"


//! \file


//! Arbitrary quadrilateral location interpolation
void Vertex_Map(
    //! [in] Cell coordinates
    const Vect2d P[4],
    //! [out] Output X grid coordinate
    double * const LX,
    //! [out] Output Y grid coordinate
    double * const LY,
    //! [in] Longitude
    const double WX,
    //! [in] Latitude
    const double WY
) {
    //! The unit cell L(x, y) is oriented as:
    //! L0(x=0, y=0), L1(0, 1), L2(1, 1), L3(1, 0).  The order matters.
    //! Reference: http://math.stackexchange.com/questions/13404/mapping-irregular-quadrilateral-to-a-rectangle

    if (P[0][1] == P[1][1]) P[0][1] -= 1e-5;
    if (P[0][0] == P[3][0]) P[0][0] -= 1e-5;

    const double wydy = (WY      - P[0][1]) / (P[0][1] - P[1][1]);
    const double y3dy = (P[3][1] - P[0][1]) / (P[0][1] - P[1][1]);
    const double y2dy = (P[2][1] - P[0][1]) / (P[0][1] - P[1][1]);
    const double wxdx = (WX      - P[0][0]) / (P[0][0] - P[3][0]);
    const double x1dx = (P[1][0] - P[0][0]) / (P[0][0] - P[3][0]);
    const double x2dx = (P[2][0] - P[0][0]) / (P[0][0] - P[3][0]);
    const double wxy = WX * WY;
    const double xy0 = P[0][0] * P[0][1];
    const double xy1 = P[1][0] * P[1][1];
    const double xy2 = P[2][0] * P[2][1];
    const double xy3 = P[3][0] * P[3][1];

          double ax  = (WX      - P[0][0]) + (P[1][0] - P[0][0]) * wydy;
    const double a3x = (P[3][0] - P[0][0]) + (P[1][0] - P[0][0]) * y3dy;
          double a2x = (P[2][0] - P[0][0]) + (P[1][0] - P[0][0]) * y2dy;
          double ay  = (WY      - P[0][1]) + (P[3][1] - P[0][1]) * wxdx;
    const double a1y = (P[1][1] - P[0][1]) + (P[3][1] - P[0][1]) * x1dx;
          double a2y = (P[2][1] - P[0][1]) + (P[3][1] - P[0][1]) * x2dx;
    const double bx  = wxy - xy0 + (xy1 - xy0) * wydy;
    const double b3x = xy3 - xy0 + (xy1 - xy0) * y3dy;
    const double b2x = xy2 - xy0 + (xy1 - xy0) * y2dy;
    const double by  = wxy - xy0 + (xy3 - xy0) * wxdx;
    const double b1y = xy1 - xy0 + (xy3 - xy0) * x1dx;
    const double b2y = xy2 - xy0 + (xy3 - xy0) * x2dx;

    ax /= a3x;
    ay /= a1y;
    a2x /= a3x;
    a2y /= a1y;

    *LX = ax + (1 - a2x) * (bx - b3x * ax) / (b2x - b3x * a2x);
    *LY = ay + (1 - a2y) * (by - b1y * ay) / (b2y - b1y * a2y);
}

//! Obtenir le gradient d'un point32_t tridimentionne a l'interieur d'un voxel
void VertexGradient(
    //! [] Data definition
    TDef * const Def,
    //! [inout] Normalized gradient / Point a l'interieur du voxel
    Vect3d Nr
) {
    Vect3d v;

    Vect_Assign(v, Nr);
    Nr[0] = VertexVal(Def, -1, v[0]-0.5, v[1], v[2]) - VertexVal(Def, -1, v[0]+0.5, v[1], v[2]);
    Nr[1] = VertexVal(Def, -1, v[0], v[1]-0.5, v[2]) - VertexVal(Def, -1, v[0], v[1]+0.5, v[2]);
    Nr[2] = VertexVal(Def, -1, v[0], v[1], v[2]-0.5) - VertexVal(Def, -1, v[0], v[1], v[2]+0.5);

    // Vect_Mul(Nr, Nr, v);
}


//! Linear interpolation of the cutting position of an isosurface on one side between two vertices
 void VertexInterp(
    //! [out] Resulting point
    Vect3d Pi,
    //! [in] First point
    const Vect3d P0,
    //! [in] Second point
    const Vect3d P1,
    //! [in] Value at point P0
    const double V0,
    //! [in] Value at point P1
    const double V1,
    //! [in] Level
    const double Level
) {
    if (ABS(Level - V0) < TINY_VALUE) {
        Vect_Assign(Pi, P0);
        return;
    }

    if (ABS(Level - V1) < TINY_VALUE) {
        Vect_Assign(Pi, P1);
        return;
    }
    if (ABS(V0 - V1) < TINY_VALUE) {
        Vect_Assign(Pi, P0);
        return;
    }

    const double mu = (Level - V0) / (V1 - V0);

    Pi[0] = P0[0] + mu * (P1[0] - P0[0]);
    Pi[1] = P0[1] + mu * (P1[1] - P0[1]);
    Pi[2] = P0[2] + mu * (P1[2] - P0[2]);
}


//! Interpolate the position of a grid point in X, Y, Z
int32_t VertexLoc(
    Vect3d **Pos,
    //! [in] Data definition
    const TDef * const Def,
    //! [out] Resulting vertex
    Vect3d Vr,
    //! [in] X coordinate ([0, NI - 1])
    double X,
    //! [in] Y coordinate ([0, NJ - 1])
    double Y,
    //! [in] Z coordinate ([0, NK - 1])
    double Z,
    //! [in] Set to true if the grid wraps in X (-180, 180 or 0.360)
    const int32_t Wrap
) {
    if ((X > Def->NI - 1 && !Wrap) || Y > Def->NJ - 1 || Z > Def->NK - 1 || X < 0 || Y < 0 || Z < 0) {
        return 0;
    }

    const unsigned long i = X; X -= i;
    const unsigned long j = Y; Y -= j;
    const unsigned long k = Z; Z -= k;

    // Get gridpoint32_t indexes
    unsigned long idx0 = Def->Idx + j * Def->NI + i % Def->NI;
    unsigned long idx1 = idx0 + 1;
    unsigned long idx3 = (j == Def->NJ - 1) ? idx0 : idx0 + Def->NI;
    unsigned long idx2 = idx3 + 1;

    if (i >= Def->NI - 1) {
        idx1 -= Def->NI;
        idx2 -= Def->NI;
    }

    Vect3d v00, v01, v10, v11, v0, v1;

    // 3D Interpolation case
    if (Z > TINY_VALUE) {
        unsigned long k1 = k + 1;

        Vect_InterpC(v00, Pos[k][idx0], Pos[k1][idx0], Z);
        Vect_InterpC(v10, Pos[k][idx1], Pos[k1][idx1], Z);
        Vect_InterpC(v11, Pos[k][idx2], Pos[k1][idx2], Z);
        Vect_InterpC(v01, Pos[k][idx3], Pos[k1][idx3], Z);
    } else {
        Vect_Assign(v00, Pos[k][idx0]);
        Vect_Assign(v10, Pos[k][idx1]);
        Vect_Assign(v01, Pos[k][idx3]);
        Vect_Assign(v11, Pos[k][idx2]);
    }

    // Interpolate over X
    if (X > TINY_VALUE) {
        Vect_InterpC(v0, v00, v10, X);
        Vect_InterpC(v1, v01, v11, X);
    }  else {
        Vect_Assign(v0, v00);
        Vect_Assign(v1, v01);
    }

    // Interpolate over Y
    if (Y > TINY_VALUE) {
        Vect_InterpC(Vr, v0, v1, Y);
    } else {
        Vect_Assign(Vr, v0);
    }

    return 1;
}

//! Moyenne la valeur d'un point32_t de grille masque avec ses voisins
static inline double VertexAvg(
    //! [in] Data definition
    const TDef *const Def,
    //! [in] Components
    int32_t Idx,
    //! [in] X position
    int32_t X,
    //! [in] Y position
    int32_t Y,
    //! [in] Z position. If less than 0, 2D interpolation only
    int32_t Z
) {
    unsigned long n = 0;
    double val = 0.0;
    const unsigned long k = Z * Def->NIJ;
    double v;
    for(unsigned long j = Y - 1; j <= Y + 1; j++) {
        for(unsigned long i = X - 1; i <= X + 1; i++) {
            unsigned long idx = j * Def->NI + i;
            if (idx >= 0 && idx < Def->NIJ && Def->Mask[k + idx]) {
                if (Idx == -1) {
                    Def_GetMod(Def, k + idx, v);
                } else {
                    Def_Get(Def, Idx, k + idx, v)
                }
                if (DEFVALID(Def, v)) {
                    val += v;
                    n++;
                }
            }
        }
    }
    return n ? val / n : Def->NoData;
}


//! Linear interpolation of the value of a point inside a voxel for one of the two components (UV = speed)
float VertexVal(
    //! [in] Data definition
    const TDef * const Def,
    //! [in] Components
    int32_t Idx,
    //! [in] X position
    double X,
    //! [in] Y position
    double Y,
    //! [in] Z position. If less than 0, 2D interpolation only
    double Z
) {
    if (X > Def->NI - 1 || Y > Def->NJ - 1 || Z > Def->NK - 1 || X < 0 || Y < 0 || Z < 0) return 0;

    unsigned long i = X; X -= i;
    unsigned long j = Y; Y -= j;
    unsigned long k = Z; Z -= k;

    // Get gridpoint32_t indexes
    unsigned long idxk = Def->NIJ;

    unsigned long idx[4];
    idx[0] = (k ? idxk * k : 0) + j * Def->NI + i;
    idx[1] = idx[0] + 1;
    idx[3] = (j == Def->NJ - 1) ? idx[0] : idx[0] + Def->NI;
    idx[2] = idx[3] + 1;

    double cube[2][4];
    if (Idx == -1) {
        Def_GetQuadMod(Def, idx, cube[0]);
    } else {
        Def_GetQuad(Def, Idx, idx, cube[0]);
    }

    // If a mask exists average the masked points with the ones around
    if (Def->Mask) {
        if (!Def->Mask[idx[0]]) cube[0][0] = VertexAvg(Def, Idx, i,  j,  k);
        if (!Def->Mask[idx[1]]) cube[0][1] = VertexAvg(Def, Idx, i+1, j,  k);
        if (!Def->Mask[idx[2]]) cube[0][2] = VertexAvg(Def, Idx, i+1, j+1, k);
        if (!Def->Mask[idx[3]]) cube[0][3] = VertexAvg(Def, Idx, i,  j+1, k);
    }

    // If either value is nodata then interpolation will be nodata as well
    if (!DEFVALID(Def, cube[0][0]) || !DEFVALID(Def, cube[0][1]) || !DEFVALID(Def, cube[0][2]) || !DEFVALID(Def, cube[0][3])) {
        return Def->NoData;
    }

    // 3D Interpolation case
    if (Z > TINY_VALUE) {
        idx[0] += idxk;
        idx[1] += idxk;
        idx[3] += idxk;
        idx[2] += idxk;
        k++;

        if (Idx == -1) {
            Def_GetQuadMod(Def, idx, cube[1]);
        } else {
            Def_GetQuad(Def, Idx, idx, cube[1]);
        }

        // Check for masked values
        if (Def->Mask) {
            if (!Def->Mask[idx[0]]) cube[1][0] = VertexAvg(Def, Idx, i,  j,  k);
            if (!Def->Mask[idx[1]]) cube[1][1] = VertexAvg(Def, Idx, i+1, j,  k);
            if (!Def->Mask[idx[2]]) cube[1][2] = VertexAvg(Def, Idx, i+1, j+1, k);
            if (!Def->Mask[idx[3]]) cube[1][3] = VertexAvg(Def, Idx, i,  j+1, k);
        }

        // If either value is nodata then interpolation will be nodata as well
        if (!DEFVALID(Def, cube[1][0]) || !DEFVALID(Def, cube[1][1]) || !DEFVALID(Def, cube[1][2]) || !DEFVALID(Def, cube[1][3])) {
            return Def->NoData;
        }

        cube[0][0] = ILIN(cube[0][0], cube[1][0], Z);
        cube[0][1] = ILIN(cube[0][1], cube[1][1], Z);
        cube[0][2] = ILIN(cube[0][2], cube[1][2], Z);
        cube[0][3] = ILIN(cube[0][3], cube[1][3], Z);
    }

    // Interpolate over X
    if (X > TINY_VALUE) {
        cube[0][0] = ILIN(cube[0][0], cube[0][1], X);
        cube[0][3] = ILIN(cube[0][3], cube[0][2], X);
    }

    // Interpolate over Y
    if (Y > TINY_VALUE) {
        cube[0][0] = ILIN(cube[0][0], cube[0][3], Y);
    }

    return cube[0][0];
}


double VertexValV(
    //! [in] Data definition
    const TDef * const Def,
    //! [in]
    double X,
    //! [in]
    double Y,
    //! [in]
    double Z,
    //! [out]
    Vect3d V
) {
    if (X > Def->NI - 1 || Y > Def->NJ - 1 || Z > Def->NK - 1 || X < 0 || Y < 0 || Z < 0) return 0;

    unsigned long i = X; X -= i;
    unsigned long j = Y; Y -= j;
    unsigned long k = Z; Z -= k;

    // Get gridpoint32_t indexes
    unsigned long idxk = Def->NIJ;

    unsigned long idx[4];
    idx[0]=(k?idxk*k:0)+j*Def->NI+i;
    idx[1]=idx[0]+1;
    idx[3]=(j==Def->NJ-1)?idx[0]:idx[0]+Def->NI;
    idx[2]=idx[3]+1;

    double cube[3][2][4];
    Def_GetQuad(Def, 0, idx, cube[0][0]);
    Def_GetQuad(Def, 1, idx, cube[1][0]);
    if (Def->Data[2]) Def_GetQuad(Def, 2, idx, cube[2][0]);

    // If either value is nodata then interpolation will be nodata as well
    if (!DEFVALID(Def, cube[0][0][0]) || !DEFVALID(Def, cube[0][0][1]) || !DEFVALID(Def, cube[0][0][2]) || !DEFVALID(Def, cube[0][0][3])) {
        return(Def->NoData);
    }

    // 3D Interpolation case
    if (Z > TINY_VALUE) {
        idx[0] += idxk;
        idx[1] += idxk;
        idx[3] += idxk;
        idx[2] += idxk;
        Def_GetQuad(Def, 0, idx, cube[0][1]);
        Def_GetQuad(Def, 1, idx, cube[1][1]);
        if (Def->Data[2]) Def_GetQuad(Def, 2, idx, cube[2][1]);

        // If either value is nodata then interpolation will be nodata as well
        if (!DEFVALID(Def, cube[0][1][0]) || !DEFVALID(Def, cube[0][1][1]) || !DEFVALID(Def, cube[0][1][2]) || !DEFVALID(Def, cube[0][1][3])) {
            return(Def->NoData);
        }

        cube[0][0][0] = ILIN(cube[0][0][0], cube[0][1][0], Z);
        cube[0][0][1] = ILIN(cube[0][0][1], cube[0][1][1], Z);
        cube[0][0][2] = ILIN(cube[0][0][2], cube[0][1][2], Z);
        cube[0][0][3] = ILIN(cube[0][0][3], cube[0][1][3], Z);
        cube[1][0][0] = ILIN(cube[1][0][0], cube[1][1][0], Z);
        cube[1][0][1] = ILIN(cube[1][0][1], cube[1][1][1], Z);
        cube[1][0][2] = ILIN(cube[1][0][2], cube[1][1][2], Z);
        cube[1][0][3] = ILIN(cube[1][0][3], cube[1][1][3], Z);
        if (Def->Data[2]) {
            cube[2][0][0] = ILIN(cube[2][0][0], cube[2][1][0], Z);
            cube[2][0][1] = ILIN(cube[2][0][1], cube[2][1][1], Z);
            cube[2][0][2] = ILIN(cube[2][0][2], cube[2][1][2], Z);
            cube[2][0][3] = ILIN(cube[2][0][3], cube[2][1][3], Z);
        }
    }

    V[0] = cube[0][0][0];
    V[1] = cube[1][0][0];
    V[2] = Def->Data[2] ? cube[2][0][0] : 0.0;

    // Interpolate over X
    if (X > TINY_VALUE) {
        V[0] = cube[0][0][0] = ILIN(cube[0][0][0], cube[0][0][1], X);
               cube[0][0][3] = ILIN(cube[0][0][3], cube[0][0][2], X);
        V[1] = cube[1][0][0] = ILIN(cube[1][0][0], cube[1][0][1], X);
               cube[1][0][3] = ILIN(cube[1][0][3], cube[1][0][2], X);
        if (Def->Data[2]) {
            V[2] = cube[2][0][0] = ILIN(cube[2][0][0], cube[2][0][1], X);
                   cube[2][0][3] = ILIN(cube[2][0][3], cube[2][0][2], X);
        }
    }

    // Interpolate over Y
    if (Y > TINY_VALUE) {
        V[0] = ILIN(cube[0][0][0], cube[0][0][3], Y);
        V[1] = ILIN(cube[1][0][0], cube[1][0][3], Y);
        if (Def->Data[2]) {
            V[2] = ILIN(cube[2][0][0], cube[2][0][3], Y);
        }
    }

    return 0;
}


static inline float VertexAvgS(
    const double * const Data,
    const char * const Mask,
    const int32_t NI,
    const int32_t NJ,
    const int32_t X,
    const int32_t Y
) {
    unsigned long n = 0;
    float val = 0.0;
    unsigned long i = X;
    unsigned long j = Y;
    unsigned long idx;
    j-=1;i-=1; idx=j*NI+i; if (j>=0 && i>=0 && Mask[idx]) { val += Data[idx]; n++; }
    i++;       idx++;      if (j>=0         && Mask[idx]) { val += Data[idx]; n++; }
    i++;       idx++;      if (j>=0 && i<NI && Mask[idx]) { val += Data[idx]; n++; }
    j++;       idx+=NI;    if (        i<NI && Mask[idx]) { val += Data[idx]; n++; }
    i-=2;      idx-=2;     if (        i>=0 && Mask[idx]) { val += Data[idx]; n++; }
    j++;       idx+=NI;    if (j<NJ && i>=0 && Mask[idx]) { val += Data[idx]; n++; }
    i++;       idx++;      if (j<NJ         && Mask[idx]) { val += Data[idx]; n++; }
    i++;       idx++;      if (j<NJ && i<NI && Mask[idx]) { val += Data[idx]; n++; }

    return n ? val / n : 0.0;
}


float VertexValS(
    const double * const Data,
    const char * const Mask,
    const int32_t NI,
    const int32_t NJ,
    double X,
    double Y,
    const char Geo
) {
    if (Geo) {
        X = CLAMP(X, 0, NI - 1);
        Y = CLAMP(Y, 0, NJ - 1);
    }

    if (!Data || X > NI - 1 || Y > NJ - 1 || X < 0 || Y < 0) return 0;

    unsigned long i = X; X -= i;
    unsigned long j = Y; Y -= j;

    // Get gridpoint32_t indexes
    unsigned long idx[4];
    idx[0] = j * NI + i;
    idx[1] = idx[0] + 1;
    idx[3] = (j == NJ - 1) ? idx[0] : idx[0] + NI;
    idx[2] = idx[3] + 1;

    double cell[4];
    cell[0] = Data[idx[0]];
    cell[1] = Data[idx[1]];
    cell[2] = Data[idx[2]];
    cell[3] = Data[idx[3]];

    if (Mask) {
        if (!Mask[idx[0]]) cell[0] = VertexAvgS(Data, Mask, NI, NJ, i,  j);
        if (!Mask[idx[1]]) cell[1] = VertexAvgS(Data, Mask, NI, NJ, i+1, j);
        if (!Mask[idx[2]]) cell[2] = VertexAvgS(Data, Mask, NI, NJ, i+1, j+1);
        if (!Mask[idx[3]]) cell[3] = VertexAvgS(Data, Mask, NI, NJ, i,  j+1);
    }

    // Interpolate over X
    if (X > TINY_VALUE) {
        if (Geo) {
            // If interpolation on coordinates, check for wrap-around cell[0] and cell[3] contain longitude
            double d = cell[0] - cell[1];
            if (d >= 180) cell[1] += 360.0;
            if (d < -180) cell[1] -= 360.0;
            d = cell[3] - cell[2];
            if (d >= 180) cell[2] += 360.0;
            if (d < -180) cell[2] -= 360.0;
        }
        cell[0]=ILIN(cell[0], cell[1], X);
        cell[3]=ILIN(cell[3], cell[2], X);
    }

    // Interpolate over Y
    if (Y > TINY_VALUE) {
        if (Geo) {
            // If interpolation on coordinates, check for wrap-around cell[0] and cell[3] contain longitude
            const double d = cell[0] - cell[3];
            if (d >= 180) cell[3] += 360.0;
            if (d < -180) cell[3] -= 360.0;
        }
        cell[0] = ILIN(cell[0], cell[3], Y);
    }

    return cell[0];
}
