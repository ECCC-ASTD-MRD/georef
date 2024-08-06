#ifndef _Vertex_h
#define _Vertex_h

#include "GeoRef.h"
#include "Def.h"

void   Vertex_Map(Vect2d P[4],double *LX,double *LY,double WX,double WY);
void   VertexGradient(TDef *Def,Vect3d Nr);
float  VertexVal(TDef *Def,int Idx,double X,double Y,double Z);
double VertexValV(TDef *Def,double X,double Y,double Z,Vect3d V);
float  VertexValS(double *Data,char *Mask,int NI,int NJ,double X,double Y,char Geo);
int    VertexLoc(Vect3d **Pos,TDef *Def,Vect3d Vr,double X,double Y,double Z,int Wrap);
void   VertexInterp(Vect3d Pi,Vect3d P0,Vect3d P1,double V0,double V1,double Level);

#endif
