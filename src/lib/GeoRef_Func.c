#include <App.h>
#include "GeoRef.h"

/**----------------------------------------------------------------------------
 * @brief  Calculate distance between to location in gridpoint32_t coordinates
 * @date   March 2007
 *    @param[in]  Ref   Pointer to geo reference
 *    @param[in]  X0    First  X coordinate 
 *    @param[in]  Y0    First  Y coordinate 
 *    @param[in]  X1    Second X coordinate 
 *    @param[in]  Y1    Second Y coordinate 
 *
 *    @return           Distance in meters
*/
double GeoRef_GridDistance(TGeoRef *Ref,double X0,double Y0,double X1,double Y1) {

   char *unit,geo;
   double i[2],j[2],lat[2],lon[2],u;

   // Check for unit type 
   geo=1;
#ifdef HAVE_GDAL
   if (Ref->Spatial) {
      u=OSRGetLinearUnits(Ref->Spatial,&unit);
      geo=(unit[0]!='M' && unit[0]!='m');
   }
#endif

   if (geo) {
      // Grid is geographic (latlon)
      i[0]=X0+1.0;
      j[0]=Y0+1.0;
      i[1]=X1+1.0;
      j[1]=Y1+1.0;

      GeoRef_XY2LL(Ref,lat,lon,i,j,2,TRUE);

      X0=DEG2RAD(lon[0]);
      X1=DEG2RAD(lon[1]);
      Y0=DEG2RAD(lat[0]);
      Y1=DEG2RAD(lat[1]);

      return(DIST(0.0,Y0,X0,Y1,X1));
   } else {
      // Grid is meters (utm)
      if (Ref->Transform) {
         i[0]=Ref->Transform[0]+Ref->Transform[1]*X0+Ref->Transform[2]*Y0;
         j[0]=Ref->Transform[3]+Ref->Transform[4]*X0+Ref->Transform[5]*Y0;
         i[1]=Ref->Transform[0]+Ref->Transform[1]*X1+Ref->Transform[2]*Y1;
         j[1]=Ref->Transform[3]+Ref->Transform[4]*X1+Ref->Transform[5]*Y1;
      } else {
         i[0]=X0;
         j[0]=Y0;
         i[1]=X1;
         j[1]=Y1;
      }
      return(hypot(j[1]-j[0],i[1]-i[0])*u);
   }
   return(hypot(X1-X0,Y1-Y0));
}

/**----------------------------------------------------------------------------
 * @brief  Interpolate the position of a point32_t on a great circle
 * @date   February 2008
 *    @param[in]  C1       First point32_t coordinate
 *    @param[in]  C2       Second point32_t coordinate
 *    @param[in]  C3       Point on great circle to locate on great circle
 *
 *    @return              Ratio of distances
*/
double GeoFunc_RadialPointRatio(TCoord C1,TCoord C2,TCoord C3) {

   TCoord cr;
   double d0,d1,d2;

   GeoFunc_RadialPointOn(C1,C2,C3,&cr);

   d0=DIST(0,C1.Lat,C1.Lon,C2.Lat,C2.Lon);
   d1=DIST(0,C1.Lat,C1.Lon,cr.Lat,cr.Lon);
   d2=DIST(0,C2.Lat,C2.Lon,cr.Lat,cr.Lon);

   if(d2>d0) {
      return(-(d2-d0)/d0);
   } else {
      return(d1/d0);
   }
}

/**----------------------------------------------------------------------------
 * @brief  Calculate intersection point32_t of a right angle from a point32_t on a great circle
 * @date   February 2008
 *    @param[in]  C1       First point32_t coordinate
 *    @param[in]  C2       Second point32_t coordinate
 *    @param[in]  C3       Point on great circle to locate on great circle
 *    @param[out] CR       Located point32_t coordinate
 *
 *    @return              Intersection exists (1=Yes, 0=No)
*/
int32_t GeoFunc_RadialPointOn(TCoord C1,TCoord C2,TCoord C3,TCoord *CR) {

   double crs12,crs13,crs3x;

   /*Calculates 90 degree course crossing*/
   crs12=COURSE(C1.Lat,C1.Lon,C2.Lat,C2.Lon);
   crs13=COURSE(C1.Lat,C1.Lon,C3.Lat,C3.Lon);
   crs3x=crs13>crs12?crs12-M_PI2:crs12+M_PI2;

   return(GeoFunc_RadialIntersect(C1,C3,crs12,crs3x,CR));
}

/**----------------------------------------------------------------------------
 * @brief  Calculate intersection point32_t of 2 great circle
 * @date   February 2008
 *    @param[in]  C1       First point32_t coordinate
 *    @param[in]  C2       Second point32_t coordinate
 *    @param[in]  CRS13    Direction between first and third point
 *    @param[in]  CRS23    Direction between second and third point
 *    @param[out] C3       Intersection poiny
 *
 *    @return              Intersection exists (1=Yes, 0=No)
*/
int32_t GeoFunc_RadialIntersect(TCoord C1,TCoord C2,double CRS13,double CRS23,TCoord *C3) {

   double dst13,dst12,crs12,crs21,ang1,ang2,ang3;
   TCoord sinc2,cosc2,sinc1,cosc1;

   sinc1.Lat=sin(C1.Lat);sinc1.Lon=sin(C1.Lon);
   cosc1.Lat=cos(C1.Lat);cosc1.Lon=cos(C1.Lon);
   sinc2.Lat=sin(C2.Lat);sinc2.Lon=sin(C2.Lon);
   cosc2.Lat=cos(C2.Lat);cosc2.Lon=cos(C2.Lon);

   dst12=2*asin(sqrt(pow((sin((C1.Lat-C2.Lat)/2)),2) + cosc1.Lat*cosc2.Lat*pow(sin((C1.Lon-C2.Lon)/2),2)));

   if (sin(C2.Lon-C1.Lon)<0) {
      crs12=acos((sinc2.Lat-sinc1.Lat*cos(dst12))/(sin(dst12)*cosc1.Lat));
   } else {
      crs12=2.0*M_PI-acos((sinc2.Lat-sinc1.Lat*cos(dst12))/(sin(dst12)*cosc1.Lat));
   }

   if (sin(C1.Lon-C2.Lon)<0) {
      crs21=acos((sinc1.Lat-sinc2.Lat*cos(dst12))/(sin(dst12)*cosc2.Lat));
   } else {
      crs21=M_2PI-acos((sinc1.Lat-sinc2.Lat*cos(dst12))/(sin(dst12)*cosc2.Lat));
   }

   ang1=fmod(CRS13-crs12+M_PI,M_2PI)-M_PI;
   ang2=fmod(crs21-CRS23+M_PI,M_2PI)-M_PI;

   if (sin(ang1)*sin(ang2)<=sqrt(10e-15)) {
      /*No intersection*/
      return(0);
   } else {
      ang1=fabs(ang1);
      ang2=fabs(ang2);
      ang3=acos(-cos(ang1)*cos(ang2)+sin(ang1)*sin(ang2)*cos(dst12));
      dst13=asin(sin(ang2)*sin(dst12)/sin(ang3));
      C3->Lat=asin(sinc1.Lat*cos(dst13)+cosc1.Lat*sin(dst13)*cos(CRS13));
      C3->Lon=fmod(C1.Lon-asin(sin(CRS13)*sin(dst13)/cos(C3->Lat))+M_PI,M_2PI)-M_PI;
   }

   return(1);
}

//TODO: include following ezscint funcs into GeoFunc logic 
void c_ez_calcdist2(double *distance, float lat1, float lon1, float lat2, float lon2)
{
   double degre_a_radian;
   double radlat1, radlat2, radlon1, radlon2;
   double dist;
   double earth_radius = 6370997.;
   double a,c,dlat,dlon,sindlat, sindlon;

   degre_a_radian = M_PI / 180.0;

   radlat1 = lat1 * degre_a_radian;
   radlat2 = lat2 * degre_a_radian;
   radlon1 = lon1 * degre_a_radian;
   radlon2 = lon2 * degre_a_radian;

   dlon = radlon2 - radlon1;
   dlat = radlat2 - radlat1;
   sindlat = sin(dlat*0.5);
   sindlon = sin(dlon*0.5);
   a = sindlat*sindlat + cos(radlat1) * cos(radlat2) * sindlon*sindlon;
   c = 2.0 * atan2(sqrt(a),sqrt(1.0-a));
   c = 2.0 * asin(sqrt(a));
   *distance = (float)(c*earth_radius);

/*   *distance = (float)(earth_radius*acos(cos(radlat1)*cos(radlat2)*cos(radlon1-radlon2)+sin(radlat1)*sin(radlat2)));*/
}
/*
  c_ez_calcarea :
   This function computes the area of the solid rectangle formed by  2 latlon points on the sphere.
   Source of the formula : http://mathworld.wolfram.com/GreatCircle.html
                           http://mathworld.wolfram.com/SphericalTrigonometry.html
                           http://mathworld.wolfram.com/SphericalTriangle.html
   Latitudes and longitudes are in degrees.

   The computation is done by splitting the solid rectangle formed by the latlon points into 2 triangles,
   compute the area of each triangle and add the two areas

   If we have
                                    seg_d
                     *--------------------------------------* (lat2, lon2)
                     *                         -----------  *
             seg_e   *          ---------------             * seg_b
                     *----------     seg_c                  *
      (lat1, lon1)   *--------------------------------------*
                                     seg_a

  We have the corresponding angles for triangle 1

                                                          a * (lat2, lon2)
                                               -----------  *
                                ---------------             *
                     *----------                         c  *
      (lat1, lon1)   *-b -----------------------------------*

  and for triangle 2


                     *------------------------------------e-* (lat2, lon2)
                     * c                       -----------  *
                     *          ---------------
                     *-d--------
      (lat1, lon1)   *
*/



void c_ez_calcarea_rect(float *area, float lat1, float lon1, float lat2, float lon2)
{
   double degre_a_radian;
   double radlat1, radlat2, radlon1, radlon2;
   double dist;
   double earth_radius = 6370997.;
   double a,b,c,d,e;
   double seg_a, seg_b, seg_c, seg_d, seg_e;
   double area1, area2;

   degre_a_radian = M_PI / 180.0;

   radlat1 = lat1 * degre_a_radian;
   radlat2 = lat2 * degre_a_radian;
   radlon1 = lon1 * degre_a_radian;
   radlon2 = lon2 * degre_a_radian;

   /*
   seg_a = distance pt(1,1) - pt(2,1)
   seg_b = distance pt(2,1) - pt(2,2)
   seg_c = distance pt(1,1) - pt(2,2)
   seg_d = distance pt(1,2) - pt(2,2)
   seg_e = distance pt(1,1) - pt(1,2)
   */

   /* These computations have not been optimized */

   seg_a = acos(cos(radlat1)*cos(radlat1)*cos(radlon1-radlon2)+sin(radlat1)*sin(radlat1));
   seg_b = acos(cos(radlat1)*cos(radlat2)*cos(radlon2-radlon2)+sin(radlat1)*sin(radlat2));
   seg_c = acos(cos(radlat1)*cos(radlat2)*cos(radlon1-radlon2)+sin(radlat1)*sin(radlat2));
   seg_d = acos(cos(radlat2)*cos(radlat2)*cos(radlon1-radlon2)+sin(radlat2)*sin(radlat2));
   seg_e = acos(cos(radlat1)*cos(radlat2)*cos(radlon1-radlon1)+sin(radlat1)*sin(radlat2));
/*   printf("seg_a:%f seg_b:%f seg_c:%f seg_d:%f seg_e:%f\n", seg_a, seg_b, seg_c, seg_d, seg_e);*/

   c = acos((cos(seg_c)-cos(seg_b)*cos(seg_a))/(sin(seg_b)*sin(seg_a)));
   b = asin(sin(seg_b)*sin(c)/sin(seg_c));
   a = asin(sin(seg_a)*sin(c)/sin(seg_c));
   area1 = ((a + b + c) - M_PI) * earth_radius * earth_radius;
/*   printf("a:%f b:%f c:%f area:%f\n", a, b, c, area1);*/

/*   printf("seg_a:%f seg_b:%f seg_c:%f seg_d:%f seg_e:%f\n", seg_a, seg_b, seg_c, seg_d, seg_e);*/
   c = acos((cos(seg_c)-cos(seg_e)*cos(seg_d))/(sin(seg_e)*sin(seg_d)));
   e = asin(sin(seg_e)*sin(c)/sin(seg_c));
   d = asin(sin(seg_d)*sin(c)/sin(seg_c));
   area2 = ((d + e + c) - M_PI) * earth_radius * earth_radius;
/*   printf("d:%f e:%f c:%f area:%f\n", d, e, c, area2);*/
   *area = (float)(area1 + area2);

}



void c_ez_calcarea(float *area, float lats[], float lons[])
{
   double degre_a_radian;
   double radlat[4], radlon[4];
   double dist;
   double earth_radius = 6370997.;
   double a,b,c,d,e;
   double seg_a, seg_b, seg_c, seg_d, seg_e;
   double area1, area2;
   int32_t i;

   degre_a_radian = M_PI / 180.0;

   for (i=0; i < 4; i++)
      {
      radlat[i] = lats[i] * degre_a_radian;
      radlon[i] = lons[i] * degre_a_radian;
      }


   /*
   seg_a = distance pt(1,1) - pt(2,1)
   seg_b = distance pt(2,1) - pt(2,2)
   seg_c = distance pt(1,1) - pt(2,2)
   seg_d = distance pt(1,2) - pt(2,2)
   seg_e = distance pt(1,1) - pt(1,2)
   */

   /* These computations have not been optimized */

   seg_a = acos(cos(radlat[0])*cos(radlat[1])*cos(radlon[0]-radlon[1])+sin(radlat[0])*sin(radlat[1]));
   seg_b = acos(cos(radlat[1])*cos(radlat[2])*cos(radlon[1]-radlon[2])+sin(radlat[1])*sin(radlat[2]));
   seg_c = acos(cos(radlat[0])*cos(radlat[2])*cos(radlon[0]-radlon[2])+sin(radlat[0])*sin(radlat[2]));
   seg_d = acos(cos(radlat[2])*cos(radlat[3])*cos(radlon[2]-radlon[3])+sin(radlat[2])*sin(radlat[3]));
   seg_e = acos(cos(radlat[0])*cos(radlat[3])*cos(radlon[0]-radlon[3])+sin(radlat[0])*sin(radlat[3]));
/*   printf("seg_a:%f seg_b:%f seg_c:%f seg_d:%f seg_e:%f\n", seg_a, seg_b, seg_c, seg_d, seg_e);*/

   c = acos((cos(seg_c)-cos(seg_b)*cos(seg_a))/(sin(seg_b)*sin(seg_a)));
   b = asin(sin(seg_b)*sin(c)/sin(seg_c));
   a = asin(sin(seg_a)*sin(c)/sin(seg_c));
   area1 = ((a + b + c) - M_PI) * earth_radius * earth_radius;

/*   printf("seg_a:%f seg_b:%f seg_c:%f seg_d:%f seg_e:%f\n", seg_a, seg_b, seg_c, seg_d, seg_e);*/
   c = acos((cos(seg_c)-cos(seg_e)*cos(seg_d))/(sin(seg_e)*sin(seg_d)));
   e = asin(sin(seg_e)*sin(c)/sin(seg_c));
   d = asin(sin(seg_d)*sin(c)/sin(seg_c));
   area2 = ((d + e + c) - M_PI) * earth_radius * earth_radius;
   *area = (float)(area1 + area2);
   if (*area < 0.0)
      {
//      printf("a:%f b:%f c:%f area:%f\n", a, b, c, area1);
//      printf("d:%f e:%f c:%f area:%f\n", d, e, c, area2);
//      printf("area:%f\n", *area);
      }

}

void c_ez_calcarea2(double *area, float lats[], float lons[])
{
   double degre_a_radian;
   double radlat[4], radlon[4];
   double dist;
   double earth_radius = 6370997.;
   double a,b,c,d,e;
   double seg_a, seg_b, seg_c, seg_d, seg_e;
   double area1, area2;
   double s, excess;
   int32_t i;

   degre_a_radian = M_PI / 180.0;

   for (i=0; i < 4; i++)
      {
      radlat[i] = lats[i] * degre_a_radian;
      radlon[i] = lons[i] * degre_a_radian;
      }


   /*
   seg_a = distance pt(1,1) - pt(2,1)
   seg_b = distance pt(2,1) - pt(2,2)
   seg_c = distance pt(1,1) - pt(2,2)
   seg_d = distance pt(1,2) - pt(2,2)
   seg_e = distance pt(1,1) - pt(1,2)
   */

   /* These computations have not been optimized */

   c_ez_calcdist2(&seg_a, lats[0], lons[0], lats[1], lons[1]);
   c_ez_calcdist2(&seg_b, lats[1], lons[1], lats[2], lons[2]);
   c_ez_calcdist2(&seg_c, lats[0], lons[0], lats[2], lons[2]);
   c_ez_calcdist2(&seg_d, lats[2], lons[2], lats[3], lons[3]);
   c_ez_calcdist2(&seg_e, lats[3], lons[3], lats[0], lons[0]);
/*
   printf("seg_a:%f seg_b:%f seg_c:%f seg_d:%f seg_e:%f\n", seg_a, seg_b, seg_c, seg_d, seg_e);
*/

   seg_a /= earth_radius;
   seg_b /= earth_radius;
   seg_c /= earth_radius;
   seg_d /= earth_radius;
   seg_e /= earth_radius;

   s = 0.5*(seg_a+seg_b+seg_c);
   excess = sqrt(tan(0.5*s)*tan(0.5*(s-seg_a))*tan(0.5*(s-seg_b))*tan(0.5*(s-seg_c)));
   excess = atan(excess);
   excess = 4.0 * excess;
   area1 = earth_radius * earth_radius * excess;
/*   printf("s: %f excess: %f area1: %f\n", s, excess, area1);*/

   s = 0.5*(seg_c+seg_d+seg_e);
   excess = sqrt(tan(0.5*s)*tan(0.5*(s-seg_c))*tan(0.5*(s-seg_d))*tan(0.5*(s-seg_e)));
   excess = atan(excess);
   excess *= 4.0;
   area2 = earth_radius * earth_radius * excess;
/*   printf("s: %f excess: %f area2: %f\n", s, excess, area2);*/

   *area = (float)(area1 + area2);
   if (*area < 0.0)
      {
//      printf("area1:%f\n", area1);
//      printf("area2:%f\n", area2);
      }
}

