! =========================================================
! Environnement Canada
! Centre Meteorologique Canadien
! 2100 Trans-Canadienne
! Dorval, Quebec
!
! Projet       : Lecture et traitements de fichiers raster
! Fichier      : GeoRef.h.f
! Creation     : July 2020 - J.P. Gauthier
!
! Description  : Fonctions de manipulations de projections.
!
! Remarques    : Fortran header file ISO_C_BINDINGS
!
! License      :
!    This library is free software; you can redistribute it and/or
!    modify it under the terms of the GNU Lesser General Public
!    License as published by the Free Software Foundation,
!    version 2.1 of the License.
!
!    This library is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!    Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public
!    License along with this library; if not, write to the
!    Free Software Foundation, Inc., 59 Temple Place - Suite 330,
!    Boston, MA 02111-1307, USA.
!
!=========================================================

#ifndef _GeoRef_hf
#define _GeoRef_hf

module libgeoref
   use, intrinsic :: ISO_C_BINDING
   implicit none

   enum, bind(C) !< TDef_Type
      enumerator :: TD_Unknown = 0
      enumerator :: TD_Binary  = 1
      enumerator :: TD_UByte   = 2
      enumerator :: TD_Byte    = 3
      enumerator :: TD_UInt16  = 4
      enumerator :: TD_Int16   = 5
      enumerator :: TD_UInt32  = 6
      enumerator :: TD_Int32   = 7
      enumerator :: TD_UInt64  = 8
      enumerator :: TD_Int64   = 9
      enumerator :: TD_Float32 = 10
      enumerator :: TD_Float64 = 11
   end enum

   ! Raster interpolation modes
   enum, bind(C)  !<TDef_InterpR
      enumerator :: IR_NEAREST                        = 0
      enumerator :: IR_LINEAR                         = 1
      enumerator :: IR_CUBIC                          = 2
      enumerator :: IR_NORMALIZED_CONSERVATIVE        = 3
      enumerator :: IR_CONSERVATIVE                   = 4
      enumerator :: IR_MAXIMUM                        = 5
      enumerator :: IR_MINIMUM                        = 6
      enumerator :: IR_SUM                            = 7
      enumerator :: IR_AVERAGE                        = 8
      enumerator :: IR_VARIANCE                       = 9
      enumerator :: IR_SQUARE                         = 10
      enumerator :: IR_NORMALIZED_COUNT               = 11
      enumerator :: IR_COUNT                          = 12
      enumerator :: IR_VECTOR_AVERAGE                 = 13
      enumerator :: IR_NOP                            = 14
      enumerator :: IR_ACCUM                          = 15
      enumerator :: IR_BUFFER                         = 16
      enumerator :: IR_SUBNEAREST                     = 17
      enumerator :: IR_SUBLINEAR                      = 18
   end enum

   ! Raster extrapolation modes
   enum, bind(C) !< TDef_ExtrapR
      enumerator :: ER_UNDEF   = 0
      enumerator :: ER_MAXIMUM = 1
      enumerator :: ER_MINIMUM = 2
      enumerator :: ER_VALUE   = 3
      enumerator :: ER_ABORT   = 4
      enumerator :: ER_NEAREST = 5
      enumerator :: ER_LINEAR  = 6
      enumerator :: ER_CUBIC   = 7
   end enum

!TODO: If we include TDef
!   type, bind(C) :: c_tdef
!      type(C_PTR) :: Buffer,Aux        !< Buffer temporaire
!      type(C_PTR) :: Accum             !< Accumulation Buffer temporaire
!      type(C_PTR) :: Mask              !< Masque a appliquer au traitement sur le champs
!      type(C_PTR) :: Data(4)           !< Composantes du champs (Pointeurs sur les donnees)
!      type(C_PTR) :: Mode              !< Module des champs Data is vectoriel
!      type(C_PTR) :: Dir               !< Direction si vectoriel
!      type(C_PTR) :: Pres,Height       !< Pression au sol
!      type(C_PTR) :: Sub               !< Sub grid resolutions values
!      type(C_PTR) :: Pick,Poly         !< Geometry used in various interpolation method
!      type(C_PTR) :: Segments          !< Liste d'objets de rendue

!      real(C_DOUBLE) :: NoData            !< Valeur de novalue
!      integer(C_INT) :: Type              !< Type de donnees du champs (enum TDef_Type)
!      integer(C_INT) :: NI,NJ,NK,NC,NIJ   !< Dimensions du champs
!      integer(C_INT) :: Idx               !< Index displacement into supergrid

!      integer(C_INT) :: CellDim           !< Defined grid point coverage, point=1 or area=2
!      real(C_DOUBLE) :: CoordLimits(2,2)  !< Limits of processing in latlon
!      integer(C_INT) :: Limits(2,3)       !< Limits of processing in grid points
!      integer(C_INT) :: Level             !< Niveau courant
!      integer(C_INT) :: Sample,SubSample  !< Sample interval in grid points
!      integer(C_SIGNED_CHAR) :: Alias     !< Alias d'un autre TDef (Pointe sur d'autres donnees)
!   end type c_tdef

   type, bind(C) :: c_tgeooptions
      integer(C_INT) :: InterpDegree;   !< Interpolation degree
      integer(C_INT) :: ExtrapDegree;   !< Extrapolation method
      real(C_DOUBLE) :: ExtrapValue;    !< Value to use for extrapolation in ER_VALUE mode
      integer(C_INT) :: SubGrid;        !< Subgrid to use (0=all)
      integer(C_INT) :: Transform;      !< Apply transformation or stay within master referential
      integer(C_INT) :: CIndex;         !< C Indexing (starts st 0)
      integer(C_INT) :: Symmetric;      !< 
      integer(C_INT) :: WeightNum;      !<
      integer(C_CHAR) :: PolarCorrect;  !< Apply polar corrections
      integer(C_CHAR) :: VectorMode;    !< Process data as vector
      real(C_FLOAT) ::  DistTreshold;   !< Distance treshold for point clouds
      real(C_FLOAT) :: LonRef;          !<Longitude referential (-180.0,0.0)
   end type c_tgeooptions

   type(c_tgeooptions), bind(C,name="GeoRef_Options") :: c_georef_options

   interface

!TODO: If we include TDef
!      function c_def_new(ni,nj,nk,dim,type) result(def) bind(C,name='Def_New')
!         import :: C_PTR, C_CHAR, C_INT, c_tdef
!         integer(C_INT), intent(IN), value :: ni
!         integer(C_INT), intent(IN), value :: nj
!         integer(C_INT), intent(IN), value :: nk
!         integer(C_INT), intent(IN), value :: dim
!         integer(C_INT), intent(IN), value :: type  !< enum TDef_Type
!         type(c_tdef) :: def
!      end function c_def_new

!      function c_def_free(def) result(status) bind(C,name='Def_Free')
!         import :: C_INT, c_tdef
!         type(c_tdef), intent(IN) :: def
!         integer(C_INT)           :: status
!      end function c_def_free

!      function c_def_gridinterp(grefout,defout,grefin,defin,interp,extrap,mask,index) result(status) bind(C,name='Def_GridInterp')
!         import :: C_PTR, C_CHAR, C_INT, c_tdef
!         type(C_PTR), intent(IN) :: grefout
!         type(c_tdef), intent(IN) :: defout
!         type(C_PTR), intent(IN) :: grefin
!         type(c_tdef), intent(IN) :: defin
!         integer(C_INT), intent(IN), value :: interp  !< enum TDef_InterpR
!         integer(C_INT), intent(IN), value :: extrap
!         integer(C_INT), intent(IN), value :: mask
!         type(C_PTR), intent(INOUT) :: index
!         integer(C_INT)                :: status
!      end function c_def_gridinterp
      
!      void     GeoRef_GridGetExpanded(TGeoRef *Ref, float *zout, float *zin);                                                           // gdxpngd
!      int      GeoRef_AxisGetExpanded(TGeoRef* Ref, float *AX, float *AY);                                                              // gdgxpndaxes
!      void     GeoRef_AxisDefine(TGeoRef* Ref,float *AX,float *AY);
!      void     GeoRef_AxisCalcExpandCoeff(TGeoRef* Ref);
!      void     GeoRef_AxisCalcNewtonCoeff(TGeoRef* Ref);
      
!      TGeoRef* GeoRef_RPNCreate(int NI,int NJ,char *GRTYP,int IG1,int IG2,int IG3,int IG4,int FID);
      function c_georef_rpncreate(ni,nj,grtyp,ig1,ig2,ig3,ig4,iunit) result(gref) bind(C,name='GeoRef_RPNCreate')
         import :: C_PTR, C_CHAR, C_INT
         integer(C_INT), intent(IN), value :: ni
         integer(C_INT), intent(IN), value :: nj
         character(C_CHAR), dimension(*), intent(IN) :: grtyp
         integer(C_INT), intent(IN), value :: ig1
         integer(C_INT), intent(IN), value :: ig2
         integer(C_INT), intent(IN), value :: ig3
         integer(C_INT), intent(IN), value :: ig4
         integer(C_INT), intent(IN), value :: iunit
         type(C_PTR) :: gref
      end function c_georef_rpncreate

!      TGeoRef* GeoRef_RPNCreateInMemory(int NI,int NJ,char* GRTYP,char* GRREF,int IG1,int IG2,int IG3,int IG4,float* AX,float* AY)
      function c_georef_rpncreateinmemory(ni,nj,grtyp,grref,ig1,ig2,ig3,ig4,ax,ay) result(gref) bind(C,name='GeoRef_RPNCreateInMemory')
         import :: C_PTR, C_CHAR, C_INT
         integer(C_INT), intent(IN), value :: ni
         integer(C_INT), intent(IN), value :: nj
         character(C_CHAR), dimension(*), intent(IN) :: grtyp
         character(C_CHAR), dimension(*), intent(IN) :: grref
         integer(C_INT), intent(IN), value :: ig1
         integer(C_INT), intent(IN), value :: ig2
         integer(C_INT), intent(IN), value :: ig3
         integer(C_INT), intent(IN), value :: ig4
         type(C_PTR), intent(IN), value :: ax
         type(C_PTR), intent(IN), value :: ay
         type(C_PTR) :: gref
      end function c_georef_rpncreateinmemory

!      int      GeoRef_Interp(TGeoRef *RefTo,TGeoRef *RefFrom,float *zout, float *zin);                                                  // c_ezsint
      function c_georef_interp(grefout,grefin,zout,zin) result(status) bind(C,name='GeoRef_Interp') 
         import :: C_PTR, C_INT, C_FLOAT
         type(C_PTR), intent(IN), value :: grefout
         type(C_PTR), intent(IN), value :: grefin
         type(C_PTR), intent(IN), value :: zout
         type(C_PTR), intent(IN), value :: zin
         integer(C_INT)                :: status
      end function c_georef_interp

!      int      GeoRef_InterpUV(TGeoRef *RefTo,TGeoRef *RefFrom,float *uuout,float *vvout,float *uuin,float *vvin);                      // c_ezuvint
      function c_georef_interpuv(grefout,grefin,uuout,vvout,uuin,vvin) result(status) bind(C,name='GeoRef_InterpUV') 
         import :: C_PTR, C_INT
         type(C_PTR), intent(IN), value :: grefout
         type(C_PTR), intent(IN), value :: grefin
         type(C_PTR), intent(IN), value :: uuout
         type(C_PTR), intent(IN), value :: vvout
         type(C_PTR), intent(IN), value :: uuin
         type(C_PTR), intent(IN), value :: vvin
         integer(C_INT)                :: status
      end function c_georef_interpuv

!      int      GeoRef_InterpWD(TGeoRef *RefTo,TGeoRef *RefFrom,float *uuout,float *vvout,float *uuin,float *vvin);                      // c_ezwdint
      function c_georef_interpwd(grefout,grefin,uuout,vvout,uuin,vvin) result(status) bind(C,name='GeoRef_InterpWD') 
         import :: C_PTR, C_INT
         type(C_PTR), intent(IN), value :: grefout
         type(C_PTR), intent(IN), value :: grefin
         type(C_PTR), intent(IN), value :: uuout
         type(C_PTR), intent(IN), value :: vvout
         type(C_PTR), intent(IN), value :: uuin
         type(C_PTR), intent(IN), value :: vvin
         integer(C_INT)                :: status
      end function c_georef_interpwd

!      int      GeoRef_XYVal(TGeoRef *Ref,float *zout,float *zin,double *x,double *y,int n);                                             // c_gdxysval
      function c_georef_xyval(gref,out,in,x,y,n) result(status) bind(C,name='GeoRef_XYVal') 
         import :: C_PTR, C_INT
         type(C_PTR), intent(IN), value :: gref
         type(C_PTR), intent(IN), value :: out
         type(C_PTR), intent(IN), value :: in
         type(C_PTR), intent(IN), value :: x
         type(C_PTR), intent(IN), value :: y
         integer(C_INT), intent(IN), value :: n
         integer(C_INT)                :: status
      end function c_georef_xyval

!      int      GeoRef_XYUVVal(TGeoRef *Ref,float *uuout,float *vvout,float *uuin,float *vvin,double *x,double *y,int n);                // c_gdxyvval
      function c_georef_xyuvval(gref,uuout,vvout,uuin,vvin,x,y,n) result(status) bind(C,name='GeoRef_XYUVVal') 
         import :: C_PTR, C_INT
         type(C_PTR), intent(IN), value :: gref
         type(C_PTR), intent(IN), value :: uuout
         type(C_PTR), intent(IN), value :: vvout
         type(C_PTR), intent(IN), value :: uuin
         type(C_PTR), intent(IN), value :: vvin
         type(C_PTR), intent(IN), value :: x
         type(C_PTR), intent(IN), value :: y
         integer(C_INT), intent(IN), value :: n
         integer(C_INT)                :: status
      end function c_georef_xyuvval

!      int      GeoRef_XYWDVal(TGeoRef *Ref,float *uuout,float *vvout,float *uuin,float *vvin,double *x,double *y,int n);                // c_gdxywdval
      function c_georef_xywdval(gref,uuout,vvout,uuin,vvin,x,y,n) result(status) bind(C,name='GeoRef_XYWDVal') 
         import :: C_PTR, C_INT
         type(C_PTR), intent(IN), value :: gref
         type(C_PTR), intent(IN), value :: uuout
         type(C_PTR), intent(IN), value :: vvout
         type(C_PTR), intent(IN), value :: uuin
         type(C_PTR), intent(IN), value :: vvin
         type(C_PTR), intent(IN), value :: x
         type(C_PTR), intent(IN), value :: y
         integer(C_INT), intent(IN), value :: n
         integer(C_INT)                :: status
      end function c_georef_xywdval

!      int      GeoRef_LLVal(TGeoRef *Ref,float *zout,float *zin,double *lat,double *lon,int n);                                         // c_gdllsval
      function c_georef_llval(gref,out,in,lat,lon,n) result(status) bind(C,name='GeoRef_LLVal') 
         import :: C_PTR, C_INT
         type(C_PTR), intent(IN), value :: gref
         type(C_PTR), intent(IN), value :: out
         type(C_PTR), intent(IN), value :: in
         type(C_PTR), intent(IN), value :: lat
         type(C_PTR), intent(IN), value :: lon
         integer(C_INT), intent(IN), value :: n
         integer(C_INT)                :: status
      end function c_georef_llval

!      int      GeoRef_LLUVVal(TGeoRef *Ref,float *uuout,float *vvout,float *uuin,float *vvin,double *Lat,double *Lon,int n);            // c_gdllvval
      function c_georef_lluvval(gref,uuout,vvout,uuin,vvin,lat,lon,n) result(status) bind(C,name='GeoRef_LLUVVal') 
         import :: C_PTR, C_INT
         type(C_PTR), intent(IN), value :: gref
         type(C_PTR), intent(IN), value :: uuout
         type(C_PTR), intent(IN), value :: vvout
         type(C_PTR), intent(IN), value :: uuin
         type(C_PTR), intent(IN), value :: vvin
         type(C_PTR), intent(IN), value :: lat
         type(C_PTR), intent(IN), value :: lon
         integer(C_INT), intent(IN), value :: n
         integer(C_INT)                :: status
      end function c_georef_lluvval

!      int      GeoRef_LLWDVal(TGeoRef *Ref,float *uuout,float *vvout,float *uuin,float *vvin,double *Lat,double *Lon,int n);            // c_gdllwdval
      function c_georef_llwdval(gref,uuout,vvout,uuin,vvin,lat,lon,n) result(status) bind(C,name='GeoRef_LLWDVal') 
         import :: C_PTR, C_INT
         type(C_PTR), intent(IN), value :: gref
         type(C_PTR), intent(IN), value :: uuout
         type(C_PTR), intent(IN), value :: vvout
         type(C_PTR), intent(IN), value :: uuin
         type(C_PTR), intent(IN), value :: vvin
         type(C_PTR), intent(IN), value :: lat
         type(C_PTR), intent(IN), value :: lon
         integer(C_INT), intent(IN), value :: n
         integer(C_INT)                :: status
      end function c_georef_llwdval

!      int      GeoRef_UV2WD(TGeoRef *Ref,float *spd_out,float *wd_out,float *uuin,float *vvin,double *Lat,double *Lon,int npts);        // c_gdwdfuv
      function c_georef_uv2wd(gref,spdout,wdout,uuin,vvin,lat,lon,npts) result(status) bind(C,name='GeoRef_UV2WD') 
         import :: C_PTR, C_INT
         type(C_PTR), intent(IN), value :: gref
         type(C_PTR), intent(IN), value :: spdout
         type(C_PTR), intent(IN), value :: wdout
         type(C_PTR), intent(IN), value :: uuin
         type(C_PTR), intent(IN), value :: vvin
         type(C_PTR), intent(IN), value :: lat
         type(C_PTR), intent(IN), value :: lon
         integer(C_INT), intent(IN), value :: npts
         integer(C_INT)                :: status
      end function c_georef_uv2wd

!      int      GeoRef_LL2XY(TGeoRef *Ref,double *X,double *Y,double *Lat,double *Lon,int n);                                            // c_gdxyfll
      function c_georef_ll2xy(gref,x,y,lat,lon,npts) result(status) bind(C,name='GeoRef_LL2XY') 
         import :: C_PTR, C_INT
         type(C_PTR), intent(IN), value :: gref
         type(C_PTR), intent(IN), value :: x
         type(C_PTR), intent(IN), value :: y
         type(C_PTR), intent(IN), value :: lat
         type(C_PTR), intent(IN), value :: lon
         integer(C_INT), intent(IN), value :: npts
         integer(C_INT)                :: status
      end function c_georef_ll2xy

!      int      GeoRef_CalcLL(TGeoRef* Ref);                                                                                             // ez_calclatlon
      function c_georef_calcll(gref) result(status) bind(C,name='GeoRef_CalcLL') 
         import :: C_PTR, C_INT
         type(C_PTR), intent(IN), value :: gref
         integer(C_INT)                :: status
      end function c_georef_calcll

      !      int      GeoRef_GetLL(TGeoRef *Ref,double *Lat,double *Lon);                                                                      // gdll
      function c_georef_getll(gref,lat,lon) result(status) bind(C,name='GeoRef_GetLL') 
         import :: C_PTR, C_INT
         type(C_PTR), intent(IN), value :: gref
         type(C_PTR), intent(IN), value :: lat
         type(C_PTR), intent(IN), value :: lon
         integer(C_INT)                :: status
      end function c_georef_getll

!      int      GeoRef_Equal(TGeoRef* __restrict const Ref0,TGeoRef* __restrict const Ref1);
      function c_georef_equal(gref0,gref1) result(status) bind(C,name="GeoRef_Equal")
         import :: C_PTR, C_INT
         type(C_PTR), intent(IN), value :: gref0
         type(C_PTR), intent(IN), value :: gref1
         integer(C_INT)                :: status
      end function c_georef_equal

!      int      GeoRef_GridGetParams(TGeoRef *Ref,int *NI,int *NJ,char *GRTYP,int *IG1,int *IG2,int *IG3,int *IG4,char *GRREF,int *IG1REF,int *IG2REF,int *IG3REF,int *IG4REF)  //c_ezgxprm
      function c_georef_gridgetparams(gref,ni,nj,grtyp,ig1,ig2,ig3,ig4,grref,ig1ref,ig2ref,ig3ref,ig4ref) result(status) bind(C,name='GeoRef_GridGetParams') 
         import :: C_PTR, C_INT, C_CHAR
         type(C_PTR), intent(IN), value :: gref
         integer(C_INT), intent(OUT) :: ni
         integer(C_INT), intent(OUT) :: nj
         character(C_CHAR), intent(OUT) :: grtyp(*)
         integer(C_INT), intent(OUT) :: ig1
         integer(C_INT), intent(OUT) :: ig2
         integer(C_INT), intent(OUT) :: ig3
         integer(C_INT), intent(OUT) :: ig4
         character(C_CHAR), intent(OUT) :: grref(*)
         integer(C_INT), intent(OUT) :: ig1ref
         integer(C_INT), intent(OUT) :: ig2ref
         integer(C_INT), intent(OUT) :: ig3ref
         integer(C_INT), intent(OUT) :: ig4ref
         integer(C_INT)                :: status
      end function c_georef_gridgetparams

!      int      GeoRef_AxisGet(TGeoRef *Ref,float *AX,float *AY);                                                                        // gdaxes
      function c_georef_axisget(gref,ax,ay) result(status) bind(C,name="GeoRef_AxisGet")
         import :: C_PTR, C_INT
         type(C_PTR), intent(IN), value :: gref
         type(C_PTR), intent(IN), value :: ax
         type(C_PTR), intent(IN), value :: ay
         integer(C_INT)                :: status
      end function c_georef_axisget

!      int      GeoRef_WD2UV(TGeoRef *Ref,float *uugdout,float *vvgdout,float *uullin,float *vvllin,double *Lat,double *Lon,int npts);   // c_gduvfwd
      function c_georef_wd2uv(gref,uuout,vvout,uuin,vvin,lat,lon,n) result(status) bind(C,name='GeoRef_WD2UV') 
         import :: C_PTR, C_INT
         type(C_PTR), intent(IN), value :: gref
         type(C_PTR), intent(IN), value :: uuout
         type(C_PTR), intent(IN), value :: vvout
         type(C_PTR), intent(IN), value :: uuin
         type(C_PTR), intent(IN), value :: vvin
         type(C_PTR), intent(IN), value :: lat
         type(C_PTR), intent(IN), value :: lon
         integer(C_INT), intent(IN), value :: n
         integer(C_INT)                :: status
      end function c_georef_wd2uv

   end interface 

end module libgeoref
