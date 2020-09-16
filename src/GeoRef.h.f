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

   type, bind(C) :: c_tdef
      type(C_PTR) :: Buffer,Aux        !< Buffer temporaire
      type(C_PTR) :: Accum             !< Accumulation Buffer temporaire
      type(C_PTR) :: Mask              !< Masque a appliquer au traitement sur le champs
      type(C_PTR) :: Data(4)           !< Composantes du champs (Pointeurs sur les donnees)
      type(C_PTR) :: Mode              !< Module des champs Data is vectoriel
      type(C_PTR) :: Dir               !< Direction si vectoriel
      type(C_PTR) :: Pres,Height       !< Pression au sol
      type(C_PTR) :: Sub               !< Sub grid resolutions values
      type(C_PTR) :: Pick,Poly         !< Geometry used in various interpolation method
      type(C_PTR) :: Segments          !< Liste d'objets de rendue

      real(C_DOUBLE) :: NoData            !< Valeur de novalue
      integer(C_INT) :: Type              !< Type de donnees du champs (enum TDef_Type)
      integer(C_INT) :: NI,NJ,NK,NC,NIJ   !< Dimensions du champs
      integer(C_INT) :: Idx               !< Index displacement into supergrid

      integer(C_INT) :: CellDim           !< Defined grid point coverage, point=1 or area=2
      real(C_DOUBLE) :: CoordLimits(2,2)  !< Limits of processing in latlon
      integer(C_INT) :: Limits(2,3)       !< Limits of processing in grid points
      integer(C_INT) :: Level             !< Niveau courant
      integer(C_INT) :: Sample,SubSample  !< Sample interval in grid points
      integer(C_SIGNED_CHAR) :: Alias     !< Alias d'un autre TDef (Pointe sur d'autres donnees)
   end type c_tdef

   interface

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

      function c_def_new(ni,nj,nk,dim,type) result(def) bind(C,name='Def_New')
         import :: C_PTR, C_CHAR, C_INT, c_tdef
         integer(C_INT), intent(IN), value :: ni
         integer(C_INT), intent(IN), value :: nj
         integer(C_INT), intent(IN), value :: nk
         integer(C_INT), intent(IN), value :: dim
         integer(C_INT), intent(IN), value :: type  !< enum TDef_Type
         type(c_tdef) :: def
      end function c_def_new

      function c_def_free(def) result(status) bind(C,name='Def_Free')
         import :: C_INT, c_tdef
         type(c_tdef), intent(IN) :: def
         integer(C_INT)           :: status
      end function c_def_free

      function c_def_gridinterp(grefout,defout,grefin,defin,interp,extrap,mask,index) result(status) bind(C,name='Def_GridInterp')
         import :: C_PTR, C_CHAR, C_INT, c_tdef
         type(C_PTR), intent(IN) :: grefout
         type(c_tdef), intent(IN) :: defout
         type(C_PTR), intent(IN) :: grefin
         type(c_tdef), intent(IN) :: defin
         integer(C_INT), intent(IN), value :: interp  !< enum TDef_InterpR
         integer(C_INT), intent(IN), value :: extrap
         integer(C_INT), intent(IN), value :: mask
         type(C_PTR), intent(INOUT) :: index
         integer(C_INT)                :: status
      end function c_def_gridinterp

   end interface 

end module libgeoref
