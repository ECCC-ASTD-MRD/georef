!/* RMNLIB - Library of useful routines for C and FORTRAN programming
! * Copyright (C) 1975-2001  Division de Recherche en Prevision Numerique
! *                          Environnement Canada
! *
! * This library is free software; you can redistribute it and/or
! * modify it under the terms of the GNU Lesser General Public
! * License as published by the Free Software Foundation,
! * version 2.1 of the License.
! *
! * This library is distributed in the hope that it will be useful,
! * but WITHOUT ANY WARRANTY; without even the implied warranty of
! * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! * Lesser General Public License for more details.
! *
! * You should have received a copy of the GNU Lesser General Public
! * License along with this library; if not, write to the
! * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
! * Boston, MA 02111-1307, USA.
! */
!**s/r ez8_llll2gd - computes the grid co-ordinates of a point
!
      subroutine ez8_llll2gd(x,y,dlat,dlon,npts,xlat0,xlon0, dellat, dellon, lonref)
      implicit none

!*--------------------------------------------------------------------

      integer npts
      real*8 x(npts), y(npts), dlat(npts), dlon(npts), tmplon
      real xlat0, xlon0, dellat, dellon, lonref
      integer i

      do 10 i=1,npts
         tmplon=dlon(i)
         if (lonref.eq.-180.0) then
            if (tmplon.gt.180.0) then
               tmplon = tmplon - 360.0
            endif
         else if (tmplon.lt.0.0) then
               tmplon = tmplon + 360.0
         endif
         x(i) = (tmplon - xlon0)/dellon + 1.0
         y(i) = (dlat(i) - xlat0)/dellat + 1.0
 10   continue

      return
      end
!-----------------------------------------------------------------

