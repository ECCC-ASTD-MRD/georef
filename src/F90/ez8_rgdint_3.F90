      subroutine ez8_rgdint_3(zo,px,py,npts,z,ni,j1,j2,wrap,nodata)
!*******
!Auteur: Y.Chartier, drpn
!        Fevrier 1991
!
!Objet:  Interpolation bi-cubique de points a partir d'une grille
!        source reguliere.
!        
!*******
      implicit none

      integer npts,ni,j1,j2,wrap
      real    zo(npts)
      real(kind=8)  px(npts),py(npts)
      real    z(ni,j1:j2),nodata
!
!  npts   : nombre de points a interpoler
!  i1:i2  : dimension de la grille source selon x
!  j1:nj  : dimension de la grille source selon y
!  zo     : vecteur de sortie contenant les valeurs interpolees
!  px     : vecteur contenant la position x des points que l'on veut interpoler
!  py     : vecteur contenant la position y des points que l'on veut interpoler
!  z      : valeurs de la grille source.
!
!  wrap est est le facteur de "wrap around" dans le cadre d'une grille globale
!  pour une grille de type 'A' ou 'G', wrap = 2
!  pour une grille de type 'B', wrap = 1
!  dans tous les autres case wrap = 0
!
!===========================================
!
!     *   *   *   *
!     
!     *   *   *   *
!           #        ==>   pt (x,y)
!     *  (=)  *   *  ==> = pt (iind, jind)
!
!     *   *   *   *
!
!===========================================
      real(kind=8)  y1,y2,y3,y4
      integer n,i,j
      integer imoins1,iplus1,iplus2,limite

#include "cubic8.cdk"

      limite  = ni+2-wrap
 
      !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(n,i,j,imoins1,iplus1,iplus2,dx,dy,y1,y2,y3,y4) SHARED(npts,j1,j2,ni,px,py,wrap,z,zo,nodata,limite)
      do n=1,npts
         i = min(ni-2+wrap,max(1,max(2-wrap,int(px(n)))))
         j = min(j2-2,max(j1+1,int(py(n))))
         
         if (wrap.gt.0.and.(i.le.1).or.i.ge.(ni-1)) then
            imoins1 = mod(limite+i-1,limite)
            iplus1  = mod(limite+i+1,limite)
            iplus2  = mod(limite+i+2,limite)
            
            if (imoins1.eq.0) imoins1 = ni
            if (i.eq.0) i = ni
            if (iplus1.eq.0) iplus1 = ni
            if (iplus2.eq.0) iplus2 = ni
            
            if (wrap.eq.1) then
               if (iplus2.eq.ni) iplus2 = 2
               if (imoins1.eq.ni) imoins1=ni-1
            endif
         else
            imoins1 = i-1
            iplus1  = i+1
            iplus2  = i+2
         endif
         dx = px(n) - i
         dy = py(n) - j
         
         zo(n)=nodata

         if (defvalid(z(imoins1,j-1),nodata) .and. defvalid(z(i ,j-1),nodata) .and. defvalid(z(iplus1,j-1),nodata) .and. defvalid(z(iplus2,j-1),nodata)) then        
            y1=cubic(dble(z(imoins1,j-1)),dble(z(i ,j-1)),dble(z(iplus1,j-1)),dble(z(iplus2,j-1)),dx)
            if (defvalid(z(imoins1,j),nodata) .and. defvalid(z(i,j),nodata) .and. defvalid(z(iplus1,j),nodata) .and. defvalid(z(iplus2,j),nodata)) then        
               y2=cubic(dble(z(imoins1,j)),dble(z(i ,j  )),dble(z(iplus1,j  )),dble(z(iplus2,j  )),dx)
               if (defvalid(z(imoins1,j+1),nodata) .and. defvalid(z(i ,j+1),nodata) .and. defvalid(z(iplus1,j+1),nodata) .and. defvalid(z(iplus2,j+1),nodata)) then        
                  y3=cubic(dble(z(imoins1,j+1)),dble(z(i ,j+1)),dble(z(iplus1,j+1)),dble(z(iplus2,j+1)),dx)
                  if (defvalid(z(imoins1,j+2),nodata) .and. defvalid(z(i ,j+2),nodata) .and. defvalid(z(iplus1,j+2),nodata) .and. defvalid(z(iplus2,j+2),nodata)) then        
                     y4=cubic(dble(z(imoins1,j+2)),dble(z(i ,j+2)),dble(z(iplus1,j+2)),dble(z(iplus2,j+2)),dx)
                  endif
               endif
            endif
         endif
         
         zo(n)=cubic(y1,y2,y3,y4,dy)
      enddo
      !$OMP END PARALLEL DO 
      
      return
      end

      subroutine ez8_rgd_index_3(index,px,py,npts,ni,j1,j2,wrap)
         implicit none
   
         integer       npts,ni,j1,j2,wrap
         real(kind=4)  index(npts,10)
         real(kind=8)  px(npts),py(npts)
         real(kind=8)  dx,dy,y1,y2,y3,y4
         integer       n,i,j
         integer       imoins1,iplus1,iplus2,limite
      
         limite  = ni+2-wrap
      
         !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(n,i,j,imoins1,iplus1,iplus2,dx,dy,y1,y2,y3,y4) SHARED(npts,index,j1,j2,ni,px,py,wrap,limite)
         do n=1,npts
            i = min(ni-2+wrap,max(1,max(2-wrap,int(px(n)))))
            j = min(j2-2,max(j1+1,int(py(n))))
            
            if (wrap.gt.0) then
               imoins1 = mod(limite+i-1,limite)
               iplus1  = mod(limite+i+1,limite)
               iplus2  = mod(limite+i+2,limite)
               
               if (imoins1.eq.0) imoins1 = ni
               if (i.eq.0) i = ni
               if (iplus1.eq.0) iplus1 = ni
               if (iplus2.eq.0) iplus2 = ni
               
               if (wrap.eq.1) then
                  if (iplus2.eq.ni) iplus2 = 2
                  if (imoins1.eq.ni) imoins1=ni-1
               endif
            else
               imoins1 = i-1
               iplus1  = i+1
               iplus2  = i+2
            endif
            dx = px(n) - i
            dy = py(n) - j
               
            index(n,1)=dx
            index(n,2)=dy
            index(n,3)=imoins1
            index(n,4)=i
            index(n,5)=iplus1
            index(n,6)=iplus2
            index(n,7)=j-1
            index(n,8)=j
            index(n,9)=j+1
            index(n,10)=j+2
         enddo
         !$OMP END PARALLEL DO 
         
         return
         end

         subroutine ez8_apply_3(index,zo,npts,z,ni,j1,j2,nodata)
            implicit none
            
            integer       npts,ni,j1,j2,i,j,n,iplus1,jplus1,imoins1,jmoins1,iplus2,jplus2
            real          zo(npts)
            real(kind=4)  index(npts,6)
            real          z(ni,j1:j2),nodata
            real(kind=8)  y1,y2,y3,y4
      
#include "cubic8.cdk"
            
            !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(n,i,j,imoins1,jmoins1,iplus1,jplus1,iplus2,jplus2,y1,y2,y3,y4,dx,dy) SHARED(npts,index,z,zo,nodata)
            do n=1,npts            
               dx     =index(n,1)
               dy     =index(n,2)
               imoins1=index(n,3)
               i      =index(n,3)
               iplus1 =index(n,5)
               iplus2 =index(n,6)
               jmoins1=index(n,7)
               j      =index(n,8)
               jplus1 =index(n,9)
               jplus2 =index(n,10)
      
               zo(n)=nodata

               if (defvalid(z(imoins1,jmoins1),nodata) .and. defvalid(z(i ,jmoins1),nodata) .and. defvalid(z(iplus1,jmoins1),nodata) .and. defvalid(z(iplus2,jmoins1),nodata)) then        
                  y1=cubic(dble(z(imoins1,jmoins1)),dble(z(i ,jmoins1)),dble(z(iplus1,jmoins1)),dble(z(iplus2,jmoins1)),dx)
                  if (defvalid(z(imoins1,j),nodata) .and. defvalid(z(i,j),nodata) .and. defvalid(z(iplus1,j),nodata) .and. defvalid(z(iplus2,j),nodata)) then        
                     y2=cubic(dble(z(imoins1,j)),dble(z(i ,j  )),dble(z(iplus1,j  )),dble(z(iplus2,j  )),dx)
                     if (defvalid(z(imoins1,jplus1),nodata) .and. defvalid(z(i ,jplus1),nodata) .and. defvalid(z(iplus1,jplus1),nodata) .and. defvalid(z(iplus2,jplus1),nodata)) then        
                        y3=cubic(dble(z(imoins1,jplus1)),dble(z(i ,jplus1)),dble(z(iplus1,jplus1)),dble(z(iplus2,jplus1)),dx)
                        if (defvalid(z(imoins1,jplus2),nodata) .and. defvalid(z(i ,jplus2),nodata) .and. defvalid(z(iplus1,jplus2),nodata) .and. defvalid(z(iplus2,jplus2),nodata)) then        
                           y4=cubic(dble(z(imoins1,jplus2)),dble(z(i ,jplus2)),dble(z(iplus1,jplus2)),dble(z(iplus2,jplus2)),dx)
                        endif
                     endif
                  endif
               endif
               
               zo(n)=cubic(y1,y2,y3,y4,dy)
            enddo
      
            return 
            end
      