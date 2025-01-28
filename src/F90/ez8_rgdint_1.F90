
      subroutine ez8_rgdint_1(zo,px,py,npts,z,ni,j1,j2,wrap,nodata)
         implicit none
         
         integer npts,ni,j1,j2,wrap,i,j,n,limite,iplus1
         real    zo(npts)
         real(kind=8)  px(npts),py(npts)
         real    z(ni,j1:j2),nodata
         real(kind=8)  dx,dy,y2,y3

#include "zlin8.cdk"
         
         limite = ni+1-wrap
         !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(n,i,j,iplus1,y2,y3,dx,dy) SHARED(npts,j1,j2,ni,px,py,wrap,z,zo,nodata,limite)
         do n=1,npts
            if (wrap.gt.0.) then
               i = min(ni-2+wrap,max(1,int(px(n))))
            else 
               i = min(ni-1+wrap,max(1,int(px(n))))
            endif
            j = min(j2-1,max(j1,int(py(n))))
            
            iplus1 = i + 1
            if (wrap.gt.0.and.(i.eq.(ni-2+wrap))) then
               iplus1  = mod(limite+i+1,limite)
            endif

            dx = px(n) - i
            dy = py(n) - j
            
            if (defvalid(z(i,j),nodata) .and.  defvalid(z(iplus1,j),nodata) .and. defvalid(z(i,j+1),nodata) .and. defvalid(z(iplus1,j+1),nodata)) then
               y2=zlin(dble(z(i,j  )),dble(z(iplus1,j  )),dx)
               y3=zlin(dble(z(i,j+1)),dble(z(iplus1,j+1)),dx)  
               zo(n)=zlin(y2,y3,dy)
            else
               zo(n)=nodata
            endif
         enddo

         return 
      end

      subroutine ez8_rgd_index_1(index,px,py,npts,ni,j1,j2,wrap)
         implicit none
         
         integer       npts,ni,j1,j2,wrap,i,j,n,limite,iplus1
         real(kind=4)  index(npts,6)
         real(kind=8)  px(npts),py(npts)
         real(kind=8)  dx,dy
            
         limite = ni+1-wrap
         !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(n,i,j,iplus1,dx,dy) SHARED(npts,j1,j2,ni,px,py,wrap,index,limite)
         do n=1,npts
            if (wrap.gt.0.) then
               i = min(ni-2+wrap,max(1,int(px(n))))
            else 
               i = min(ni-1+wrap,max(1,int(px(n))))
            endif
            j = min(j2-1,max(j1,int(py(n))))
            
            iplus1 = i + 1
            if (wrap.gt.0.and.(i.eq.(ni-2+wrap))) then
               iplus1  = mod(limite+i+1,limite)
            endif
   
            dx = px(n)-i
            dy = py(n)-j

            index(n,1)=dx
            index(n,2)=dy
            index(n,3)=i
            index(n,4)=j
            index(n,5)=iplus1
            index(n,6)=j+1
         enddo
   
      return 
      end

      subroutine ez8_apply_1(index,zo,npts,z,ni,j1,j2,nodata)
         implicit none
         
         integer       npts,ni,j1,j2,i,j,n,iplus1,jplus1
         real          zo(npts)
         real(kind=4)  index(npts,6)
         real          z(ni,j1:j2),nodata
         real(kind=8)  dx,dy,y2,y3
   
#include "zlin8.cdk"
         
         !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(n,i,j,iplus1,jplus1,y2,y3,dx,dy) SHARED(npts,index,z,zo,nodata)
         do n=1,npts            
            dx    =index(n,1)
            dy    =index(n,2)
            i     =index(n,3)
            j     =index(n,4)
            iplus1=index(n,5)
            jplus1=index(n,6)

            if (defvalid(z(i,j),nodata) .and.  defvalid(z(iplus1,j),nodata) .and. defvalid(z(i,jplus1),nodata) .and. defvalid(z(iplus1,jplus1),nodata)) then
               y2=zlin(dble(z(i,j     )),dble(z(iplus1,j     )),dx)
               y3=zlin(dble(z(i,jplus1)),dble(z(iplus1,jplus1)),dx)  
               zo(n)=zlin(y2,y3,dy)
            else
               zo(n)=nodata
            endif
         enddo
   
         return 
         end
   