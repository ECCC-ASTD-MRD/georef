
      subroutine ez8_irgdint_1(zo,px,py,npts,z,ni,j1,j2,ax,ay,wrap,nodata)
         implicit none
         
         integer       npts,ni,wrap,limite
         integer       i,j,j1,j2,n,iplus1
         real          zo(npts)
         real(kind=8)  px(npts),py(npts)
         real(kind=8)  ax(ni),ay(j1:j2)
         real          z(ni,j1:j2),nodata
         real(kind=8)  x,y,x1,x2,y1,y2,dx,dy

#include "zlin8.cdk"

         limite = ni+2-wrap
      
         !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(n,i,j,iplus1,x,x1,x2,y,y1,y2,dx,dy) SHARED(npts,j1,j2,ni,px,py,wrap,ax,ay,z,zo,nodata,limite)
         do n=1,npts
            if (wrap.gt.0) then
               i = min(ni-1,max(1,int(px(n))))
            else 
               i = min(ni-2+wrap,max(1,int(px(n))))
            endif
            j = min(j2-1,max(j1+1, int(py(n))))
            
            if (j.lt.0) then
               j = j-1
            endif
            
            iplus1 = i+1
            
            x1=ax(i)
            x2=ax(iplus1)

            if (wrap.gt.0.and.(i.eq.(ni-2+wrap))) then
               iplus1 = mod(limite+i+1,limite)
               x2=ax(2)+ax(ni)
            endif
            
            x = x1 + (x2-x1)*(px(n)-i)
            y = ay(j) + (ay(j+1)-ay(j))*(py(n)-j)
            
            dx = (x - x1)/(x2-x1)
            dy = (y - ay(j))/(ay(j+1)-ay(j))
            
            if (defvalid(z(i,j),nodata) .and. defvalid(z(iplus1,j),nodata) .and. defvalid(z(i,j+1),nodata) &
               .and. defvalid(z(iplus1,j+1),nodata)) then
               y1 = zlin(dble(z(i,j)),dble(z(iplus1,j)),dx)
               y2 = zlin(dble(z(i,j+1)),dble(z(iplus1,j+1)),dx)
               zo(n) = zlin(y1,y2,dy)
            else
               zo(n)=nodata
            endif
         enddo
      
      return
      end
      
      subroutine ez8_irgd_index_1(index,px,py,npts,ni,j1,j2,ax,ay,wrap)
         implicit none
         
         integer       npts,ni,wrap,limite
         integer       i,j,j1,j2,n,iplus1
         real(kind=4)  index(npts,6)
         real          zo(npts)
         real(kind=8)  px(npts),py(npts)
         real(kind=8)  ax(ni),ay(j1:j2)
         real          z(ni,j1:j2),nodata
         real(kind=8)  x,y,x1,x2,y1,y2,dx,dy

         limite = ni+2-wrap
      
         !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(n,i,j,iplus1,x,x1,x2,y,y1,y2,dx,dy) SHARED(npts,j1,j2,ni,px,py,wrap,ax,ay,index,limite)
         do n=1,npts
            if (wrap.gt.0) then
               i = min(ni-1,max(1,int(px(n))))
            else 
               i = min(ni-2+wrap,max(1,int(px(n))))
            endif
            j = min(j2-1,max(j1+1, int(py(n))))
            
            if (j.lt.0) then
               j = j-1
            endif
            
            iplus1 = i+1
            
            x1=ax(i)
            x2=ax(iplus1)

            if (wrap.gt.0.and.(i.eq.(ni-2+wrap))) then
               iplus1 = mod(limite+i+1,limite)
               x2=ax(2)+ax(ni)
            endif
            
            x = x1 + (x2-x1)*(px(n)-i)
            y = ay(j) + (ay(j+1)-ay(j))*(py(n)-j)
            
            dx = (x - x1)/(x2-x1)
            dy = (y - ay(j))/(ay(j+1)-ay(j))
            
            index(n,1)=dx
            index(n,2)=dy
            index(n,3)=i
            index(n,4)=j
            index(n,5)=iplus1
            index(n,6)=j+1
         enddo
   
      return 
      end
