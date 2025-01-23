
      subroutine ez8_rgdint_0(zo,px,py,npts,z,ni,j1,j2,nodata)
         implicit none
         
         integer npts,ni,j1,j2,i,j,n
         real    zo(npts)
         real(kind=8)  px(npts),py(npts)
         real    z(ni,j1:j2),nodata
         
         !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(n,i,j) SHARED(npts,j1,j2,ni,px,py,z,zo)
         do n=1,npts
            i = min(ni,max(1,nint(px(n))))
            j = min(j2,max(j1,nint(py(n))))
            
            zo(n)=z(i,j)
         enddo
         
         return
      end
      
      subroutine ez8_rgd_index_0(index,px,py,npts,ni,j1,j2)
         implicit none
         
         integer       npts,ni,j1,j2,i,j,n
         real(kind=4)  index(npts,2)
         real(kind=8)  px(npts),py(npts)
         
         !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(n,i,j) SHARED(npts,index,ni,j1,j2,px,py)
         do n=1,npts
            i = min(ni,max(1,nint(px(n))))
            j = min(j2,max(j1,nint(py(n))))
            
            index(n,1)=i
            index(n,2)=j
         enddo

      return 
      end

      subroutine ez8_apply_0(index,zo,npts,z,ni,j1,j2,nodata)
         implicit none
         
         integer       npts,ni,j1,j2,i,j,n
         real(kind=4)  index(npts,2)
         real          zo(npts),nodata
         real(kind=8)  px(npts),py(npts)
         real          z(ni,j1:j2)
         
         !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(n,i,j) SHARED(npts,index,z,zo)
         do n=1,npts
            i=index(n,1)
            j=index(n,2)
            
            zo(n)=z(i,j)
         enddo
  
         return 
         end
   