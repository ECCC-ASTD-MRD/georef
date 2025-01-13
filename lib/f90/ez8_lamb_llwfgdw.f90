      subroutine ez8_lamb_llwfgdw(z1,z2,xlon,li,lj,grtyp,ig1,ig2,ig3,ig4,x1,y1,xlat)
      implicit none

      integer li,lj
      character grtyp
      integer ig1,ig2,ig3,ig4   
      real z1(li*lj), z2(li*lj)
      real(kind=8) xlon(li*lj),xlat(li*lj,2)
      real(kind=8) x1(li*lj,2),y1(li*lj,2)

      real delx,dely,uuu,vvv,alpha
      integer i

!     
!     * 1.866025=(1+sin60),   6.371e+6=earth radius in meters.      
!     
!     rdtodg = 180/pie, dgtord = pie/180        
      real spd0,dir0
      real pie,rdtodg,dgtord
      
      data pie    /3.1415926535898/     
      data rdtodg /57.295779513082/
      data dgtord /1.7453292519943e-2/

      do i=1,li*lj
         xlat(i,1) = 45.0
         xlat(i,2) = 50.0
      enddo

      call ez8_lambfll(x1(1,1),y1(1,1),xlat(1,1),xlon(1),li*lj,grtyp,ig1,ig2,ig3,ig4)
      call ez8_lambfll(x1(1,2),y1(1,2),xlat(1,2),xlon(1),li*lj,grtyp,ig1,ig2,ig3,ig4)

      do i=1,li*lj
         delx = x1(i,2) - x1(i,1)
         dely = y1(i,2) - y1(i,1)
         alpha = pie*0.50 - atan2(dely,delx)
         uuu = z1(i)
         vvv = z2(i)
         z1(i) = uuu*cos(alpha)-vvv*sin(alpha)
         z2(i) = uuu*sin(alpha)+vvv*cos(alpha)
         
         spd0=sqrt(z1(i)*z1(i)+z2(i)*z2(i))
         if( (spd0.eq. 0.0))then
            dir0= 0.0
         else 
            if( (z1(i).eq. 0.0))then
               if( (z2(i).ge. 0.0))then
                  dir0= 180.0
               else 
                  dir0= 0.0
               endif 
            else 
               dir0=270.0 - rdtodg*atan2(z2(i),z1(i))
            endif 
         endif 
         dir0=amod(amod(dir0,360.0)+360.0,360.0)
         z1(i)=spd0
         z2(i)=dir0
      enddo
      
      return
      end
