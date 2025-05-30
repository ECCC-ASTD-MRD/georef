      subroutine ez8_lambxyfll99(x,y,lat,lon,n,      latin1,latin2,yaxislat,yaxislon)
      implicit none
      integer n
      real(kind=8) x(n),y(n),lat(n),lon(n)
      real latin1,latin2,yaxislat,yaxislon
      
      real pi,pisur4,d2r,rn,rphi1,rphi2,r,f,rtan,rho,theta,rhozero,dlon,tmplat
      integer i
      
      pisur4 = atan(1.0)
      pi     = 4.0 * pisur4
      d2r    = pi / 180.0
      r = 6370997.0
      
      rphi1 = d2r * latin1
      rphi2 = d2r * latin2

      if (rphi1.eq.rphi2) then
         rn = sin(rphi1)
      else
         rn = log(cos(rphi1)/cos(rphi2))
         rn = rn / log((tan(pisur4+0.5*rphi2))/tan(pisur4+0.5*rphi1))
      endif
      
      rtan = tan(pisur4 + rphi1 * 0.5)
      
      f = (cos(rphi1)*(rtan ** rn))/rn
      rhozero = r * f / (tan(pisur4+yaxislat*d2r*.5)**rn)
      
      do i = 1,n
         tmplat = lat(i)
         if (tmplat.gt.90.0) tmplat = 89.95
         rho   = r * f / (tan(pisur4+tmplat*0.5*d2r)**rn)
         dlon = lon(i) - yaxislon

         if (dlon.lt.-180) then
            dlon = dlon + 360.
         else if (dlon.gt.180.) then
            dlon = dlon - 360.
         endif

         theta = rn * (d2r * dlon)
         
         x(i) = rho * sin(theta)
         y(i) = rhozero - rho*cos(theta)
      enddo
      
      return
      end

