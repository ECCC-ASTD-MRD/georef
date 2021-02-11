subroutine ez_uvfllc2d(u, v, x, y, x0, y0, x1, y1, x2, y2, x3, y3)
   implicit none
   real :: u, v, lat, lon, x, y

   real a,b,c,d,e, f, g,h, i
   real dx1, dx2, dy1, dy2, som_x ,som_y
   real x0, y0, x1, y1, x2, y2, x3, y3, q

   real, dimension(3, 3) ::  invmat

   som_x = x0 - x1 + x2 - x3
   som_y = y0 - y1 + y2 - y3

   dx1 = x1 - x2
   dx2 = x3 - x2
   dy1 = y1 - y2
   dy2 = y3 - y2

   g = (dy2*som_x-som_y*dx2)/(dx1*dy2-dx2*dy1)
   h = (dx1*som_y-som_x*dy1)/(dx1*dy2-dx2*dy1)
   i = 1.0

   a = x1 - x0 + g * x1
   b = x3 - x0 + h * x3
   c = x0
   d = y1 - y0 + g * y1
   e = y3 - y0 + h * y3
   f = y0

   invmat(1,3) = e * i - f * h
   invmat(2,3) = c * h - b * i
   invmat(3,3) = b * f - c * e
   invmat(1,2) = f * g - d * i
   invmat(2,2) = a * i - c * g
   invmat(3,2) = c * d - a * f
   invmat(1,1) = d * h - e * g
   invmat(2,1) = b * g - a * h
   invmat(3,1) = a * e - b * d

   q = invmat(1,1) * x + invmat(2,1) * y + invmat(3,1)
   if (q /= 0.0) then
      u = (invmat(1,3) * x + invmat(2,3) * y + invmat(3,3)) / q
      v = (invmat(1,2) * x + invmat(2,2) * y + invmat(3,2)) / q
      if (abs(u) < 0.01) u = abs(u)
      if (abs(v) < 0.01) v = abs(u)
   else
      u = -1.0
      v = -1.0
   endif

   end subroutine ez_uvfllc2d
   
   logical function pt_in_quad(x,y,x1,y1,x2,y2,x3,y3,x4,y4)
   implicit none

   real x,y,x1,y1,x2,y2,x3,y3,x4,y4
   logical pt_in_triangle
   external pt_in_triangle

   pt_in_quad = .false.
   if (pt_in_triangle(x,y,x1,y1,x2,y2,x3,y3)) then
      pt_in_quad = .true.
   else
      pt_in_quad = pt_in_triangle(x,y,x1,y1,x3,y3,x4,y4)
   endif

   return
   end function pt_in_quad

   logical function pt_in_triangle(x,y,x1,y1,x2,y2,x3,y3)
   implicit none

   real x,y,x1,y1,x2,y2,x3,y3
   real a,b,c,d, det, lambda1, lambda2, lambda3
   pt_in_triangle = .false.

   a = x1 - x3
   b = x2 - x3
   c = y1 - y3
   d = y2 - y3

   det = 1.0 / (a * d - b * c)

   lambda1 = (d * (x - x3) - b * (y - y3)) * det
   lambda2 = (a * (y - y3) - c * (x - x3)) * det
   lambda3 = 1.0 - lambda1 - lambda2

   if (lambda1 < 0.0 .or. lambda1 > 1.0) then
      return
   else if (lambda2 < 0.0 .or. lambda2 > 1.0) then
      return
   else if (lambda3 < 0.0 .or. lambda3 > 1.0) then
      return
   else
      pt_in_triangle =  .true.
   endif

   return
   end function pt_in_triangle


   subroutine inside_or_outside(masque,x, y, out_lat,out_lon, gdin_lat, gdin_lon, ni_src, nj_src, wgts, idxs, num_wgts)
   implicit none

   integer :: masque, ni_src, nj_src, num_wgts
   real out_lat, out_lon, gdin_lat(ni_src, nj_src), gdin_lon(ni_src, nj_src)
   real :: wgts(num_wgts), x, y
   integer :: idxs(num_wgts,2), i,j, ix
   logical pt_in_quad
   external pt_in_quad

   ix = minloc(wgts, 1)
   i = min(max(2,idxs(ix,1)),ni_src-1)
   j = min(max(2,idxs(ix,2)),nj_src-1)

   !  Cas 1 - Coin superieur droit

   if (pt_in_quad(out_lon, out_lat, gdin_lon(i-1,j-1), gdin_lat(i-1,j-1), &
         gdin_lon(i,j-1), gdin_lat(i,j-1), gdin_lon(i,j), gdin_lat(i,j), &
         gdin_lon(i-1,j), gdin_lat(i-1,j))) then
         masque = 1
         call ez_uvfllc2d(x, y,out_lon, out_lat, &
            gdin_lon(i-1,j-1), gdin_lat(i-1,j-1), &
            gdin_lon(i,j-1), gdin_lat(i,j-1), &
            gdin_lon(i,j), gdin_lat(i,j), &
            gdin_lon(i-1,j), gdin_lat(i-1,j))
         x = x + (i-1)
         y = y + (j-1)

   !  Cas 2 - Coin superieur gauche

   else if (pt_in_quad(out_lon, out_lat, gdin_lon(i,j-1), gdin_lat(i,j-1), &
         gdin_lon(i+1,j-1), gdin_lat(i+1,j-1), gdin_lon(i+1,j), gdin_lat(i+1,j), &
         gdin_lon(i,j), gdin_lat(i,j))) then
         masque = 1
         call ez_uvfllc2d(x, y,out_lon, out_lat, &
            gdin_lon(i,j-1), gdin_lat(i,j-1), &
            gdin_lon(i+1,j-1), gdin_lat(i+1,j-1), &
            gdin_lon(i+1,j), gdin_lat(i+1,j), &
            gdin_lon(i,j), gdin_lat(i,j))
         x = x + (i)
         y = y + (j-1)

!  Cas 3 - Coin inferieur droit

   else if (pt_in_quad(out_lon, out_lat, gdin_lon(i-1,j), gdin_lat(i-1,j), &
         gdin_lon(i,j), gdin_lat(i,j), gdin_lon(i,j+1), gdin_lat(i,j+1), &
         gdin_lon(i-1,j+1), gdin_lat(i-1,j+1))) then
         masque = 1
         call ez_uvfllc2d(x, y,out_lon, out_lat, &
            gdin_lon(i-1,j), gdin_lat(i-1,j), &
            gdin_lon(i,j), gdin_lat(i,j), &
            gdin_lon(i,j+1), gdin_lat(i,j+1), &
            gdin_lon(i-1,j+1), gdin_lat(i-1,j+1))
         x = x + (i-1)
         y = y + (j)

!  Cas 4 - Coin inferieur gauche

   else if (pt_in_quad(out_lon, out_lat, gdin_lon(i,j), gdin_lat(i,j), &
         gdin_lon(i+1,j), gdin_lat(i+1,j), gdin_lon(i+1,j+1), gdin_lat(i+1,j+1), &
         gdin_lon(i,j+1), gdin_lat(i,j+1))) then
         masque = 1
         call ez_uvfllc2d(x, y,out_lon, out_lat, &
            gdin_lon(i,j), gdin_lat(i,j), &
            gdin_lon(i+1,j), gdin_lat(i+1,j), &
            gdin_lon(i+1,j+1), gdin_lat(i+1,j+1), &
            gdin_lon(i,j+1), gdin_lat(i,j+1))
         x = x + (i)
         y = y + (j)
   else
      masque = 0
      x = -1.0
      y = -1.0
   endif

   return

   end subroutine inside_or_outside
