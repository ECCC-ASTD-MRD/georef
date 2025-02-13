subroutine ez8_testprojt
    implicit none

    integer i,j,ni,nj
    real clat, clon, d60, dgrw

    real(kind=8) tlat(5,5), tlon(5,5)
    real(kind=8) tlat2(5,5), tlon2(5,5)
    real(kind=8) x(5,5),y(5,5)

    ni = 11
    nj = 11
    clat = 45.0
    clon = 260.0
    d60 = 100000.0
    dgrw = 0.0

    do j=1,5
        do i=1,5
            tlat(i,j) = 35.0 + (j-1)*5
            tlon(i,j) = 250.0 +(i-1)*5
        enddo
    enddo

    call ez8_vtxyfll(x, y, tlat, tlon, clat, clon, d60, dgrw, ni, nj, 25)

    do j=1,5
        do i=1,5
            print *, i,j,tlat(i,j),tlon(i,j),x(i,j),y(i,j)
        enddo
    enddo

    print *, '***************************************************'

    call ez8_vtllfxy(tlat2, tlon2, x, y, clat, clon, d60, dgrw, ni, nj, 25)

    do j=1,5
        do i=1,5
            print *, i,j,tlat2(i,j)-tlat(i,j),tlon2(i,j) - tlon(i,j)
        enddo
    enddo

    stop
end


subroutine ez8_vtxyfll(x, y, lat, lon, clat, clon, d60, dgrw, ni, nj, n)
    implicit none

    #include "pi.inc"

    integer i, n, ni, nj
    real(kind=8) x(n), y(n), lat(n), lon(n)
    real clat, clon, d60, dgrw
    real r
    real(kind=8) k
    real offsetx, offsety,sinclat,cosclat

    r = 6371000.0
    sinclat = sin (clat * dgtord)
    cosclat = cos (clat * dgtord)

    offsetx = (-ni-1) * 0.5
    offsety = (-nj-1) * 0.5

    !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,k) SHARED(n,x,y,r,dgtord,sinclat,cosclat,lat,lon,clon,offsetx,offsety,d60)
    do i=1,n
        k = 2.0 / (1.0 + sinclat*sin(lat(i) * dgtord)+ cosclat* cos(lat(i)*dgtord)*cos(dgtord*(lon(i)-clon)))
        x(i) = r * k * cos(lat(i)*dgtord) * sin(dgtord*(lon(i)-clon))
        y(i) = r * k * (cosclat*sin(lat(i) * dgtord) -  (sinclat * cos(lat(i)*dgtord)*cos(dgtord*(lon(i)-clon))))
        x(i) = x(i) / d60 - offsetx
        y(i) = y(i) / d60 - offsety
    enddo
end


subroutine ez8_vtllfxy(lat, lon, x, y, clat, clon, d60, dgrw, ni, nj, n)
      implicit none

#include "pi.inc"

    integer i, n, ni, nj
    real(kind=8) x(n), y(n), lat(n), lon(n)
    real clat, clon, d60, dgrw
    real r
    real offsetx, offsety,sinclat,cosclat
    real(kind=8) rho, c, a, b, temp


    r = 6371000.0
    sinclat = sin (clat * dgtord)
    cosclat = cos (clat * dgtord)

    offsetx = (-ni-1) * 0.5
    offsety = (-nj-1) * 0.5

    !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,rho,temp,a,b) SHARED(n,x,y,r,rdtodg,c,sinclat,cosclat,lat,lon,clat,clon,offsetx,offsety,d60)
    do i=1,n
        x(i) = (x(i) + offsetx) * d60
        y(i) = (y(i) + offsety) * d60

        rho = sqrt(x(i)*x(i) + y(i)*y(i))
        if (rho.eq.0.0) then
            lat(i) = clat
            lon(i) = clon
        else
            c = 2.0 * atan(rho/(2.0*r))

            temp = cos(c)*sinclat + (y(i) * sin(c) * cosclat)/rho
            lat(i) = rdtodg * asin(max(-1.0D0, min(temp, 1.0D0)))
            lon(i) = clon + rdtodg * atan(((x(i)*sin(c))/            (rho*cosclat*cos(c)-y(i)*sinclat*sin(c))))
            a = x(i) * sin(c)
            b = rho*cosclat*cos(c)-y(i)*sinclat*sin(c)
            lon(i) = clon + rdtodg * atan2(a,b)


        endif
        lon(i) = mod(mod(lon(i),360.0d0)+360.0,360.0d0)
    enddo
end
