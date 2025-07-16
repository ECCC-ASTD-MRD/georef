module GeoRef_internal_def
    implicit none

    real, parameter :: pie = 3.1415926535898
    ! rdtodg = 180/pie
    real, parameter :: rdtodg = 57.295779513082
    ! dgtord = pie/180
    real, parameter :: dgtord = 1.7453292519943e-2

    integer, parameter :: global = 0
    integer, parameter :: nord = 1
    integer, parameter :: sud = 2
    integer, parameter :: sudnord = 0
    integer, parameter :: nordsud = 1

    integer, parameter :: absolu = 0
    integer, parameter :: relatif = 1
end module