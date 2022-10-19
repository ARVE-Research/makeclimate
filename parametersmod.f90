module parametersmod

use iso_fortran_env, only : int16,int32,real32,real64

implicit none

public  :: ndaymon
private :: leapyear

integer, parameter :: i2 = int16
integer, parameter :: i4 = int32
integer, parameter :: sp = real32
integer, parameter :: dp = real64

integer(i2), parameter :: imissing = -32768
real(sp), parameter    :: rmissing = -9999.

integer, parameter :: tlen_blk = 12 * 30  ! number of months for calculation and incremental output

integer :: xlen
integer :: ylen
integer :: tlen_out    ! number of months to output
integer :: tlen_anom   ! number of months in the anomaly input
integer :: bfid
integer :: ofid
integer :: numcyc

integer, allocatable, dimension(:) :: offset

contains

!----------------------------

function ndaymon(year)

implicit none

integer, parameter, dimension(12) :: nd0 = [ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 ]  !number of days in each month
integer, parameter, dimension(12) :: nd1 = [ 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 ]  !number of days in each month (leap year)

integer, intent(in) :: year  !requires year AD as input

integer, dimension(12) :: ndaymon

!---

if (leapyear(year)) then
  ndaymon = nd1
else
  ndaymon = nd0
end if

end function ndaymon

!----------------------------

logical function leapyear(year)

implicit none

integer, intent(in) :: year  !requires year AD as input

!---

if ((mod(year,4) == 0 .and. mod(year,100) /= 0) .or. mod(year,400) == 0) then
  leapyear = .true.
else
  leapyear = .false.
end if

end function leapyear

!----------------------------

end module parametersmod
