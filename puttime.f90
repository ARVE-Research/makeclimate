program puttime

! gfortran -o puttime calendarmod.f90 puttime.f90 -ftree-vectorize -Wall -I/home/public/easybuild/software/netCDF-Fortran/4.6.1-gompi-2023a/include -lnetcdff

use iso_fortran_env
use netcdf
use calendarmod,   only : timestruct,ymdt2jd

implicit none

integer, parameter :: sp = real32
integer, parameter :: dp = real64
integer, parameter :: i2 = int16

character(100) :: infile
character(15)  :: cyr
character(2)   :: bpflag = 'bp'

integer :: status
integer :: ncid
integer :: dimid
integer :: varid
integer :: tlen
integer :: t
integer :: y
integer :: m
integer :: y0
integer :: y1

integer :: yr
integer :: yrbp
integer :: bcad
integer :: nyrs

real(dp), allocatable, dimension(:) :: time

real(dp), dimension(2) :: time_range

type(timestruct) :: bt
type(timestruct) :: ts

! -------------------------------------------------------

call getarg(1,infile)
call getarg(2,cyr)
call getarg(3,bpflag)

! -------------------------------------------------------

status = nf90_open(infile,nf90_write,ncid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inq_dimid(ncid,'time',dimid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_inquire_dimension(ncid,dimid,len=tlen)
if (status /= nf90_noerr) call handle_err(status)

allocate(time(tlen))

nyrs = tlen / 12

! -------------------------------------------------------

read(cyr,*)yr

if (bpflag == 'bp') then

  ! the input year will equal the last year of the timeseries

  yrbp = yr

  if (yrbp >= 1950) then
    bcad = 1950 - yrbp - 1
  else
    bcad = 1950 - yrbp
  end if

  y1 = bcad
  y0 = y1 - nyrs + 1

else

 ! the input year is the first year of the timeseries

  bcad = yr
  
  y0 = bcad
  y1 = y0 + nyrs - 1

end if


write(0,*)'writing: ',bcad,nyrs,y0,y1

! -------------------------------------------------------

bt = timestruct(1950,1,1,0,0,0.)   ! reference all time to 1950-01-01

call ymdt2jd(bt)

t = 1

do y = y0,y1
  do m = 1,12
  
    ts = timestruct(y,m,1,0,0,0.)
    
    call ymdt2jd(ts)
    
    time(t) = ts%jd - bt%jd
    
    t = t + 1

  end do
end do

time_range = [minval(time),maxval(time)]

! -------------------------------------------------------

status = nf90_inq_varid(ncid,'time',varid)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_put_var(ncid,varid,time)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'units','days since 1950-01-01 00:00:00')
if (status /= nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'calendar','proleptic_gregorian')
if (status /= nf90_noerr) call handle_err(status)

status = nf90_put_att(ncid,varid,'actual_range',time_range)
if (status /= nf90_noerr) call handle_err(status)

status = nf90_close(ncid)
if (status /= nf90_noerr) call handle_err(status)

! -------------------------------------------------------

contains

subroutine handle_err(status)

! Internal subroutine - checks error status after each netcdf call,
! prints out text message each time an error code is returned. 

integer, intent (in) :: status

if(status /= nf90_noerr) then 
  print *, trim(nf90_strerror(status))
  stop
end if

end subroutine handle_err

! -------------------------------------------------------

end program puttime