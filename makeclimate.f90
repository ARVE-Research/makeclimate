program makeclimate

! ====================================================================================================
! Utility program to generate spinup climate data for input to LPJ-LMfire
! based on a climatological monthly mean baseline climate (12 months)
! and an interannual variability file containing monthly anomalies for
! several climate variables. Typically interannual variability is provided by the
! 20th Century Reanalysis in the form of anomalies reydimvarive to the period covered
! by the baseline climatology (varies per variable)

! ===== prerequisite input files both covering the target region and at target resolution =====
! --- baseline climate file ---
! elevation (elv) (m)
! temperature at 2m (tmp) (degC)
! diurnal temperature range at 2m (dtr) (degC)
! total accumuydimvared precipitation (pre) (mm)
! days with measurable precipitation (wet) (days)
! cloud cover full column (cld) (percent)
! windspeed at 10m (wnd) (m s-1)
! cloud-to-ground lightning stroke density (lght) (strokes km-2 d-1)

! --- interannual variability anomalies (separate files) ---
! temperature at 2m (tmp) (degC)
! diurnal temperature range (dtr) (degC)
! accumuydimvared precipitation (apcp) (mm)
! cloud cover full column (tdcd) (percent)
! windspeed at 10m (wspd) (m s-1)
! lightning strokes (lght) (strokes km-2 d-1)

! use ./configure and make to compile
! usage: ./makeclimate <jobfile> <outfile>
! jobfile is a fortran90 .namelist file containing the anomaly file root path and names

! major revision by JOK January 2021

! ====================================================================================================

use parametersmod, only : i2,sp,dp,ndaymon,imissing,rmissing,xlen,ylen,tlen_out,tlen_anom,bfid,ofid,numcyc,offset
use randomdistmod, only : randomstate,ran_seed,ranur
use netcdfmod,     only : ncstat,handle_err,netcdf_create
use calcvarmod,    only : calcvar
use calcwetmod,    only : calcwetf
use netcdf

implicit none

!------------------------------
! arguments

character(200) :: jobfile
character(200) :: outfile

! job options

character(200) :: basefile
character(200) :: anompath
character(200) :: tmpfile
character(200) :: dtrfile
character(200) :: prefile
character(200) :: cldfile
character(200) :: wndfile
character(200) :: wetVprefile
character(200) :: lghtfile

integer :: startyr =   0  !year (BP) of the first year of the output climatology
integer :: numyrs  = 1020  !default value in case it is not specified in the namelist

logical :: timeslice = .true.

namelist  / joboptions / basefile,anompath,tmpfile,dtrfile,prefile,cldfile,wndfile,wetVprefile,lghtfile,startyr,numyrs,timeslice

! local variables

! counters
integer :: a,b,c
integer :: i,j
integer :: yr
integer :: m

! netcdf id's
integer :: afid
integer :: dimid
integer :: varid

! lengths
integer, parameter :: ycyc = 30
integer :: anomyrs
integer :: nblocks

! other variables
integer :: seed
integer :: val

integer, allocatable, dimension(:) :: seg
integer, allocatable, dimension(:) :: offset30 ! index of the starting value of the 30-year climate blocks
logical, allocatable, dimension(:) :: used

integer, dimension(13) :: nd
integer, dimension(12) :: ndm

type(randomstate) :: rndst

real(dp),    allocatable, dimension(:)   :: xdimvar
real(dp),    allocatable, dimension(:)   :: ydimvar
real(dp),    allocatable, dimension(:)   :: time
integer(i2), allocatable, dimension(:,:) :: var2d

character(10) :: xdimname = 'lon'
character(10) :: ydimname = 'lat'



! integer :: x,y,k

! integer :: day

! integer :: climmonths

! character(80) :: status_msg

! integer(i2) :: missing

! real(sp),   dimension(2) :: actual_range

! integer(i2), allocatable, dimension(:,:,:) :: base
! integer(i2), allocatable, dimension(:,:,:) :: vout

! real(sp),    allocatable, dimension(:,:,:) :: anom
! real(sp),    allocatable, dimension(:,:,:) :: rbase
! real(sp),    allocatable, dimension(:,:,:) :: rvout

! real(sp) :: scale_factor
! real(sp) :: add_offset

!-----------------------------------------------------------------------------------------------------
! program starts here
! read job options

call getarg(1,jobfile)
call getarg(2,outfile)

open(10,file=jobfile,status='old')
read(10,nml=joboptions)
close(10)

!------------------------------------------------------------------------

! inquire about the length of the anomaly timeseries and create the offset vector

ncstat = nf90_open(trim(anompath)//tmpfile,nf90_nowrite,afid)
if (ncstat /= nf90_noerr) call handle_err(ncstat)

ncstat = nf90_inq_dimid(afid,'time',dimid)
if (ncstat /= nf90_noerr) call handle_err(ncstat)

ncstat = nf90_inquire_dimension(afid,dimid,len=tlen_anom)
if (ncstat /= nf90_noerr) call handle_err(ncstat)

ncstat = nf90_close(afid)
if (ncstat /= nf90_noerr) call handle_err(ncstat)

anomyrs = tlen_anom / 12

! work out the 30-year offsets

a = anomyrs / 30
b = (anomyrs - 10) / 30
c = (anomyrs - 20) / 30

nblocks = a + b + c

write(0,*)'anomalies have',anomyrs,'years of data'
write(0,*)tlen_anom
write(0,*)a,b,c,nblocks

allocate(offset30(nblocks))
allocate(used(nblocks))

offset30 = 0

m = 1  ! start at month one
do i = 1,a
  offset30(i) = m
  m = m + 360
end do

m = 121  ! 10 years offset
do i = a+1,a+b
  offset30(i) = m
  m = m + 360
end do

m = 241 ! 20 years offset
do i = a+b+1,nblocks
  offset30(i) = m
  m = m + 360
end do

!------------------------------------------------------------------------
! calcuydimvare the random sequence of blocks for timeslice (spinup) output

if (timeslice) then

  seed = -(32768 + startyr)

  call ran_seed(seed,rndst)   !initialize the random seed with the start year

  !------
  !select periods for construction of climatology

  numcyc = numyrs / 30  ! this will truncate downwards if the number of years requested is not an even multiple

  if (mod(numyrs,30) /= 0) then
    write(0,*)'Warning: for spinup climatologies an even multiple of 30 years should be selected. Adjusting...'
    numcyc = numcyc + 1
  end if

  write(0,'(i5,a)')numyrs,' years of climate requested'
  write(0,'(a,2i5,a,i5,a)')' will generate',numcyc,ycyc,'-year climate cycles = ',numcyc*ycyc,' years of climate'
  write(0,*)'generating sequence order'

  !generate a pseudo-random sequence of 30 year climate blocks
  !avoid repeating some blocks frequently and never using others by setting a flag vector
  !to keep choosing random numbers until one that has not already been used is selected
  !once all values have been used once, reset the flag vector and start again

  open(20,file='sequence.dat',status='unknown')

  allocate(seg(numcyc))
  allocate(offset(numcyc))

  used = .false.

  do i = 1,numcyc

    do
      val = nint(1. + real(nblocks-1) * ranur(rndst))  !generates a random integer between 1 and nblocks

      if (.not.used(val)) exit
    end do

    seg(i) = val
    used(val) = .true.

    if (all(used)) used = .false.

    write(20,*)i,seg(i)

    offset(i) = offset30(seg(i))

    ! write(0,*)i,seg(i),offset(i)

  end do

  close(20)

else

  !transient climatology; use the transient interannual variability from the input anomalies

!   ycyc = 20
!   tlen = 12 * ycyc
! 
!   numcyc = 7
! 
!   allocate(seg(numcyc))
!   allocate(offset(numcyc))
! 
!   do i = 1,numcyc
!     seg(i) = i
!     offset(i) = 1 + tlen * (i - 1)
! 
!     write(0,*)seg(i),offset(i)
! 
!   end do
! 
!   write(0,'(a,2i5,a,i5,a)')' transient climatology:',numcyc,ycyc,'-year climate cycles = ',numcyc*ycyc,' years of climate'

end if

!------------------------------------------------------------------------
! get the dimensions of the input array

ncstat = nf90_open(basefile,nf90_nowrite,bfid)
if (ncstat /= nf90_noerr) call handle_err(ncstat)

ncstat = nf90_inq_dimid(bfid,xdimname,dimid)
if (ncstat /= nf90_noerr) then  ! try x instead
  xdimname = 'x'
  ncstat = nf90_inq_dimid(bfid,xdimname,dimid)
end if
if (ncstat /= nf90_noerr) call handle_err(ncstat)

ncstat = nf90_inquire_dimension(bfid,dimid,len=xlen)
if (ncstat /= nf90_noerr) call handle_err(ncstat)

ncstat = nf90_inq_dimid(bfid,ydimname,dimid)
if (ncstat /= nf90_noerr) then  ! try y instead
  ydimname = 'y'
  ncstat = nf90_inq_dimid(bfid,ydimname,dimid)
end if
if (ncstat /= nf90_noerr) call handle_err(ncstat)

ncstat = nf90_inquire_dimension(bfid,dimid,len=ylen)
if (ncstat /= nf90_noerr) call handle_err(ncstat)

tlen_out = numyrs * 12

!------------------------------------------------------------------------
! generate the output file
! easier to generate the output file using an ncdump of the baseline input file

! call netcdf_create(outfile,xdimname,ydimname,xlen,ylen,tlen_out,ofid)

ncstat = nf90_open(outfile,nf90_write,ofid)
if (ncstat /= nf90_noerr) call handle_err(ncstat)

!------------------------------------------------------------------------
! write the dimension variables

! this will all be done through copy of the baseline input file

! write(0,*)'write dimension variables'
! write(0,*)'NB, dimension variable names'
! write(0,*)'x ',xdimname
! write(0,*)'y ',ydimname
! 
! ! ---
! ! xdimvar
! 
! allocate(xdimvar(xlen))
! 
! ncstat = nf90_inq_varid(bfid,xdimname,varid)
! if (ncstat /= nf90_noerr) call handle_err(ncstat)
! 
! ncstat = nf90_get_var(bfid,varid,xdimvar)
! if (ncstat /= nf90_noerr) call handle_err(ncstat)
! 
! ncstat = nf90_inq_varid(ofid,xdimname,varid)
! if (ncstat /= nf90_noerr) call handle_err(ncstat)
! 
! ncstat = nf90_put_var(ofid,varid,xdimvar)
! if (ncstat /= nf90_noerr) call handle_err(ncstat)
! 
! deallocate(xdimvar)
! 
! ! ---
! ! ydimvar
! 
! allocate(ydimvar(ylen))
! 
! ncstat = nf90_inq_varid(bfid,ydimname,varid)
! if (ncstat /= nf90_noerr) call handle_err(ncstat)
! 
! ncstat = nf90_get_var(bfid,varid,ydimvar)
! if (ncstat /= nf90_noerr) call handle_err(ncstat)
! 
! ncstat = nf90_inq_varid(ofid,ydimname,varid)
! if (ncstat /= nf90_noerr) call handle_err(ncstat)
! 
! ncstat = nf90_put_var(ofid,varid,ydimvar)
! if (ncstat /= nf90_noerr) call handle_err(ncstat)
! 
! deallocate(ydimvar)

! ---
! time

allocate(time(tlen_out))

! --

yr = 1

ndm = ndaymon(yr)  ! initialize the days per month vector

nd = 0

do j = 1,12
  nd = eoshift(nd,1,ndm(j))
end do

time(1) = nd(1)

do i = 2,tlen_out
  m = mod(i,12)
  
  if (m == 0) m = 12

  if (m == 1) then  ! refill the days per month vector
  
    yr = yr + 1

    ndm = ndaymon(yr)

    do j = 1,12
      nd = eoshift(nd,1,ndm(j))
    end do
    
  end if

  time(i) = time(i-1) + nd(m)
  
end do 

! --

write(0,*)'puttime',tlen_out

ncstat = nf90_inq_varid(ofid,'time',varid)
if (ncstat /= nf90_noerr) call handle_err(ncstat)

ncstat = nf90_put_var(ofid,varid,time)
if (ncstat /= nf90_noerr) call handle_err(ncstat)

!------------------------------------------------------------------------
! write the regular variables

!--------------------------------------------------
! elevation

allocate(var2d(xlen,ylen))

ncstat = nf90_inq_varid(bfid,'elv',varid)
if (ncstat /= nf90_noerr) call handle_err(ncstat)

ncstat = nf90_get_var(bfid,varid,var2d)
if (ncstat /= nf90_noerr) call handle_err(ncstat)

ncstat = nf90_inq_varid(ofid,'elv',varid)
if (ncstat /= nf90_noerr) call handle_err(ncstat)

ncstat = nf90_put_var(ofid,varid,var2d)
if (ncstat /= nf90_noerr) call handle_err(ncstat)

deallocate(var2d)

!--------------------------------------------------
! temperature

call calcvar('tmp','air',trim(anompath)//tmpfile)

!--------------------------------------------------
! dtr

call calcvar('dtr','dtr',trim(anompath)//dtrfile,llimit=0.)

!--------------------------------------------------
! precipitation

call calcvar('pre','apcp',trim(anompath)//prefile,llimit=0.)

!--------------------------------------------------
! cloud

call calcvar('cld','tcdc',trim(anompath)//cldfile,llimit=0.,ulimit=100.)

!--------------------------------------------------
! wind

call calcvar('wnd','wspd',trim(anompath)//wndfile,llimit=0.)

!--------------------------------------------------
! lightning

call calcvar('lght','lght',trim(anompath)//lghtfile,llimit=0.)

!--------------------------------------------------
! wet days

call calcwetf(wetVprefile)

!------------------------------------------------------------------------
! close files 

100 continue

write(0,*)
write(0,*)'clean up'

ncstat = nf90_close(bfid)
if (ncstat /= nf90_noerr) call handle_err(ncstat)

ncstat = nf90_close(ofid)
if (ncstat /= nf90_noerr) call handle_err(ncstat)

!------------------------------------------------------------------------------------------------------------------------

end program makeclimate
