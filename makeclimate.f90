program makeclimate

! ====================================================================================================
! Utility program to generate transient climate data for input to LPJ-LMfire
! based on a climatological monthly mean baseline climate (12 months)
! and an interannual variability file containing monthly anomalies for
! several climate variables. Typically interannual variability is provided by the
! 20th Century Reanalysis in the form of anomalies relative to the period covered
! by the baseline climatology (varies per variable)

! ===== prerequisite input files both covering the target region and at target resolution =====
! --- baseline climate file ---
! elevation (elv) (m)
! temperature at 2m (tmp) (degC)
! diurnal temperature range at 2m (dtr) (degC)
! total accumulated precipitation (pre) (mm)
! days with measurable precipitation (wet) (days)
! cloud cover full column (cld) (percent)
! windspeed at 10m (wnd) (m s-1)
! cloud-to-ground lightning stroke density (lght) (strokes km-2 d-1)

! --- interannual variability anomalies (separate files) ---
! temperature at 2m (tmp) (degC)
! diurnal temperature range (dtr) (degC)
! accumulated precipitation (apcp) (mm)
! cloud cover full column (tdcd) (percent)
! windspeed at 10m (wspd) (m s-1)
! lightning strokes (lght) (strokes km-2 d-1)

! use ./configure and make to compile
! usage: ./makeclimate <jobfile> <outfile>
! jobfile is a fortran90 .namelist file containing the anomaly file root path and names

! major revision by JOK January 2021
! further revisions to support both transient and spinup climate, January 2024

! ====================================================================================================

use parametersmod, only : i2,sp,dp,ndaymon,imissing,rmissing,xlen,ylen,tlen_out,tlen_anom,bfid,ofid,numcyc,offset
use randomdistmod, only : randomstate,ran_seed,ranur
use netcdfmod,     only : ncstat,handle_err,netcdf_create
use calcvarmod,    only : calcvar
use calcwetmod,    only : calcwetf
use netcdf

implicit none

! ------------------------------
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

! default values for job options in case they are not specified in the namelist

integer :: baseyr    =   1900   ! first year of the transient climate or year after the final year of spinup (standard calendar year)
integer :: numyrs    =   1020   ! number of spinup years to calculate if not specified in the namelist
logical :: transient = .false.  ! generate spinup or transient climate

namelist  / joboptions / basefile,anompath,tmpfile,dtrfile,prefile,cldfile,wndfile,wetVprefile,lghtfile,baseyr,numyrs,transient

! local variables

! counters
integer :: a,b,c
integer :: i,j
integer :: yr
integer :: m
integer :: t
integer :: t0
integer :: t1

! netcdf id's
integer :: afid
integer :: dimid
integer :: varid

! lengths
integer, parameter :: ycyc = 30  ! time block size parameter in years - data will be written out in blocks of this length
integer :: anomyrs
integer :: nblocks
integer :: tlen

! other variables
integer :: seed
integer :: val

integer, allocatable, dimension(:) :: seg
integer, allocatable, dimension(:) :: offset30 ! index of the starting value of the 30-year climate blocks
logical, allocatable, dimension(:) :: used

integer, dimension(13) :: nd
integer, dimension(12) :: ndm
integer, dimension(12) :: cumdays

type(randomstate) :: rndst

real(dp),    allocatable, dimension(:)   :: xdimvar
real(dp),    allocatable, dimension(:)   :: ydimvar
real(dp),    allocatable, dimension(:)   :: time
integer(i2), allocatable, dimension(:,:) :: var2d

character(10) :: xdimname = 'lon'
character(10) :: ydimname = 'lat'

! -----------------------------------------------------------------------------------------------------
! program starts here
! read job options

call getarg(1,jobfile)
call getarg(2,outfile)

open(10,file=jobfile,status='old')
read(10,nml=joboptions)
close(10)

! ------------------------------------------------------------------------
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

write(0,'(a,i0,a)')'anomalies have ',anomyrs,'years of data'

! ------------------------------------------------------------------------

if (transient) then

  ! -----------------------------    
  ! transient climate: use the transient interannual variability from the input anomalies
  
  tlen = 12 * ycyc   ! ycyc time block size parameter (set to 30 yrs)

  numcyc = tlen_anom / tlen  ! numcyc number of cycles that will be calculated

  allocate(seg(numcyc))
  allocate(offset(numcyc))

  do i = 1,numcyc
  
    seg(i) = i
    offset(i) = 1 + tlen * (i - 1)

    write(0,*)seg(i),offset(i)

  end do

  write(0,'(a,2i0,a,i0,a)')' transient climatology:',numcyc,ycyc,'-year climate cycles = ',numcyc*ycyc,' years of climate'

else

  ! -----------------------------    
  ! spinup climate: calculate the random sequence of blocks for repeating output
  
  ! work out the 30-year offsets

  a = anomyrs / 30
  b = (anomyrs - 10) / 30
  c = (anomyrs - 20) / 30

  nblocks = a + b + c
  
  write(0,*)a,b,c,nblocks
  
  allocate(offset30(nblocks))
  allocate(used(nblocks))

  offset30 = 0

  m = 1  ! start at year one month one
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

  ! -----------------------------    

  seed = -(32768 + baseyr)

  call ran_seed(seed,rndst)   !initialize the random seed with the start year

  ! ------
  ! select periods for construction of climatology

  numcyc = numyrs / 30  ! this will truncate downwards if the number of years requested is not an even multiple

  if (mod(numyrs,30) /= 0) then
    write(0,*)'Warning: for spinup climatologies an even multiple of 30 years should be selected. Adjusting...'
    numcyc = numcyc + 1
  end if

  write(0,'(i0,a)')numyrs,' years of climate requested'
  write(0,'(a,i0,a,i0,a,i0,a)')' will generate ',numcyc," ",ycyc,'-year climate cycles = ',numcyc*ycyc,' years of climate'
  write(0,*)'generating sequence order'

  ! generate a pseudo-random sequence of 30 year climate blocks
  ! avoid repeating some blocks frequently and never using others by setting a flag vector
  ! to keep choosing random numbers until one that has not already been used is selected
  ! once all values have been used once, reset the flag vector and start again

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
  
  ! reset the base year for the first year of the spinup
  
  baseyr = baseyr - numyrs

end if

! ------------------------------------------------------------------------
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

! ------------------------------------------------------------------------
! the output file must be generated using an ncdump of the baseline input file 
! in the wrapper shell script

ncstat = nf90_open(outfile,nf90_write,ofid)
if (ncstat /= nf90_noerr) call handle_err(ncstat)

! ------------------------------------------------------------------------
! write the dimension variables

allocate(xdimvar(xlen))
allocate(ydimvar(ylen))
allocate(var2d(xlen,ylen))

! ----------------------
! x

ncstat = nf90_inq_varid(bfid,'x',varid)
if (ncstat /= nf90_noerr) call handle_err(ncstat)

ncstat = nf90_get_var(bfid,varid,xdimvar)
if (ncstat /= nf90_noerr) call handle_err(ncstat)

ncstat = nf90_inq_varid(ofid,'x',varid)
if (ncstat /= nf90_noerr) call handle_err(ncstat)

ncstat = nf90_put_var(ofid,varid,xdimvar)
if (ncstat /= nf90_noerr) call handle_err(ncstat)

! ----------------------
! y

ncstat = nf90_inq_varid(bfid,'y',varid)
if (ncstat /= nf90_noerr) call handle_err(ncstat)

ncstat = nf90_get_var(bfid,varid,ydimvar)
if (ncstat /= nf90_noerr) call handle_err(ncstat)

ncstat = nf90_inq_varid(ofid,'y',varid)
if (ncstat /= nf90_noerr) call handle_err(ncstat)

ncstat = nf90_put_var(ofid,varid,ydimvar)
if (ncstat /= nf90_noerr) call handle_err(ncstat)

! ----------------------
! lon

ncstat = nf90_inq_varid(bfid,'lon',varid)
if (ncstat /= nf90_noerr) call handle_err(ncstat)

ncstat = nf90_get_var(bfid,varid,var2d)
if (ncstat /= nf90_noerr) call handle_err(ncstat)

ncstat = nf90_inq_varid(ofid,'lon',varid)
if (ncstat /= nf90_noerr) call handle_err(ncstat)

ncstat = nf90_put_var(ofid,varid,var2d)
if (ncstat /= nf90_noerr) call handle_err(ncstat)

! ----------------------
! lat

ncstat = nf90_inq_varid(bfid,'lat',varid)
if (ncstat /= nf90_noerr) call handle_err(ncstat)

ncstat = nf90_get_var(bfid,varid,var2d)
if (ncstat /= nf90_noerr) call handle_err(ncstat)

ncstat = nf90_inq_varid(ofid,'lat',varid)
if (ncstat /= nf90_noerr) call handle_err(ncstat)

ncstat = nf90_put_var(ofid,varid,var2d)
if (ncstat /= nf90_noerr) call handle_err(ncstat)

! ----------------------
! time

allocate(time(tlen_out))

t0 = 0

write(0,*)'base year',baseyr

yr = baseyr

do t = 1,tlen_out,12

  ! write(0,*)yr,ndaymon(yr)

  ndm = eoshift(ndaymon(yr),-1,0)
  
  cumdays = [(sum(ndm(1:m)),m=1,12)]

  time(t:t+11) = cumdays + t0
  
  t0 = time(t+11) + 31.

  yr = yr + 1
  
end do

if (.not.transient) then

  ! invert the time array so it counts down towards the last day of the spinup
  
  t1 = time(tlen_out) + 31  ! set this to 1 January of the first year after the end of the spinup

  do i = 1,tlen_out
    time(i) = time(i) - t1
  end do

end if

write(0,*)'writing time',tlen_out

ncstat = nf90_inq_varid(ofid,'time',varid)
if (ncstat /= nf90_noerr) call handle_err(ncstat)

ncstat = nf90_put_var(ofid,varid,time)
if (ncstat /= nf90_noerr) call handle_err(ncstat)

! ------------------------------------------------------------------------
! write the regular variables

! --------------------------------------------------
! elevation

write(0,*)'writing elevation'

ncstat = nf90_inq_varid(bfid,'elv',varid)
if (ncstat /= nf90_noerr) call handle_err(ncstat)

ncstat = nf90_get_var(bfid,varid,var2d)
if (ncstat /= nf90_noerr) call handle_err(ncstat)

ncstat = nf90_inq_varid(ofid,'elv',varid)
if (ncstat /= nf90_noerr) call handle_err(ncstat)

ncstat = nf90_put_var(ofid,varid,var2d)
if (ncstat /= nf90_noerr) call handle_err(ncstat)

deallocate(var2d)

! -----------


! i=1
! do yr = 1,numyrs
!   do m = 1,12
!     write(0,*)yr,m,time(i)
!     i = i + 1
!   end do
! end do

! ncstat = nf90_close(bfid)
! if (ncstat /= nf90_noerr) call handle_err(ncstat)
! 
! ncstat = nf90_close(ofid)
! if (ncstat /= nf90_noerr) call handle_err(ncstat)
! 
! stop

! --------------------------------------------------
! temperature

call calcvar('tmp','air',trim(anompath)//tmpfile)

! --------------------------------------------------
! dtr

call calcvar('dtr','dtr',trim(anompath)//dtrfile,llimit=0.)

! --------------------------------------------------
! precipitation

call calcvar('pre','apcp',trim(anompath)//prefile,llimit=0.)

! --------------------------------------------------
! cloud

call calcvar('cld','tcdc',trim(anompath)//cldfile,llimit=0.,ulimit=100.)

! --------------------------------------------------
! wind

call calcvar('wnd','wspd',trim(anompath)//wndfile,llimit=0.)

! --------------------------------------------------
! lightning

call calcvar('lght','lght',trim(anompath)//lghtfile,llimit=0.)

! --------------------------------------------------
! wet days

call calcwetf(wetVprefile)

! ------------------------------------------------------------------------
! close files 

100 continue

write(0,*)
write(0,*)'clean up'

ncstat = nf90_close(bfid)
if (ncstat /= nf90_noerr) call handle_err(ncstat)

ncstat = nf90_close(ofid)
if (ncstat /= nf90_noerr) call handle_err(ncstat)

! ------------------------------------------------------------------------------------------------------------------------

end program makeclimate
