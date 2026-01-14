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

! ---         interannual variability anomalies (can be single file or separate)        ---
! --- variable names default to below but can also be specified in the jobfile namelist ---
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

use parametersmod, only : i2,sp,dp,ndaymon,imissing,rmissing,xlen,ylen,tlen_out,tlen_anom,bfid,ofid,numcyc,offset,tlen_blk
use randomdistmod, only : randomstate,ran_seed,ranur
use netcdfmod,     only : ncstat,handle_err,netcdf_create
use calcvarmod,    only : varinfotype,calcvar
use calcwetmod,    only : calcwetf
use calendarmod,   only : timestruct,ymdt2jd
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
character(200) :: wetfile
character(200) :: lghtfile
character(200) :: varinfofile

type(varinfotype), dimension(7) :: varinfo

! default values for job options in case they are not specified in the namelist

integer :: baseyr    =   1900   ! first year of the transient climate or year after the final year of spinup (standard calendar year)
integer :: numyrs    =   1020   ! number of spinup years to calculate if not specified in the namelist
logical :: transient = .false.  ! generate spinup or transient climate

namelist  / joboptions /                      &
  baseyr,numyrs,transient,basefile,anompath,  &
  tmpfile,dtrfile,prefile,cldfile,wndfile,wetfile,lghtfile,varinfofile

! local variables

! counters
integer :: a,b,c
integer :: i,j
integer :: yr
integer :: m
integer :: t
integer :: t0
integer :: t1
integer :: y0
integer :: y1

! netcdf id's
integer :: afid
integer :: dimid
integer :: varid

! lengths
integer, parameter :: ycyc = 30  ! time block size parameter in years - data will be written out in blocks of this length (or shorter in case of transient climate)
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

logical :: projgrid = .false.

type(timestruct) :: bt
type(timestruct) :: ts

! -----------------------------------------------------------------------------------------------------
! program starts here

! initialize variables with default values (anomaly name, baseline name)

varinfo%ll = -9999.
varinfo%ul = -9999.

varinfo(1) = varinfotype('air','tmp',op='add')
varinfo(2) = varinfotype('dtr','dtr',op='add')
varinfo(3) = varinfotype('apcp','pre',op='add',ll=0.)
varinfo(4) = varinfotype('tcdc','cld',op='add',ll=0.,ul=100.)
varinfo(5) = varinfotype('wspd','wnd',op='add',ll=0.)
varinfo(6) = varinfotype('lght','lght',op='add',ll=0.)
varinfo(7) = varinfotype('wet','wet',op='add',ll=0.,ul=1.)

! -----

! read job options

call getarg(1,jobfile)
call getarg(2,outfile)

open(10,file=jobfile,status='old')
read(10,nml=joboptions)
close(10)

! if (varinfofile /= '') then
! 
!   open(10,file=varinfofile)
!   read(10,*)varinfo
!   close(10)
! 
! end if

do i = 1,size(varinfo)
  write(0,*)i,varinfo(i)
end do

! ------------------------------------------------------------------------
! inquire about the length of the anomaly timeseries and create the offset vector

ncstat = nf90_open(trim(anompath)//tmpfile,nf90_nowrite,afid)
if (ncstat /= nf90_noerr) call handle_err(ncstat)

ncstat = nf90_inq_dimid(afid,'time',dimid)
if (ncstat /= nf90_noerr) call handle_err(ncstat)

ncstat = nf90_inquire_dimension(afid,dimid,len=tlen_anom)
if (ncstat /= nf90_noerr) call handle_err(ncstat)

anomyrs = tlen_anom / 12

write(0,'(a,i0,a)')'anomalies have ',anomyrs,' years of data'

! ------------------------------------------------------------------------

if (transient) then

  ! --------------------------------------------------------------------------------------
  ! transient climate: use the transient interannual variability from the input anomalies
  
  tlen = 12 * ycyc   ! (months) ycyc time block size parameter (set to 30 yrs)

  numcyc = tlen_anom / tlen  ! numcyc number of cycles that will be calculated
  
  if (mod(tlen_anom,tlen) /= 0) numcyc = numcyc + 1

  allocate(seg(numcyc))
  allocate(offset(numcyc))
  allocate(tlen_blk(numcyc))

  do i = 1,numcyc
  
    seg(i) = i
    offset(i) = 1 + tlen * (i - 1)
    tlen_blk(i) = min(tlen,1 + tlen_anom - offset(i))

    write(0,*)seg(i),tlen_blk(i),offset(i),offset(i)+tlen_blk(i)-1

  end do
  
  tlen_out = tlen_anom

  write(0,'(a,i0,a,i0,a,i0,a)')' transient climatology: ',numcyc,' ',ycyc,'-year climate cycles = ',anomyrs,' years of climate'

else

  ! --------------------------------------------------------------------------------------
  ! spinup climate: calculate the random sequence of blocks for repeating output
  
  ! work out the 30-year offsets

  a = anomyrs / ycyc
  b = (anomyrs - 10) / ycyc
  c = (anomyrs - 20) / ycyc

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

  numcyc = numyrs / ycyc  ! this will truncate downwards if the number of years requested is not an even multiple

  if (mod(numyrs,ycyc) /= 0) then
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
  allocate(tlen_blk(numcyc))
  
  tlen_blk = ycyc * 12

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
  
!   baseyr = baseyr - numyrs
  
  tlen_out = numyrs * 12

end if  ! transient or spinup

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

if (xdimname == 'x') then 
  projgrid = .true.
  write(0,*)'NOTE input data are on a projected grid'
end if

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

if (projgrid) then

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

else

  ! ----------------------
  ! lon
  
  ncstat = nf90_inq_varid(bfid,'lon',varid)
  if (ncstat /= nf90_noerr) call handle_err(ncstat)
  
  ncstat = nf90_get_var(bfid,varid,xdimvar)
  if (ncstat /= nf90_noerr) call handle_err(ncstat)
  
  ncstat = nf90_inq_varid(ofid,'lon',varid)
  if (ncstat /= nf90_noerr) call handle_err(ncstat)
  
  ncstat = nf90_put_var(ofid,varid,xdimvar)
  if (ncstat /= nf90_noerr) call handle_err(ncstat)
  
  ! ----------------------
  ! lat
  
  ncstat = nf90_inq_varid(bfid,'lat',varid)
  if (ncstat /= nf90_noerr) call handle_err(ncstat)
  
  ncstat = nf90_get_var(bfid,varid,ydimvar)
  if (ncstat /= nf90_noerr) call handle_err(ncstat)
  
  ncstat = nf90_inq_varid(ofid,'lat',varid)
  if (ncstat /= nf90_noerr) call handle_err(ncstat)
  
  ncstat = nf90_put_var(ofid,varid,ydimvar)
  if (ncstat /= nf90_noerr) call handle_err(ncstat)

end if

! ----------------------
! time

allocate(time(tlen_out))

if (transient) then

  ! just copy the input time vector to the output

  ncstat = nf90_inq_varid(afid,'time',varid)
  if (ncstat /= nf90_noerr) call handle_err(ncstat)
  
  ncstat = nf90_get_var(afid,varid,time)
  if (ncstat /= nf90_noerr) call handle_err(ncstat)

else

  ! start numyrs before the first year of the spinup (baseyr)
  ! calculate day value for start and end year and number of years
      
  bt = timestruct(1950,1,1,0,0,0.)   ! reference all time to 1950-01-01
  
  call ymdt2jd(bt)
  
  y0 = baseyr - numyrs
  
  y1 = baseyr - 1

  if (y0 <= 0 .and. y1 > 0) y0 = y0 - 1  ! need to start an additional year earlier to account for no year zero
  
  t = 1

  do yr = y0,y1
  
    if (yr == 0) cycle
  
    do m = 1,12
    
      ts = timestruct(yr,m,1,0,0,0.)
      
      call ymdt2jd(ts)
      
      time(t) = ts%jd - bt%jd
      
      t = t + 1

    end do
  end do
end if

! ---

ncstat = nf90_close(afid)
if (ncstat /= nf90_noerr) call handle_err(ncstat)

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

! --------------------------------------------------
! temperature

call calcvar(varinfo(1),trim(anompath)//tmpfile)

! --------------------------------------------------
! dtr

call calcvar(varinfo(2),trim(anompath)//dtrfile)

! --------------------------------------------------
! precipitation

call calcvar(varinfo(3),trim(anompath)//prefile)

! --------------------------------------------------
! cloud

call calcvar(varinfo(4),trim(anompath)//cldfile)

! --------------------------------------------------
! wind

call calcvar(varinfo(5),trim(anompath)//wndfile)

! --------------------------------------------------
! lightning

call calcvar(varinfo(6),trim(anompath)//lghtfile)

! --------------------------------------------------
! wet days

call calcvar(varinfo(7),trim(anompath)//wetfile)

! call calcwetf(wetVprefile)

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
