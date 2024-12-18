module calcvarmod

use parametersmod, only : sp

implicit none

type varinfotype
  character(10) :: aname        ! name of the variable in the anomaly file
  character(10) :: bname        ! name of the variable in the baseline file
  character(10) :: op = 'add'   ! operation to perform; defaults to add
  real :: ll = -9999.           ! lower limit of the variable value, initialize as undefined
  real :: ul = -9999.           ! upper limit of the variable value, initialize as undefined
end type varinfotype

contains

! --------------------------------------------------------------------------------

! function op(a,b,opname)
! 
! real(sp), dimension(:,:,:), intent(in) :: a
! real(sp), dimension(:,:,:), intent(in) :: b
! character(*),               intent(in) :: opname
! 
! real(sp), allocatable, dimension(:,:,:) :: op
! 
! allocate(op(size(a)))
! 
! if (opname == 'add') then
! 
!   op = a + b
! 
! else if (opname == 'mul') then
! 
!   op = a * b
! 
! else
!   
!   stop 'error: operator not defined'
! 
! end if
! 
! end function op  

! --------------------------------------------------------------------------------

subroutine calcvar(varinfo,anomfile)

use parametersmod, only : i2,sp,dp,ndaymon,imissing,rmissing,xlen,ylen,tlen_blk,tlen_anom,bfid,ofid,numcyc,offset
use netcdfmod,     only : ncstat,handle_err
use netcdf

implicit none

type(varinfotype), intent(in) :: varinfo    ! 
character(*),      intent(in) :: anomfile   ! anomaly filename

! --

character(10) :: bvarname   ! baseline file variable name
character(10) :: avarname   ! anomaly file variable name

integer :: a,b
integer :: i,j,k,m

real(sp) :: scale_factor
real(sp) :: add_offset

real(sp), dimension(2) :: actual_range
 
integer(i2), allocatable, dimension(:,:,:) :: base
integer(i2), allocatable, dimension(:,:,:) :: vout

real(sp),    allocatable, dimension(:,:,:) :: anom
real(sp),    allocatable, dimension(:,:,:) :: rbase
real(sp),    allocatable, dimension(:,:,:) :: rvout

integer :: afid

integer(i2) :: missing

character(80) :: status_msg
integer :: varid

! real(sp) :: llimit
! real(sp) :: ulimit

! --------------------------------------------------------------

avarname = varinfo%aname
bvarname = varinfo%bname

write(0,*)
write(0,'(a,a,a)')'=== working on ',trim(bvarname),' ==='

allocate(base(xlen,ylen,12))
allocate(rbase(xlen,ylen,12))

allocate(anom(xlen,ylen,tlen_anom))

! ---
! baseline

ncstat = nf90_inq_varid(bfid,bvarname,varid)
if (ncstat /= nf90_noerr) call handle_err(ncstat)

ncstat = nf90_get_att(bfid,varid,'missing_value',missing)
if (ncstat /= nf90_noerr) call handle_err(ncstat)

ncstat = nf90_get_att(bfid,varid,'scale_factor',scale_factor)
if (ncstat /= nf90_noerr) call handle_err(ncstat)

ncstat = nf90_get_att(bfid,varid,'add_offset',add_offset)
if (ncstat /= nf90_noerr) call handle_err(ncstat)

if (trim(bvarname) == 'cld') scale_factor = 0.01

write(0,'(a,f0.5,a,f0.5)')' read baseline: ',scale_factor,' ',add_offset

ncstat = nf90_get_var(bfid,varid,base)
if (ncstat /= nf90_noerr) call handle_err(ncstat)

where (base /= missing) rbase = real(base) * scale_factor + add_offset

! ---
! anomaly

write(status_msg,'(a,i0)')' reading anomalies... ',tlen_anom
call overprint(status_msg)

ncstat = nf90_open(anomfile,nf90_nowrite,afid)
if (ncstat /= nf90_noerr) call handle_err(ncstat)

ncstat = nf90_inq_varid(afid,avarname,varid)
if (ncstat /= nf90_noerr) call handle_err(ncstat)

ncstat = nf90_get_var(afid,varid,anom)
if (ncstat /= nf90_noerr) call handle_err(ncstat)

ncstat = nf90_close(afid)
if (ncstat /= nf90_noerr) call handle_err(ncstat)

write(status_msg,'(a)')' reading anomalies... done'
call overprint(status_msg)
write(0,*)

! ---
! output

ncstat = nf90_inq_varid(ofid,bvarname,varid)
if (ncstat /= nf90_noerr) call handle_err(ncstat)

! ---
! calculate and write in 30-year blocks

actual_range = [9999.,-9999.]


!work in blocks of 30 years = 2160 months = tlen

! write(0,*)'creating 30-year-blocks'

do j = 1,numcyc

  k = offset(j)  ! starting month
  
  allocate(vout(xlen,ylen,tlen_blk(j)))
  allocate(rvout(xlen,ylen,tlen_blk(j)))
  
  rvout = rmissing
  
  ! operate on one year of data at a time

  do i = 1,tlen_blk(j),12

    a = k + i - 1  ! anomaly start month
    b = a + 11     ! anomaly end month
    
    ! write(0,*)k,i,a,b
    
!     where (base /= missing) rvout(:,:,i:11+i) = op(rbase,anom(:,:,a:b),varinfo%op)
    
    if (varinfo%op == 'mul') then
    
      anom(:,:,a:b) = min(max(anom(:,:,a:b),0.),100.)  ! limit the range because there are some spurious values
    
      where (base /= missing) rvout(:,:,i:11+i) = rbase * anom(:,:,a:b)
      
    else  ! default to add
    
      where (base /= missing) rvout(:,:,i:11+i) = rbase + anom(:,:,a:b)

    end if
    
  end do

  if (j == 1) then
    m = 1
  else
    m = 1 + sum(tlen_blk(1:j-1))  ! start value for write
  end if

  write(status_msg,'(a,i0,a,i0,a,i0,a,i0)')' working on block ',j,' ',k,' ',m,' ',tlen_blk(j)
  call overprint(status_msg)
  
  if (varinfo%ll > -9999.) then
    where (rvout /= rmissing) rvout = max(rvout,varinfo%ll)
  end if
    
  if (varinfo%ul > -9999.) then
    where (rvout /= rmissing) rvout = min(rvout,varinfo%ul)
  end if
  
  where (rvout /= rmissing)
    vout = nint((rvout - add_offset) / scale_factor,kind=i2)
  elsewhere
    vout = imissing
  end where  

  ncstat = nf90_put_var(ofid,varid,vout,start=[1,1,m],count=[xlen,ylen,tlen_blk(j)])
  if (ncstat /= nf90_noerr) call handle_err(ncstat)  

  actual_range(1) = min(actual_range(1),minval(rvout,mask = rvout /= rmissing))
  actual_range(2) = max(actual_range(2),maxval(rvout,mask = rvout /= rmissing))

  ! ---
  
  deallocate(vout)
  deallocate(rvout)

end do

write(0,*)
write(0,'(a,a,a,f0.5,a,f0.5)')' writing actual range for ',trim(bvarname),' ',actual_range(1),' ',actual_range(2)

ncstat = nf90_put_att(ofid,varid,'scale_factor',scale_factor)
if (ncstat /= nf90_noerr) call handle_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'add_offset',add_offset)
if (ncstat /= nf90_noerr) call handle_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'actual_range',actual_range)
if (ncstat /= nf90_noerr) call handle_err(ncstat)

end subroutine calcvar

!--------------------------------------------------------------------------------

end module calcvarmod
