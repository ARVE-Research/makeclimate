module calcvarmod

implicit none

contains

!--------------------------------------------------------------------------------

subroutine calcvar(bvarname,avarname,anomfile,llimit,ulimit)

use parametersmod, only : i2,sp,dp,ndaymon,imissing,rmissing,xlen,ylen,tlen_blk,tlen_anom,bfid,ofid,numcyc,offset
use netcdfmod,     only : ncstat,handle_err
use netcdf

implicit none

character(*), intent(in) :: bvarname
character(*), intent(in) :: avarname
character(*), intent(in) :: anomfile

real(sp),     intent(in), optional :: llimit
real(sp),     intent(in), optional :: ulimit

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

!--------------------------------------------------------------

write(0,*)
write(0,'(a,a,a)')'=== working on ',trim(bvarname),' ==='

allocate(base(xlen,ylen,12))
allocate(rbase(xlen,ylen,12))

allocate(anom(xlen,ylen,tlen_anom))
allocate(vout(xlen,ylen,tlen_blk))
allocate(rvout(xlen,ylen,tlen_blk))

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

write(0,'(a,2f8.5)')' read baseline',scale_factor,add_offset

ncstat = nf90_get_var(bfid,varid,base)
if (ncstat /= nf90_noerr) call handle_err(ncstat)

where (base /= missing) rbase = real(base) * scale_factor + add_offset

! ---
! anomaly

write(status_msg,'(a)')' reading anomalies... '
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

rvout = rmissing

!work in blocks of 30 years = 2160 months = tlen

! write(0,*)'creating 30-year-blocks'

do j = 1,numcyc

  k = offset(j)

  do i = 1,tlen_blk,12

    a = k + i - 1
    b = a + 11
    
    where (base /= missing) rvout(:,:,i:11+i) = rbase + anom(:,:,a:b)
    
  end do

  m = 1 + tlen_blk * (j-1)

  write(status_msg,'(a,3i7)')' working on block ',j,k,m
  call overprint(status_msg)
  
  if (present(llimit)) then
    where (rvout /= rmissing) rvout = max(rvout,llimit)
  end if
    
  if (present(ulimit)) then
    where (rvout /= rmissing) rvout = min(rvout,ulimit)
  end if
  
  where (rvout /= rmissing)
    vout = nint((rvout - add_offset) / scale_factor)
  elsewhere
    vout = imissing
  end where  

  ncstat = nf90_put_var(ofid,varid,vout,start=[1,1,m],count=[xlen,ylen,tlen_blk])
  if (ncstat /= nf90_noerr) call handle_err(ncstat)  

  ! write(0,*)j,k,m,actual_range

  actual_range(1) = min(actual_range(1),minval(rvout,mask = rvout /= rmissing))
  actual_range(2) = max(actual_range(2),maxval(rvout,mask = rvout /= rmissing))

end do

write(0,*)
write(0,'(a,a,2f12.5)')' writing actual range for ',bvarname,actual_range

ncstat = nf90_put_att(ofid,varid,'scale_factor',scale_factor)
if (ncstat /= nf90_noerr) call handle_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'add_offset',add_offset)
if (ncstat /= nf90_noerr) call handle_err(ncstat)

ncstat = nf90_put_att(ofid,varid,'actual_range',actual_range)
if (ncstat /= nf90_noerr) call handle_err(ncstat)

end subroutine calcvar

!--------------------------------------------------------------------------------

end module calcvarmod
