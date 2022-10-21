module calcwetmod

implicit none

contains

! --------------------------------------------------------------------------------


subroutine calcwetf(wetVprefile)

use parametersmod, only : i2,sp,dp,ndaymon,imissing,rmissing,xlen,ylen,tlen_blk,ofid,numcyc,offset
use netcdfmod,     only : ncstat,handle_err
use netcdf

implicit none

character(*), intent(in) :: wetVprefile

integer :: cfid
integer :: ivarid
integer :: ovarid

real(sp), allocatable, dimension(:,:,:) :: a
real(sp), allocatable, dimension(:,:,:) :: b
real(sp), allocatable, dimension(:,:,:) :: pre
real(sp), allocatable, dimension(:,:,:) :: wet

integer(i2), allocatable, dimension(:,:,:) :: ivar
integer(i2), allocatable, dimension(:,:,:) :: ovar

real(sp) :: isf,iao
real(sp) :: osf,oao

integer :: i,j,k,m

real(sp), dimension(2) :: actual_range

character(80) :: status_msg

! --------------------------------------------------------------

write(0,*)
write(0,'(a)')'=== working on wetf ==='

allocate(a(xlen,ylen,12))
allocate(b(xlen,ylen,12))

! ---
! read the precip-wetdays coefficients from an external file

ncstat = nf90_open(wetVprefile,nf90_nowrite,cfid)
if (ncstat /= nf90_noerr) call handle_err(ncstat)

ncstat = nf90_inq_varid(cfid,'a',ivarid)
if (ncstat /= nf90_noerr) call handle_err(ncstat)

ncstat = nf90_get_var(cfid,ivarid,a)
if (ncstat /= nf90_noerr) call handle_err(ncstat)

ncstat = nf90_inq_varid(cfid,'b',ivarid)
if (ncstat /= nf90_noerr) call handle_err(ncstat)

ncstat = nf90_get_var(cfid,ivarid,b)
if (ncstat /= nf90_noerr) call handle_err(ncstat)

ncstat = nf90_close(cfid)
if (ncstat /= nf90_noerr) call handle_err(ncstat)

! ---
! read and write scale factors and add offset

ncstat = nf90_inq_varid(ofid,'pre',ivarid)
if (ncstat /= nf90_noerr) call handle_err(ncstat)

ncstat = nf90_get_att(ofid,ivarid,'scale_factor',isf)
if (ncstat /= nf90_noerr) call handle_err(ncstat)

ncstat = nf90_get_att(ofid,ivarid,'add_offset',iao)
if (ncstat /= nf90_noerr) call handle_err(ncstat)

oao = 0.
osf = 0.001

ncstat = nf90_inq_varid(ofid,'wet',ovarid)
if (ncstat /= nf90_noerr) call handle_err(ncstat)

ncstat = nf90_put_att(ofid,ovarid,'scale_factor',osf)
if (ncstat /= nf90_noerr) call handle_err(ncstat)

ncstat = nf90_put_att(ofid,ovarid,'add_offset',oao)
if (ncstat /= nf90_noerr) call handle_err(ncstat)

allocate(pre(xlen,ylen,12))
allocate(wet(xlen,ylen,12))
allocate(ivar(xlen,ylen,tlen_blk))
allocate(ovar(xlen,ylen,tlen_blk))

! ---
! read transient monthly precipitation calculated previously by makeclimate
! and calculate wet fraction based on the coefficients read above

actual_range = [9999.,-9999.]

do j = 1,numcyc  ! 30 year cycles

  k = offset(j)

  m = 1 + tlen_blk * (j-1)

  write(status_msg,'(a,3i7)')' working on block ',j,k,m
  call overprint(status_msg)

  ovar = imissing

  ncstat = nf90_get_var(ofid,ivarid,ivar,start=[1,1,m],count=[xlen,ylen,tlen_blk])
  if (ncstat /= nf90_noerr) call handle_err(ncstat)  

  do i = 1,tlen_blk,12   ! one whole year at a time
    
    wet = rmissing
    
    where (ivar(:,:,i:11+i) /= imissing)
      
      pre = real(ivar(:,:,i:11+i)) * isf + iao
      wet = min(1.,max(0.,b * pre**a))
      
      ovar(:,:,i:11+i) = nint((wet - oao) / osf)
      
    end where
    
    actual_range(1) = min(actual_range(1),minval(wet,mask = wet /= rmissing))
    actual_range(2) = max(actual_range(2),maxval(wet,mask = wet /= rmissing))
    
  end do

  ncstat = nf90_put_var(ofid,ovarid,ovar,start=[1,1,m],count=[xlen,ylen,tlen_blk])
  if (ncstat /= nf90_noerr) call handle_err(ncstat)  

end do

write(0,*)
write(0,'(a,2f12.5)')' writing actual range for wet',actual_range

ncstat = nf90_put_att(ofid,ovarid,'actual_range',actual_range)
if (ncstat /= nf90_noerr) call handle_err(ncstat)  
    
! ---

end subroutine calcwetf

! --------------------------------------------------------------------------------

end module calcwetmod
