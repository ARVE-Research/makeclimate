module netcdfmod

use netcdf

implicit none

public :: handle_err
public :: netcdf_create

integer, parameter :: dl = 1

integer :: ncstat

contains

!--------------------------------------------------------------------------------

subroutine handle_err(ncstat)

! Internal subroutine - checks error ncstat after each netcdf call,
! prints out text message each time an error code is returned.

integer, intent (in) :: ncstat

if(ncstat /= nf90_noerr) then
  write(0,*)trim(nf90_strerror(ncstat))
  stop
end if

end subroutine handle_err

!--------------------------------------------------------------------------------

subroutine netcdf_create(filename,nx,ny,nt,ncid)

use parametersmod, only : sp,dp,imissing

implicit none

character(*), intent(in)  :: filename
integer,      intent(in)  :: nx
integer,      intent(in)  :: ny
integer,      intent(in)  :: nt
integer,      intent(out) :: ncid

integer, dimension(3) :: dimid
integer :: varid

integer, dimension(3) :: cs

integer :: cx
integer :: cy
integer :: ct

! ------------------------

write(0,*)'creating outputfile ',trim(filename),' with:'
write(0,*)'xlen',nx
write(0,*)'ylen',ny
write(0,*)'tlen',nt

cx = min(nx,60)
cy = min(ny,60)
ct = min(nt,360)

cs = [cx,cy,ct]

ncstat = nf90_create(filename,nf90_netcdf4,ncid)
if (ncstat /= nf90_noerr) call handle_err(ncstat)

! -------
! dimensions

ncstat = nf90_def_dim(ncid,'lon',nx,dimid(1))
if (ncstat /= nf90_noerr) call handle_err(ncstat)

ncstat = nf90_def_dim(ncid,'lat',nx,dimid(2))
if (ncstat /= nf90_noerr) call handle_err(ncstat)

ncstat = nf90_def_dim(ncid,'time',nt,dimid(3))
if (ncstat /= nf90_noerr) call handle_err(ncstat)

! -------
! dimension variables

ncstat = nf90_def_var(ncid,'lon',nf90_double,dimid(1),varid)
if (ncstat /= nf90_noerr) call handle_err(ncstat)

ncstat = nf90_put_att(ncid,varid,'long_name','longitude')
if (ncstat /= nf90_noerr) call handle_err(ncstat)

ncstat = nf90_put_att(ncid,varid,'units','degrees_east')
if (ncstat /= nf90_noerr) call handle_err(ncstat)

ncstat = nf90_put_att(ncid,varid,'actual_range',[0._dp,0._dp])
if (ncstat /= nf90_noerr) call handle_err(ncstat)

! ---

ncstat = nf90_def_var(ncid,'lat',nf90_double,dimid(2),varid)
if (ncstat /= nf90_noerr) call handle_err(ncstat)

ncstat = nf90_put_att(ncid,varid,'long_name','latitude')
if (ncstat /= nf90_noerr) call handle_err(ncstat)

ncstat = nf90_put_att(ncid,varid,'units','degrees_north')
if (ncstat /= nf90_noerr) call handle_err(ncstat)

ncstat = nf90_put_att(ncid,varid,'actual_range',[0._dp,0._dp])
if (ncstat /= nf90_noerr) call handle_err(ncstat)

! ---

ncstat = nf90_def_var(ncid,'time',nf90_double,dimid(3),varid)
if (ncstat /= nf90_noerr) call handle_err(ncstat)

ncstat = nf90_put_att(ncid,varid,'long_name','time')
if (ncstat /= nf90_noerr) call handle_err(ncstat)

ncstat = nf90_put_att(ncid,varid,'units','days since 0001-01-01')
if (ncstat /= nf90_noerr) call handle_err(ncstat)

ncstat = nf90_put_att(ncid,varid,'actual_range',[0._dp,0._dp])
if (ncstat /= nf90_noerr) call handle_err(ncstat)

ncstat = nf90_put_att(ncid,varid,'avg_period','0000-01-00')
if (ncstat /= nf90_noerr) call handle_err(ncstat)

ncstat = nf90_put_att(ncid,varid,'delta_t','0000-01-00')
if (ncstat /= nf90_noerr) call handle_err(ncstat)

ncstat = nf90_put_att(ncid,varid,'coordinate_defines','start')
if (ncstat /= nf90_noerr) call handle_err(ncstat)

ncstat = nf90_put_att(ncid,varid,'calendar','proleptic_gregorian')
if (ncstat /= nf90_noerr) call handle_err(ncstat)

! -------
! regular variables

call defvar(ncid,dimid(1:2),'elv','elevation above mean sea level','m',1._sp,0._sp,cs(1:2))
call defvar(ncid,dimid,'tmp','2m air temperature','degC',0.1_sp,0._sp,cs)
call defvar(ncid,dimid,'dtr','diurnal temperature range','degC',0.1_sp,0._sp,cs)
call defvar(ncid,dimid,'pre','accumulated precipitation','mm',0.1_sp,0._sp,cs)
call defvar(ncid,dimid,'wet','fraction of days with >0.1 mm precipitation per month','fraction',0.0001_sp,0._sp,cs)
call defvar(ncid,dimid,'cld','total column cloud cover','percent',0.01_sp,0._sp,cs)
call defvar(ncid,dimid,'wnd','10m windspeed','m s-1',0.01_sp,0._sp,cs)
call defvar(ncid,dimid,'lght','lightning strike density','km-2 d-1',0.0001_sp,3.2767_sp,cs)

! ----

end subroutine netcdf_create

!--------------------------------------------------------------------------------

subroutine defvar(ncid,dimid,varname,longname,units,sf,ao,cs)

use parametersmod, only : sp,imissing

implicit none

integer,      intent(in) :: ncid
character(*), intent(in) :: varname
character(*), intent(in) :: longname
character(*), intent(in) :: units
real(sp),     intent(in) :: sf
real(sp),     intent(in) :: ao

integer, dimension(:), intent(in) :: dimid
integer, dimension(:), intent(in) :: cs

integer :: varid


! ----

ncstat = nf90_def_var(ncid,varname,nf90_short,dimid,varid,chunksizes=cs,deflate_level=dl)
if (ncstat /= nf90_noerr) call handle_err(ncstat)

ncstat = nf90_put_att(ncid,varid,'long_name',longname)
if (ncstat /= nf90_noerr) call handle_err(ncstat)

ncstat = nf90_put_att(ncid,varid,'units',units)
if (ncstat /= nf90_noerr) call handle_err(ncstat)

ncstat = nf90_put_att(ncid,varid,'scale_factor',sf)
if (ncstat /= nf90_noerr) call handle_err(ncstat)

ncstat = nf90_put_att(ncid,varid,'add_offset',ao)
if (ncstat /= nf90_noerr) call handle_err(ncstat)

ncstat = nf90_put_att(ncid,varid,'_FillValue',imissing)
if (ncstat /= nf90_noerr) call handle_err(ncstat)

ncstat = nf90_put_att(ncid,varid,'missing_value',imissing)
if (ncstat /= nf90_noerr) call handle_err(ncstat)

ncstat = nf90_put_att(ncid,varid,'actual_range',[0._sp,0._sp])
if (ncstat /= nf90_noerr) call handle_err(ncstat)

ncstat = nf90_put_att(ncid,varid,'source','unknown')
if (ncstat /= nf90_noerr) call handle_err(ncstat)

end subroutine defvar

!--------------------------------------------------------------------------------

end module netcdfmod