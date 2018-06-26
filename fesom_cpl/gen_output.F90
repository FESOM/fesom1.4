subroutine init_output
  use o_mesh
  use o_passive_tracer_mod
  use o_age_tracer_mod
  use g_config
  use g_diag
  use g_clock
  use g_PARFE
  implicit none

#include "netcdf.inc" 

  integer                   :: status, ncid, dimid_rec, j
  integer                   :: dimid_2d, dimid_3d, dimids(2)
  integer                   :: time_varid, iter_varid
  integer                   :: ssh_varid, tra_varid(num_tracer)
  integer                   :: u_varid, v_varid, w_varid, wpot_varid
  integer                   :: area_varid, hice_varid, hsnow_varid
  integer                   :: uice_varid, vice_varid
  character(100)            :: longname
  character(2000)           :: filename
  character(1)              :: trind

  if(yearnew==yearold) return

  ! initialize the counter for saving results
  save_count=1
  save_count_restart=1
  seconds_from_yearstart=0

  if (mype/=0) return

  write(*,*) 'initialize new output files'

  ! first, snapshots

  ! ocean

  filename=trim(ResultPath)//runid//'.'//cyearnew//'.oce.nc'

  ! create a file
  status = nf_create(filename, nf_clobber, ncid)
  if (status.ne.nf_noerr) call handle_err(status)

  ! Define the dimensions
  status = nf_def_dim(ncid, 'nodes_2d', nod2d, dimid_2d)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_dim(ncid, 'nodes_3d', nod3d, dimid_3d)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_dim(ncid, 'time', NF_UNLIMITED, dimid_rec)
  if (status .ne. nf_noerr) call handle_err(status)

  ! Define the time and iteration variables
  status = nf_def_var(ncid, 'time', NF_DOUBLE, 1, dimid_rec, time_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_var(ncid, 'iter', NF_INT, 1, dimid_rec, iter_varid)
  if (status .ne. nf_noerr) call handle_err(status)

  ! Define the netCDF variables for 2D fields.
  ! In Fortran, the unlimited dimension must come
  ! last on the list of dimids.
  dimids(1) = dimid_2d
  dimids(2) = dimid_rec

  status = nf_def_var(ncid, 'ssh', NF_DOUBLE, 2, dimids, ssh_varid)
  if (status .ne. nf_noerr) call handle_err(status)

  ! Define the netCDF variables for 3D fields
  dimids(1) = dimid_3d
  dimids(2) = dimid_rec

  status = nf_def_var(ncid, 'u', NF_DOUBLE, 2, dimids, u_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_var(ncid, 'v', NF_DOUBLE, 2, dimids, v_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_var(ncid, 'w', NF_DOUBLE, 2, dimids, w_varid)
  if (status .ne. nf_noerr) call handle_err(status)
#ifndef use_non_hydrostatic
  status = nf_def_var(ncid, 'wpot', NF_DOUBLE, 2, dimids, wpot_varid)
  if (status .ne. nf_noerr) call handle_err(status)
#endif
  status = nf_def_var(ncid, 'temp', NF_DOUBLE, 2, dimids, tra_varid(1))
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_var(ncid, 'salt', NF_DOUBLE, 2, dimids, tra_varid(2))
  if (status .ne. nf_noerr) call handle_err(status)

  if(use_passive_tracer) then
     do j=1,num_passive_tracer
        write(trind,'(i1)') j
        status = nf_def_var(ncid, 'ptr'//trind, NF_DOUBLE, 2, dimids, &
             tra_varid(index_passive_tracer(j)))
        if (status .ne. nf_noerr) call handle_err(status)
     end do
  end if

  if(use_age_tracer) then
     do j=1,num_age_tracer
        write(trind,'(i1)') j
        status = nf_def_var(ncid, 'age'//trind, NF_DOUBLE, 2, dimids, &
             tra_varid(index_age_tracer(j)))
        if (status .ne. nf_noerr) call handle_err(status)
     end do
  end if

  ! Assign long_name and units attributes to variables.
  longname='model time'
  status = nf_PUT_ATT_TEXT(ncid, time_varid, 'long_name', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, time_varid, 'units', 1, 's')
  if (status .ne. nf_noerr) call handle_err(status)
  longname='iteration_count'
  status = nf_PUT_ATT_TEXT(ncid, iter_varid, 'long_name', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)

  longname='sea surface elevation'
  status = nf_put_att_text(ncid, ssh_varid, 'description', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, ssh_varid, 'units', 1, 'm')
  if (status .ne. nf_noerr) call handle_err(status)
  longname='zonal velocity'
  status = nf_put_att_text(ncid, u_varid, 'description', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, u_varid, 'units', 3, 'm/s')
  if (status .ne. nf_noerr) call handle_err(status)
  longname='meridional velocity'
  status = nf_put_att_text(ncid, v_varid, 'description', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, v_varid, 'units', 3, 'm/s')
  if (status .ne. nf_noerr) call handle_err(status)
  longname='vertical velocity'
  status = nf_put_att_text(ncid, w_varid, 'description', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, w_varid, 'units', 3, 'm/s')
  if (status .ne. nf_noerr) call handle_err(status)
#ifndef use_non_hydrostatic
  longname='vertical velocity potential'
  status = nf_put_att_text(ncid, wpot_varid, 'description', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, wpot_varid, 'units', 5, 'm.m/s')
  if (status .ne. nf_noerr) call handle_err(status)
#endif
  longname='potential temperature'
  status = nf_put_att_text(ncid, tra_varid(1), 'description', len_trim(longname), trim(longname))
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, tra_varid(1), 'units', 4, 'degC')
  if (status .ne. nf_noerr) call handle_err(status)
  longname='salinity'
  status = nf_put_att_text(ncid, tra_varid(2), 'description', len_trim(longname), longname) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, tra_varid(2), 'units', 3, 'psu')
  if (status .ne. nf_noerr) call handle_err(status)

  if(use_passive_tracer) then
     do j=1,num_passive_tracer
        write(trind,'(i1)') j
        longname='passive tracer '//trind
        status = nf_put_att_text(ncid, tra_varid(index_passive_tracer(j)), &
             'description', len_trim(longname), longname) 
        if (status .ne. nf_noerr) call handle_err(status)
        !status = nf_put_att_text(ncid, tra_varid(index_passive_tracer(j)), &
        !'units', 3, 'NaN')
        !if (status .ne. nf_noerr) call handle_err(status)
     end do
  end if

  if(use_age_tracer) then
     do j=1,num_age_tracer
        write(trind,'(i1)') j
        longname='age tracer '//trind
        status = nf_put_att_text(ncid, tra_varid(index_age_tracer(j)), &
             'description', len_trim(longname), longname) 
        if (status .ne. nf_noerr) call handle_err(status)
        status = nf_put_att_text(ncid, tra_varid(index_age_tracer(j)), &
             'units', 4, 'Year')
        if (status .ne. nf_noerr) call handle_err(status)
     end do
  end if

  status = nf_enddef(ncid)
  if (status .ne. nf_noerr) call handle_err(status)

  status=nf_close(ncid)
  if (status .ne. nf_noerr) call handle_err(status)


  ! ice

#ifdef use_ice
  filename=trim(ResultPath)//runid//'.'//cyearnew//'.ice.nc'

  ! create a file
  status = nf_create(filename, nf_clobber, ncid)
  if (status.ne.nf_noerr) call handle_err(status)

  ! Define the dimensions
  status = nf_def_dim(ncid, 'nodes_2d', nod2d, dimid_2d)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_dim(ncid, 'time', NF_UNLIMITED, dimid_rec)
  if (status .ne. nf_noerr) call handle_err(status)

  ! Define the time and iteration variables
  status = nf_def_var(ncid, 'time', NF_DOUBLE, 1, dimid_rec, time_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_var(ncid, 'iter', NF_INT, 1, dimid_rec, iter_varid)
  if (status .ne. nf_noerr) call handle_err(status)

  ! Define the netCDF variables for 2D fields.
  ! In Fortran, the unlimited dimension must come
  ! last on the list of dimids.
  dimids(1) = dimid_2d
  dimids(2) = dimid_rec

  status = nf_def_var(ncid, 'area', NF_DOUBLE, 2, dimids, area_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_var(ncid, 'hice', NF_DOUBLE, 2, dimids, hice_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_var(ncid, 'hsnow', NF_DOUBLE, 2, dimids, hsnow_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_var(ncid, 'uice', NF_DOUBLE, 2, dimids, uice_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_var(ncid, 'vice', NF_DOUBLE, 2, dimids, vice_varid)
  if (status .ne. nf_noerr) call handle_err(status)

  ! Assign long_name and units attributes to variables.
  longname='model time'
  status = nf_PUT_ATT_TEXT(ncid, time_varid, 'long_name', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, time_varid, 'units', 1, 's')
  if (status .ne. nf_noerr) call handle_err(status)
  longname='iteration_count'
  status = nf_PUT_ATT_TEXT(ncid, iter_varid, 'long_name', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)

  longname='ice concentration [0 to 1]'
  status = nf_PUT_ATT_TEXT(ncid, area_varid, 'description', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)
  longname='effective ice thickness'
  status = nf_PUT_ATT_TEXT(ncid, hice_varid, 'description', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, hice_varid, 'units', 1, 'm')
  if (status .ne. nf_noerr) call handle_err(status)
  longname='effective snow thickness'
  status = nf_PUT_ATT_TEXT(ncid, hsnow_varid, 'description', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, hsnow_varid, 'units', 1, 'm')
  if (status .ne. nf_noerr) call handle_err(status)
  longname='zonal velocity'
  status = nf_PUT_ATT_TEXT(ncid, uice_varid, 'description', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, uice_varid, 'units', 3, 'm/s')
  if (status .ne. nf_noerr) call handle_err(status)
  longname='meridional velocity'
  status = nf_PUT_ATT_TEXT(ncid, vice_varid, 'description', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, vice_varid, 'units', 3, 'm/s')
  if (status .ne. nf_noerr) call handle_err(status)

  status = nf_enddef(ncid)
  if (status .ne. nf_noerr) call handle_err(status)

  status=nf_close(ncid)
  if (status .ne. nf_noerr) call handle_err(status)

#endif


  ! second, mean fields
#if defined(allow_calcmeans) || defined(allow_diag)
  call init_output_mean
#endif

  ! third, mesh diagnose
#ifdef allow_diag
  if(diag_mesh) then
     call init_output_mesh
  end if
#endif

end subroutine init_output
!
!----------------------------------------------------------------------------
!
subroutine init_output_mean
  use o_param
  use o_mesh
  use o_passive_tracer_mod
  use o_age_tracer_mod
  use g_config
  use g_diag
  use g_clock
  use g_PARFE
  implicit none

#include "netcdf.inc" 

  integer                   :: status, ncid, dimid_rec, j
  integer                   :: dimid_2d, dimid_3d, dimids(2)
  integer                   :: time_varid, iter_varid
  integer                   :: Kv_varid
  integer                   :: sgs_u_varid, sgs_v_varid, sgs_ut_varid
  integer                   :: sgs_vt_varid, sgs_us_varid, sgs_vs_varid
  character(100)            :: longname
  character(2000)           :: filename
  character(100)            :: att_text
  character(1)              :: trind

  if (mype/=0) return

  ! diagnose
#ifdef allow_diag

  ! ocean
  if(diag_oce) then

     filename=trim(ResultPath)//runid//'.'//cyearnew//'.oce.diag.nc'

     ! create a file
     status = nf_create(filename, nf_clobber, ncid)
     if (status.ne.nf_noerr) call handle_err(status)

     ! Define the dimensions
     status = nf_def_dim(ncid, 'nodes_2d', nod2d, dimid_2d)
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_def_dim(ncid, 'nodes_3d', nod3d, dimid_3d)
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_def_dim(ncid, 'time', NF_UNLIMITED, dimid_rec)
     if (status .ne. nf_noerr) call handle_err(status)

     ! Define the time and iteration variables
     status = nf_def_var(ncid, 'time', NF_DOUBLE, 1, dimid_rec, time_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_def_var(ncid, 'iter', NF_INT, 1, dimid_rec, iter_varid)
     if (status .ne. nf_noerr) call handle_err(status)

     ! In Fortran, the unlimited dimension must come
     ! last on the list of dimids.
     ! Define the netCDF variables for 3D fields
     dimids(1) = dimid_3d
     dimids(2) = dimid_rec

     if(Redi_GM .and. diag_oce_GM_vel) then
        status = nf_def_var(ncid, 'sgs_u', NF_FLOAT, 2, dimids, sgs_u_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        status = nf_def_var(ncid, 'sgs_v', NF_FLOAT, 2, dimids, sgs_v_varid)
        if (status .ne. nf_noerr) call handle_err(status)
     endif

     if(diag_oce_SGS_transp) then
        status = nf_def_var(ncid, 'sgs_ut', NF_FLOAT, 2, dimids, sgs_ut_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        status = nf_def_var(ncid, 'sgs_vt', NF_FLOAT, 2, dimids, sgs_vt_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        status = nf_def_var(ncid, 'sgs_us', NF_FLOAT, 2, dimids, sgs_us_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        status = nf_def_var(ncid, 'sgs_vs', NF_FLOAT, 2, dimids, sgs_vs_varid)
        if (status .ne. nf_noerr) call handle_err(status)
     end if

     if(diag_oce_Kv) then
        status = nf_def_var(ncid, 'Kv', NF_FLOAT, 2, dimids, Kv_varid)
        if (status .ne. nf_noerr) call handle_err(status)
     end if

!  Assign long_name and units attributes to variables.
 
!  longname='model time'
!  status = nf_PUT_ATT_TEXT(ncid, time_varid, 'long_name', len_trim(longname), trim(longname)) 
!  if (status .ne. nf_noerr) call handle_err(status)

!  status = nf_put_att_text(ncid, time_varid, 'units', 1, 's')
!  if (status .ne. nf_noerr) call handle_err(status)

! NetCDF Climate and Forecast (CF) Metadata Conventions will be used instead of that commented above
     longname='time'
     status = nf_PUT_ATT_TEXT(ncid, time_varid, 'long_name', len_trim(longname), trim(longname))
     if (status .ne. nf_noerr) call handle_err(status)
  write(att_text, '(a14,I4.4,a1,I2.2,a1,I2.2,a6)') 'seconds since ', year_start, '-', month_start, '-', day_start, ' 0:0:0'
     status = nf_PUT_ATT_TEXT(ncid, time_varid, 'units', len_trim(att_text), trim(att_text))
     if (include_fleapyear) then
	att_text='standard'
     else
        att_text='noleap'
     end if
     status = nf_put_att_text(ncid, time_varid, 'calendar', len_trim(att_text), trim(att_text))
     if (status .ne. nf_noerr) call handle_err(status)

     longname='iteration_count'
     status = nf_PUT_ATT_TEXT(ncid, iter_varid, 'long_name', len_trim(longname), trim(longname)) 
     if (status .ne. nf_noerr) call handle_err(status)

     if(Redi_GM .and. diag_oce_GM_vel) then
        longname='SGS (GM) zonal velocity integrated from bottom (k*S_x)'
        status = nf_put_att_text(ncid, sgs_u_varid, 'description', len_trim(longname), trim(longname))
        if (status .ne. nf_noerr) call handle_err(status)
        status = nf_put_att_text(ncid, sgs_u_varid, 'units', 5, 'm^2/s')
        if (status .ne. nf_noerr) call handle_err(status)
        longname='SGS (GM) meridional velocity integrated from bottom (k*S_y)'
        status = nf_put_att_text(ncid, sgs_v_varid, 'description', len_trim(longname), trim(longname))
        if (status .ne. nf_noerr) call handle_err(status)
        status = nf_put_att_text(ncid, sgs_v_varid, 'units', 5, 'm^2/s')
        if (status .ne. nf_noerr) call handle_err(status)
     endif

     if(diag_oce_SGS_transp) then
        longname='SGS zonal temperature flux'
        status = nf_put_att_text(ncid, sgs_ut_varid, 'description', len_trim(longname), trim(longname))
        if (status .ne. nf_noerr) call handle_err(status)
        status = nf_put_att_text(ncid, sgs_ut_varid, 'units', 8, 'degC m/s')
        if (status .ne. nf_noerr) call handle_err(status)
        longname='SGS meridional temperature flux'
        status = nf_put_att_text(ncid, sgs_vt_varid, 'description', len_trim(longname), trim(longname))
        if (status .ne. nf_noerr) call handle_err(status)
        status = nf_put_att_text(ncid, sgs_vt_varid, 'units', 8, 'degC m/s')
        if (status .ne. nf_noerr) call handle_err(status)
        longname='SGS zonal salinity flux'
        status = nf_put_att_text(ncid, sgs_us_varid, 'description', len_trim(longname), trim(longname))
        if (status .ne. nf_noerr) call handle_err(status)
        status = nf_put_att_text(ncid, sgs_us_varid, 'units', 7, 'psu m/s')
        if (status .ne. nf_noerr) call handle_err(status)
        longname='SGS meridional salinity flux'
        status = nf_put_att_text(ncid, sgs_vs_varid, 'description', len_trim(longname), trim(longname))
        if (status .ne. nf_noerr) call handle_err(status)
        status = nf_put_att_text(ncid, sgs_vs_varid, 'units', 7, 'psu m/s')
        if (status .ne. nf_noerr) call handle_err(status)
     end if

     if(diag_oce_Kv) then
        longname='Instantaneous vertical diffusivity'
        status = nf_put_att_text(ncid, Kv_varid, 'description', len_trim(longname), trim(longname))
        if (status .ne. nf_noerr) call handle_err(status)
        status = nf_put_att_text(ncid, Kv_varid, 'units', 5, 'm^2/s')
        if (status .ne. nf_noerr) call handle_err(status)
     end if

     status = nf_enddef(ncid)  !end def
     if (status .ne. nf_noerr) call handle_err(status)

     status=nf_close(ncid)     !close file
     if (status .ne. nf_noerr) call handle_err(status)

  endif  !ocean

#endif

end subroutine init_output_mean


subroutine write_snapshots
  use o_array
  use o_mesh
  use o_passive_tracer_mod
  use o_age_tracer_mod
  use i_array
  use g_config
  use g_clock
  use g_PARFE
  implicit none

#include "netcdf.inc" 

  integer                   :: status, ncid, j
  integer                   :: time_varid, iter_varid, dimid_rec
  integer                   :: ssh_varid, tra_varid(num_tracer)
  integer                   :: u_varid, v_varid, w_varid, wpot_varid
  integer                   :: area_varid, hice_varid, hsnow_varid
  integer                   :: uice_varid, vice_varid
  integer                   :: start(2), count(2), n3
  character(2000)           :: filename
  character(1)              :: trind
  real(kind=8), allocatable :: aux2(:), aux3(:) 

  allocate(aux2(nod2D), aux3(nod3D)) 
  n3=myDim_nod3D+eDim_nod3D             

  if (mype==0) then 

     ! ocean
     ! open files
     filename=trim(ResultPath)//runid//'.'//cyearnew//'.oce.nc'
     status = nf_open(filename, nf_write, ncid)
     if (status .ne. nf_noerr) call handle_err(status)

     ! inquire variable id
     status=nf_inq_varid(ncid, 'time', time_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status=nf_inq_varid(ncid, 'iter', iter_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status=nf_inq_varid(ncid, 'ssh', ssh_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status=nf_inq_varid(ncid, 'u', u_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status=nf_inq_varid(ncid, 'v', v_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status=nf_inq_varid(ncid, 'w', w_varid)
     if (status .ne. nf_noerr) call handle_err(status)
#ifndef use_non_hydrostatic
     status=nf_inq_varid(ncid, 'wpot', wpot_varid)
     if (status .ne. nf_noerr) call handle_err(status)
#endif
     status=nf_inq_varid(ncid, 'temp', tra_varid(1))
     if (status .ne. nf_noerr) call handle_err(status)
     status=nf_inq_varid(ncid, 'salt', tra_varid(2))
     if (status .ne. nf_noerr) call handle_err(status)

     if(use_passive_tracer) then
        do j=1,num_passive_tracer
           write(trind,'(i1)') j
           status = nf_inq_varid(ncid, 'ptr'//trind, tra_varid(index_passive_tracer(j)))
           if (status .ne. nf_noerr) call handle_err(status)
        end do
     end if

     if(use_age_tracer) then
        do j=1,num_age_tracer
           write(trind,'(i1)') j
           status = nf_inq_varid(ncid, 'age'//trind, tra_varid(index_age_tracer(j)))
           if (status .ne. nf_noerr) call handle_err(status)
        end do
     end if
     ! write variables
     ! time and iteration
     status=nf_put_vara_double(ncid, time_varid, save_count_restart, 1, real(seconds_from_yearstart))
     if (status .ne. nf_noerr) call handle_err(status)
     status=nf_put_vara_int(ncid, iter_varid, save_count_restart, 1, istep)
     if (status .ne. nf_noerr) call handle_err(status)
  end if    !! mype==0     

  ! 2d fields
  call gather2D(ssh,aux2)  
  if(mype==0) then            
     start=(/1,save_count_restart/)
     count=(/nod2d, 1/)
     status=nf_put_vara_double(ncid, ssh_varid, start, count, aux2) 
     if (status .ne. nf_noerr) call handle_err(status)
  end if


  ! 3d fields
  call gather3D(uf(1:n3), aux3)   
  if (mype==0) then                  
     start=(/1,save_count_restart/)
     count=(/nod3d, 1/)
     status=nf_put_vara_double(ncid, u_varid, start, count, aux3) 
     if (status .ne. nf_noerr) call handle_err(status)
  end if

  call gather3D(uf(1+n3:2*n3),aux3)  
  if(mype==0) then                      
     status=nf_put_vara_double(ncid, v_varid, start, count, aux3)
     if (status .ne. nf_noerr) call handle_err(status)
  end if

#ifdef use_non_hydrostatic
  call gather3D(uf(1+2*n3:3*n3), aux3)  
  if(mype==0) then                       
     status=nf_put_vara_double(ncid, w_varid, start, count, aux3) 
     if (status .ne. nf_noerr) call handle_err(status)
  end if
#else
  call gather3D(wrhs,aux3)             
  if(mype==0) then                        
     status=nf_put_vara_double(ncid, w_varid, start, count, aux3) 
     if (status .ne. nf_noerr) call handle_err(status)
  end if
  call gather3D(w,aux3)             
  if(mype==0) then                     
     status=nf_put_vara_double(ncid, wpot_varid, start, count, aux3)
     if (status .ne. nf_noerr) call handle_err(status)
  end if
#endif
  do j=1,num_tracer
     call gather3D(tracer(:,j),aux3)    
     if(mype==0) then                      
        status=nf_put_vara_double(ncid, tra_varid(j), start, count, aux3) 
        if (status .ne. nf_noerr) call handle_err(status)
     end if
  end do

  if(mype==0) then
     status=nf_close(ncid)
     if (status .ne. nf_noerr) call handle_err(status)
  end if

  ! ice

#ifdef use_ice
  if(mype==0) then                     
     ! open files
     filename=trim(ResultPath)//runid//'.'//cyearnew//'.ice.nc'
     status = nf_open(filename, nf_write, ncid)
     if (status .ne. nf_noerr) call handle_err(status)

     ! inquire variable id
     status=nf_inq_varid(ncid, 'time', time_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status=nf_inq_varid(ncid, 'iter', iter_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status=nf_inq_varid(ncid, 'area', area_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status=nf_inq_varid(ncid, 'hice', hice_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status=nf_inq_varid(ncid, 'hsnow', hsnow_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status=nf_inq_varid(ncid, 'uice', uice_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status=nf_inq_varid(ncid, 'vice', vice_varid)
     if (status .ne. nf_noerr) call handle_err(status)

     ! write variables

     ! time and iteration
     status=nf_put_vara_double(ncid, time_varid, save_count_restart, 1, real(seconds_from_yearstart))
     if (status .ne. nf_noerr) call handle_err(status)
     status=nf_put_vara_int(ncid, iter_varid, save_count_restart, 1, istep)
     if (status .ne. nf_noerr) call handle_err(status)

     ! 2d fields
     start=(/1,save_count_restart/)
     count=(/nod2d, 1/)
  end if     !! mype=0                      
  call gather2D(a_ice,aux2)                
  if(mype==0) then                           
     status=nf_put_vara_double(ncid, area_varid, start, count, aux2)  
     if (status .ne. nf_noerr) call handle_err(status)
  end if
  call gather2D(m_ice,aux2)              
  if(mype==0) then                          
     status=nf_put_vara_double(ncid, hice_varid, start, count, aux2)   
     if (status .ne. nf_noerr) call handle_err(status)
  end if
  call gather2D(m_snow,aux2)              
  if(mype==0) then                            
     status=nf_put_vara_double(ncid, hsnow_varid, start, count, aux2)  
     if (status .ne. nf_noerr) call handle_err(status)
  end if
  call gather2D(u_ice,aux2)              
  if(mype==0) then                           
     status=nf_put_vara_double(ncid, uice_varid, start, count, aux2)    
     if (status .ne. nf_noerr) call handle_err(status)
  end if
  call gather2D(v_ice,aux2)              
  if(mype==0) then                          
     status=nf_put_vara_double(ncid, vice_varid, start, count, aux2)    
     if (status .ne. nf_noerr) call handle_err(status)

     status=nf_close(ncid)
     if (status .ne. nf_noerr) call handle_err(status)
  end if
#endif

  deallocate(aux3, aux2)
end subroutine write_snapshots
!
!--------------------------------------------------------------------------------------------
!
subroutine write_means_part1
  ! write mean arrays and diagnose variables
  ! SGS parameterizations are saved by write_means_part2
  use o_mesh
  use o_array
  use o_passive_tracer_mod
  use o_age_tracer_mod
  use g_config
  use g_diag
  use g_meanarrays
  use g_parfe
  use g_clock
  implicit none

#include "netcdf.inc" 

  integer                   :: status, ncid, j
  integer                   :: time_varid, iter_varid
  integer                   :: Kv_varid
  integer                   :: start(2), count(2), n3
  character(2000)           :: filename
  character(1)              :: trind
  real(kind=4), allocatable :: aux2(:), aux3(:) 

  n3=myDim_nod3D+eDim_nod3D         
  allocate(aux2(nod2D), aux3(nod3D))  

#ifdef allow_diag

  ! ocean
  if(diag_oce) then
     if(mype==0) then                          
        ! open files
        filename=trim(ResultPath)//runid//'.'//cyearnew//'.oce.diag.nc'
        status = nf_open(filename, nf_write, ncid)
        if (status .ne. nf_noerr) call handle_err(status)

        ! inquire variable id
        status=nf_inq_varid(ncid, 'time', time_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        status=nf_inq_varid(ncid, 'iter', iter_varid)
        if (status .ne. nf_noerr) call handle_err(status)

        if(diag_oce_Kv) then
        status=nf_inq_varid(ncid, 'Kv', Kv_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        end if

        ! write variables

        ! time and iteration
        status=nf_put_vara_double(ncid, time_varid, save_count, 1, seconds_from_yearstart)
        if (status .ne. nf_noerr) call handle_err(status)
        status=nf_put_vara_int(ncid, iter_varid, save_count, 1, istep)
        if (status .ne. nf_noerr) call handle_err(status)

        ! 2d fields
        start=(/1,save_count/)
        count=(/nod2d, 1/)
     end if    !! mype=0                        

     ! 3d fields
     start=(/1,save_count/)
     count=(/nod3d, 1/)

     ! Kv
     if(diag_oce_Kv) then
     call gather3D_real4(Kv,aux3)                  
     if(mype==0) then                          
        status=nf_put_vara_real(ncid, Kv_varid, start, count, aux3)    
        if (status .ne. nf_noerr) call handle_err(status)
        end if
     end if

     if(mype==0) then  
        status=nf_close(ncid)  !close file
        if (status .ne. nf_noerr) call handle_err(status)
     end if
     
  endif  ! diag ocean

#endif

  deallocate(aux3, aux2)
end subroutine write_means_part1
!
!--------------------------------------------------------------------------------------------
!
subroutine write_means_part2
  ! write ocean SGS parameterizations
  use o_param
  use o_mesh
  use o_elements
  use o_array
  use g_config
  use g_diag
  use g_meanarrays
  use g_parfe
  use g_clock
  implicit none

#include "netcdf.inc" 

  integer                   :: m, elem, elnodes(4)
  integer                   :: status, ncid, sgs_varid
  integer                   :: start(2), count(2)
  real(kind=8)              :: array_3d(nod3d)
  character(2000)           :: filename

  ! prepare cluster volume
  ! use wrhs as a temporary array
  wrhs=0.0
  do elem=1,myDim_elem3d                                          
     elnodes=elem3d_nodes(:,elem)
     wrhs(elnodes)=wrhs(elnodes)+voltetra(elem)
  end do

  ! 3d fields
  start=(/1,save_count/)
  count=(/nod3d, 1/)

  if(Redi_GM .and. diag_oce_GM_vel) then
     ! processing
     call process_elem2node(1,array_3d)
     if(mype==0) then
        ! open files
        filename=trim(ResultPath)//runid//'.'//cyearnew//'.oce.diag.nc'
        status = nf_open(filename, nf_write, ncid)
        if (status .ne. nf_noerr) call handle_err(status)
        ! inquire variable id
        status=nf_inq_varid(ncid, 'sgs_u', sgs_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        ! write variables
        status=nf_put_vara_real(ncid, sgs_varid, start, count,real(array_3d,4))
        if (status .ne. nf_noerr) call handle_err(status)
        ! close
        status=nf_close(ncid)
        if (status .ne. nf_noerr) call handle_err(status)
     end if

     ! processing
     call process_elem2node(2,array_3d)
     if(mype==0) then
        ! open files
        filename=trim(ResultPath)//runid//'.'//cyearnew//'.oce.diag.nc'
        status = nf_open(filename, nf_write, ncid)
        if (status .ne. nf_noerr) call handle_err(status)
        ! inquire variable id
        status=nf_inq_varid(ncid, 'sgs_v', sgs_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        ! write variables
        status=nf_put_vara_real(ncid, sgs_varid, start, count,real(array_3d,4))
        if (status .ne. nf_noerr) call handle_err(status)
        ! close
        status=nf_close(ncid)
        if (status .ne. nf_noerr) call handle_err(status)
     end if

  end if

  if(diag_oce_SGS_transp) then
     ! processing
     call process_elem2node(3,array_3d)
     if(mype==0) then
        ! open files
        filename=trim(ResultPath)//runid//'.'//cyearnew//'.oce.diag.nc'
        status = nf_open(filename, nf_write, ncid)
        if (status .ne. nf_noerr) call handle_err(status)
        ! inquire variable id
        status=nf_inq_varid(ncid, 'sgs_ut', sgs_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        ! write variables
        status=nf_put_vara_real(ncid, sgs_varid, start, count,real(array_3d,4))
        if (status .ne. nf_noerr) call handle_err(status)
        ! close
        status=nf_close(ncid)
        if (status .ne. nf_noerr) call handle_err(status)
     end if

     ! processing
     call process_elem2node(4,array_3d)
     if(mype==0) then
        ! open files
        filename=trim(ResultPath)//runid//'.'//cyearnew//'.oce.diag.nc'
        status = nf_open(filename, nf_write, ncid)
        if (status .ne. nf_noerr) call handle_err(status)
        ! inquire variable id
        status=nf_inq_varid(ncid, 'sgs_vt', sgs_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        ! write variables
        status=nf_put_vara_real(ncid, sgs_varid, start, count,real(array_3d,4))
        if (status .ne. nf_noerr) call handle_err(status)
        ! close
        status=nf_close(ncid)
        if (status .ne. nf_noerr) call handle_err(status)
     end if

     ! processing
     call process_elem2node(5,array_3d)
     if(mype==0) then
        ! open files
        filename=trim(ResultPath)//runid//'.'//cyearnew//'.oce.diag.nc'
        status = nf_open(filename, nf_write, ncid)
        if (status .ne. nf_noerr) call handle_err(status)
        ! inquire variable id
        status=nf_inq_varid(ncid, 'sgs_us', sgs_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        ! write variables
        status=nf_put_vara_real(ncid, sgs_varid, start, count,real(array_3d,4))
        if (status .ne. nf_noerr) call handle_err(status)
        ! close
        status=nf_close(ncid)
        if (status .ne. nf_noerr) call handle_err(status)
     end if

     ! processing
     call process_elem2node(6,array_3d)
     if(mype==0) then
        ! open files
        filename=trim(ResultPath)//runid//'.'//cyearnew//'.oce.diag.nc'
        status = nf_open(filename, nf_write, ncid)
        if (status .ne. nf_noerr) call handle_err(status)
        ! inquire variable id
        status=nf_inq_varid(ncid, 'sgs_vs', sgs_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        ! write variables
        status=nf_put_vara_real(ncid, sgs_varid, start, count,real(array_3d,4))
        if (status .ne. nf_noerr) call handle_err(status)
        ! close
        status=nf_close(ncid)
        if (status .ne. nf_noerr) call handle_err(status)
     end if
  end if

end subroutine write_means_part2
!
!--------------------------------------------------------------------------------------------
!
subroutine output()
  use g_config
  use g_clock
  use g_diag
  use g_PARFE
  use g_meanarrays, only : all_meandata
  implicit none
  integer :: i
  logical :: do_output=.false.
  logical :: do_output_snapshot=.false.
  logical*1 :: is_last_day_in_month = .false.
  logical*1 :: should_write_output = .false.

#ifdef allow_calcmeans
  is_last_day_in_month = day_in_month==num_day_in_month(fleapyear,month)
  do i=1, size(all_meandata)
    call should_output(all_meandata(i)%ptr%name//CHAR(0), istep, timenew, day_count, month_count, yearnew, is_last_day_in_month, should_write_output)
    if(should_write_output) then
      call all_meandata(i)%ptr%append_mean_output
    end if
  end do
#endif

  !Mean fields and diagnostics: check whether we want to do output
  call check_output(output_length_unit, do_output)
  if (do_output) then
    if(mype==0) write(*,*)'Writing mean fields (netCDF) ...'
#if defined(allow_calcmeans) || defined(allow_diag)
    call compute_means
    call write_means_part1
#ifdef allow_diag
    if(diag_oce .and. (diag_oce_SGS_transp .or. diag_oce_GM_vel)) call write_means_part2
#endif
    call clean_meanarrays
#endif
    save_count=save_count+1
  end if
  
  call should_output('restart'//CHAR(0), istep, timenew, day_count, month_count, yearnew, is_last_day_in_month, should_write_output)

  !Restart: check whether we want to do output
  call check_output(output_length_unit_restart, do_output)

  if (should_write_output) then
    if(mype==0) write(*,*)'Writing restarts (netCDF) ...'
    call write_snapshots
    save_count_restart=save_count_restart+1
  end if

end subroutine output
!
!--------------------------------------------------------------------------------------------
!
subroutine annual_output(do_output)
  !decides whether it's time to do output
  use g_clock
  implicit none

  logical :: do_output

  if ((daynew == ndpyr) .and. (timenew==86400.)) then
     do_output=.true.
  else
     do_output=.false.
  endif

end subroutine annual_output
!
!--------------------------------------------------------------------------------------------
!
subroutine monthly_output(do_output)
  !decides whether it's time to do output
  use g_clock
  implicit none

  logical :: do_output

  if (day_in_month==num_day_in_month(fleapyear,month) .and. &
       timenew==86400.) then
     do_output=.true.
  else
     do_output=.false.
  end if

end subroutine monthly_output
!
!--------------------------------------------------------------------------------------------
!
subroutine daily_output(do_output)
  !decides whether it's time to do output
  use g_clock
  implicit none

  logical :: do_output

  if (mod(daynew,output_length)==0 .and. timenew==86400.) then
     do_output=.true.
  else
     do_output=.false.
  endif

end subroutine daily_output
!
!--------------------------------------------------------------------------------------------
!
subroutine hourly_output(do_output)
  !decides whether it's time to do output
  use g_clock
  implicit none

  logical :: do_output

  if (mod(timenew, 3600.*output_length)==0) then
     do_output=.true.
  else
     do_output=.false.
  endif

end subroutine hourly_output
!
!--------------------------------------------------------------------------------------------
!
subroutine step_output(do_output)
  !decides whether it's time to do output
  use g_config
  implicit none

  logical :: do_output

  if (mod(istep, output_length)==0) then
     do_output=.true.
  else
     do_output=.false.
  endif

end subroutine step_output
!
!--------------------------------------------------------------------------------------------
!
subroutine check_output(length_unit, do_output)
  use g_config
  implicit none

  character, intent(in)  :: length_unit
  logical,   intent(out) :: do_output

  do_output=.false.
  select case(length_unit)
    case('y')
      call annual_output(do_output)
    case('m')
      call monthly_output(do_output)
    case('d')
      call daily_output(do_output)
    case('h')
      call hourly_output(do_output)
    case('s')
      call step_output(do_output)
    case default			
      write(*,*) 'You did not specify a supported outputflag.'
      write(*,*) 'The program will stop to give you opportunity to do it.'
      call par_ex
      stop
  end select
end subroutine check_output

      
subroutine assert_nf(status)
  implicit none
  integer, intent(in) :: status  
#include "netcdf.inc"

  if(status /= nf_noerr) call handle_err(status)
end subroutine assert_nf


subroutine handle_err(errcode)
  use g_parfe
  implicit none
  
#include "netcdf.inc" 
  
  integer errcode
  
  write(*,*) 'Error: ', nf_strerror(errcode)
  call par_ex
  stop
end subroutine handle_err
!
!--------------------------------------------------------------------------------------------
!
subroutine oce_out
  use o_param
  use o_MESH
  use o_array
  use g_PARFE
  use g_config
  implicit none
  !
  integer                    :: i, j, n3
  real(kind=8), allocatable  :: aux2(:), aux3(:)

  n3=myDim_nod3D+eDim_nod3D              
  allocate(aux2(nod2D), aux3(nod3D))     

  call broadcast2D(ssh,aux2)
  call broadcast3D(uf(1:n3),aux3)
  if (mype==0) then      
     write(*,*) 'writing ocean results (ASCII)'
     open(35,file='ssh.out')
     do i=1,nod3D
        write(35,'(1f9.5)') aux3(i)
        !write(35,'(1e10.3)') aux3(i)
     end do
  end if
  call broadcast3D(uf(1+n3:2*n3), aux3)   
  if(mype==0) then
     do i=1,nod3D
        write(35,'(1f9.5)') aux3(i)
        !write(35,'(1e10.3)') aux3(i)
     end do
     !
     do i=1,nod2d
        write(35,'(1f9.5)') aux2(i)
        !write(35,'(1e10.3)') aux2(i)
     end do
     close(35)
  end if

  if(mype==0) open(36,file='TS.out')
  do j=1,num_tracer
     call broadcast3D(tracer(:,j),aux3)
     if(mype==0) then
        do i=1,nod3D
           write(36,'(1f9.5)') aux3(i)
           !write(36,'(1e12.4)') aux3(i)
        end do
     end if
  end do
  if(mype==0) close(36)

#ifndef use_non_hydrostatic
  call broadcast3D(wrhs,aux3)
  if(mype==0) then
     open(37,file='vvel.out') 
     do i=1,nod3D
        write(37,'(1e10.3)') aux3(i)
     end do
     close(37)
  end if
#else
  call broadcast3D(uf(2*n3+1:3*n3), aux3)
  if(mype==0) then
     open(37,file='vvel.out') 
     do i=1,nod3D
        write(37,'(1e10.3)') aux3(i)
     end do
     close(37)
  end if
#endif

  call broadcast3D(Kv, aux3)
  if(mype==0) then 
     open(38,file='mix_coeff.out')
     do i=1,nod3D
        write(38,'(1e10.3)') aux3(i)
     end do
     close(38)
  end if

  deallocate(aux3, aux2)

end subroutine oce_out
!
!--------------------------------------------------------------------------------------------
!
subroutine ice_out
  use o_MESH
  use i_array
  use g_parfe

  implicit none
  integer :: i
  real(kind=8), allocatable   :: aux2(:) 

  allocate(aux2(nod2D))         

  call broadcast2D(m_ice, aux2) 
  if (mype==0) then 
     write(*,*) 'writing ice results (ASCII)'
     open(35,file='m_ice.out') 
     do i=1,nod2D
        !write(35,'(1f9.5)')  aux2(i)
        write(35,'(1e10.3)') aux2(i)
     end do
     close(35)
  end if

  call broadcast2D(a_ice, aux2)
  if(mype==0) then
     open(36,file='a_ice.out') 
     do i=1,nod2D
        !write(36, '(1f9.5)') aux2(i)
        write(36,'(1e10.3)') aux2(i)
     end do
     close(36)
  end if

  call broadcast2D(m_snow, aux2)
  if(mype==0) then 
     open(37,file='m_snow.out') 
     do i=1,nod2D
        !write(37, '(1f9.5)') aux2(i)
        write(37,'(1e10.3)') aux2(i)
     end do
     close(37)
  end if

  call broadcast2D(u_ice, aux2)
  if(mype==0) then
     open(38,file='u_ice.out')
     do i=1,nod2D
        write(38, '(1f9.5)') aux2(i)
     end do
  end if
  call broadcast2D(v_ice, aux2)
  if(mype==0) then
     do i=1,nod2D
        write(38, '(1f9.5)') aux2(i)
     end do
     close(38)
  end if

  call broadcast2D(net_heat_flux, aux2)
  if(mype==0) then
     open(39,file='heat_water_flux.out')
     do i=1,nod2D
        write(39, *) aux2(i)
     end do
  end if
  call broadcast2D(fresh_wa_flux,aux2)
  if(mype==0) then
     do i=1,nod2D
        write(39, *) aux2(i)
     end do
     close(39)
  end if

  deallocate(aux2)

end subroutine ice_out
!
!--------------------------------------------------------------------------------------------
!
subroutine save_MY_vara
  use o_MESH
  use o_array
  use o_mixing_my2p5_mod
  use g_PARFE
  use g_config
  implicit none

  integer :: i
  real(kind=8), allocatable  :: aux3(:)

  allocate(aux3(nod3D))  

  if (mype==0) then
     write(*,*) 'writing MY2.5 variables (ascII)'
     Write(*,*) 'Notice: Currently MY variables are only saved once at the end of a run.'
     ! Later we may update to save MY in Netcdf formate when required.

     open(35,file=trim(ResultPath)//runid//'_MY_restart.out')
  end if

  call broadcast3D(Kv,aux3)
  if(mype==0) then
     do i=1,nod3D
        write(35,'(1e11.4)') aux3(i)
     end do
  end if

  call broadcast3D(Av,aux3)
  if(mype==0) then
     do i=1,nod3D
        write(35,'(1e11.4)') aux3(i)
     end do
  end if

  call broadcast3D(Kq,aux3)
  if(mype==0) then
     do i=1,nod3D
        write(35,'(1e11.4)') aux3(i)
     end do
  end if

  call broadcast3D(q2,aux3)
  if(mype==0) then
     do i=1,nod3D
        write(35,'(1e11.4)') aux3(i)
     end do
  end if

  call broadcast3D(q2b,aux3)
  if(mype==0) then
     do i=1,nod3D
        write(35,'(1e11.4)') aux3(i)
     end do
  end if

  call broadcast3D(q2l,aux3)
  if(mype==0) then
     do i=1,nod3D
        write(35,'(1e11.4)') aux3(i)
     end do
  end if

  call broadcast3D(q2lb,aux3)
  if(mype==0) then
     do i=1,nod3D
        write(35,'(1e11.4)') aux3(i)
     end do
  end if

  if(mype==0) close(35)

  deallocate(aux3)
end subroutine Save_MY_vara
