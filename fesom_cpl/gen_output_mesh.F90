subroutine init_output_mesh
  ! Initialize output files for mesh diagnose:
  ! cluster volume
  !
  ! Coded by Qiang Wang
  ! Reviewed by ??
  !-----------------------------------------------------------  

  use o_mesh
  use g_config
  use g_clock
  use g_PARFE
  implicit none

#include "netcdf.inc" 

  integer                   :: status, ncid, dimid_rec, j
  integer                   :: dimid_3d, dimids(2)
  integer                   :: time_varid, iter_varid
  integer                   :: vol_varid
  character(100)            :: longname
  character(2000)            :: filename

  if(yearnew==yearold) return
  if (mype/=0) return

  filename=trim(ResultPath)//runid//'.'//cyearnew//'.mesh.diag.nc'

  ! create a file
  status = nf_create(filename, nf_clobber, ncid)
  if (status.ne.nf_noerr) call handle_err(status)

  ! Define the dimensions
  ! status = nf_def_dim(ncid, 'nodes_2d', nod2d, dimid_2d)
  ! if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_dim(ncid, 'nodes_3d', nod3d, dimid_3d)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_dim(ncid, 'T', NF_UNLIMITED, dimid_rec)
  if (status .ne. nf_noerr) call handle_err(status)

  ! Define the time and iteration variables
  status = nf_def_var(ncid, 'time', NF_DOUBLE, 1, dimid_rec, time_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_var(ncid, 'iter', NF_INT, 1, dimid_rec, iter_varid)
  if (status .ne. nf_noerr) call handle_err(status)

  ! Define the netCDF variables 
  ! In Fortran, the unlimited dimension must come
  ! last on the list of dimids.
  dimids(1) = dimid_3d
  dimids(2) = dimid_rec

  status = nf_def_var(ncid, 'cluster_vol', NF_FLOAT, 2, dimids, vol_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  
  ! Assign long_name and units attributes to variables.
   longname='3D mesh cluster volume'
  status = nf_put_att_text(ncid, vol_varid, 'description', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, vol_varid, 'units', 3, 'm^3')
  if (status .ne. nf_noerr) call handle_err(status)

  status = nf_enddef(ncid)
  if (status .ne. nf_noerr) call handle_err(status)

  status=nf_close(ncid)
  if (status .ne. nf_noerr) call handle_err(status)

end subroutine init_output_mesh
!
!-----------------------------------------------------------------------------------
!
subroutine write_mesh_diag
  ! write updated cluster volume
  !
  ! Coded by Qiang Wang
  ! Reviewed by ??
  !--------------------------------------------------------------  

  use o_mesh
  use o_elements
  use g_config
  use g_parfe
  use g_clock
  implicit none

#include "netcdf.inc" 

  integer                   :: status, ncid, j
  integer                   :: time_varid, iter_varid
  integer                   :: vol_varid
  integer                   :: start(2), count(2)
  real(kind=8)              :: sec_in_year
  character(2000)            :: filename
  real(kind=8), allocatable :: aux3(:) 
      
  allocate(aux3(nod3D))  

  sec_in_year=dt*istep

  if (mype==0) then 

     ! ocean

     ! open files
     filename=trim(ResultPath)//runid//'.'//cyearnew//'.mesh.diag.nc'
     status = nf_open(filename, nf_write, ncid)
     if (status .ne. nf_noerr) call handle_err(status)

     ! inquire variable id
     status=nf_inq_varid(ncid, 'time', time_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status=nf_inq_varid(ncid, 'iter', iter_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status=nf_inq_varid(ncid, 'cluster_vol', vol_varid)
     if (status .ne. nf_noerr) call handle_err(status)

     ! write variables

     ! time and iteration
     status=nf_put_vara_double(ncid, time_varid, save_count, 1, sec_in_year)
     if (status .ne. nf_noerr) call handle_err(status)
     status=nf_put_vara_int(ncid, iter_varid, save_count, 1, istep)
     if (status .ne. nf_noerr) call handle_err(status)

     ! 3d fields
     start=(/1,save_count/)
     count=(/nod3d, 1/)
  end if    !! mype=0        

  call broadcast3D(cluster_vol_3D,aux3)       
  if(mype==0) then                      
     status=nf_put_vara_real(ncid, vol_varid, start, count, real(aux3,4))  
     if (status .ne. nf_noerr) call handle_err(status)
  end if

  if(mype==0) then                       
     status=nf_close(ncid)
     if (status .ne. nf_noerr) call handle_err(status)
  end if

  deallocate(aux3)
end subroutine write_mesh_diag
!
!-------------------------------------------------------------------------
!
subroutine write_initial_mesh_diag
  ! Write initial mesh information:
  ! cluster area and volume
  !
  ! Coded by Qiang Wang
  ! Reviewed by ??
  !-----------------------------------------------------------  

  use o_mesh
  use o_elements
  use g_config
  use g_clock
  use g_PARFE
  implicit none

#include "netcdf.inc" 

  integer                   :: status, ncid, j
  integer                   :: dimid_2d, dimid_3d, dimids
  integer                   :: time_varid, iter_varid
  integer                   :: area_varid, vol_varid
  integer                   :: start, count
  character(100)            :: longname
  character(2000)            :: filename
  real(kind=8), allocatable :: aux2(:), aux3(:) 

  if (mype==0) then  ! create a file

     filename=trim(ResultPath)//runid//'.initial.mesh.diag.nc'

     ! create a file
     status = nf_create(filename, nf_clobber, ncid)
     if (status.ne.nf_noerr) call handle_err(status)

     ! Define the dimensions
     status = nf_def_dim(ncid, 'nodes_2d', nod2d, dimid_2d)
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_def_dim(ncid, 'nodes_3d', nod3d, dimid_3d)
     if (status .ne. nf_noerr) call handle_err(status)
  
     ! Define the netCDF variables for 2D 
     ! In Fortran, the unlimited dimension must come
     ! last on the list of dimids.
     dimids = dimid_2d
     
     status = nf_def_var(ncid, 'cluster_area', NF_FLOAT, 1, dimids, area_varid)
     if (status .ne. nf_noerr) call handle_err(status)

     ! Define the netCDF variables for 3D
     dimids = dimid_3d
          
     status = nf_def_var(ncid, 'cluster_vol', NF_FLOAT, 1, dimids, vol_varid)
     if (status .ne. nf_noerr) call handle_err(status)
  
     ! Assign long_name and units attributes to variables.
     longname='2D mesh cluster area'
     status = nf_put_att_text(ncid, area_varid, 'description', len_trim(longname), trim(longname)) 
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_put_att_text(ncid, area_varid, 'units', 3, 'm^2')
     if (status .ne. nf_noerr) call handle_err(status)
     longname='3D mesh cluster volume: initial'
     status = nf_put_att_text(ncid, vol_varid, 'description', len_trim(longname), trim(longname)) 
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_put_att_text(ncid, vol_varid, 'units', 3, 'm^3')
     if (status .ne. nf_noerr) call handle_err(status)
     
     status = nf_enddef(ncid)
     if (status .ne. nf_noerr) call handle_err(status)
     
     status=nf_close(ncid)
     if (status .ne. nf_noerr) call handle_err(status)

  end if  !file created

  !-------------------------------------------------------------------

  allocate(aux2(nod2d), aux3(nod3D))  

  if (mype==0) then 

     ! open files
     filename=trim(ResultPath)//runid//'.initial.mesh.diag.nc'
     status = nf_open(filename, nf_write, ncid)
     if (status .ne. nf_noerr) call handle_err(status)

     ! inquire variable id
     status=nf_inq_varid(ncid, 'cluster_area', area_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status=nf_inq_varid(ncid, 'cluster_vol', vol_varid)
     if (status .ne. nf_noerr) call handle_err(status)
  end if    !! mype=0        

  ! write variables

  call broadcast2D(cluster_area_2D,aux2)       
  if(mype==0) then                      
     ! 2d fields
     start=1
     count=nod2d
     status=nf_put_vara_real(ncid, area_varid, start, count, real(aux2,4))  
     if (status .ne. nf_noerr) call handle_err(status)
  end if

  call broadcast3D(cluster_vol_3D,aux3)       
  if(mype==0) then          
     ! 3d fields
     start=1
     count=nod3d 
     status=nf_put_vara_real(ncid, vol_varid, start, count, real(aux3,4))  
     if (status .ne. nf_noerr) call handle_err(status)
  end if

  if(mype==0) then                       
     status=nf_close(ncid)
     if (status .ne. nf_noerr) call handle_err(status)
  end if

  deallocate(aux2, aux3)

end subroutine write_initial_mesh_diag

