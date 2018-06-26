subroutine ocean_init
  !reads the initial state or the restart file for the ocean
  use o_param
  use o_array
  use o_mesh
  use o_passive_tracer_mod
  use o_age_tracer_mod
  use g_config
  use g_clock
  use g_parfe
  implicit none

  integer       :: i, node

  ! ocean dynamic fields and active tracers

  ! read ocean status
  if (.not.r_restart) then
     uf=0.0
     ssh=0.0
     if(mype==0) write(*,*) 'read ocean T&S climate data ', trim(OceClimaDataName)
     call read_init_ts
#ifdef use_cavity
     call init_cavity_ts
#endif
  else
     if(mype==0) write(*,*) 'read ocean restart file'
     call oce_input
     if(mix_scheme=='MY2p5') call read_MY_vara ! this shoulb be separated later!!
  end if
  ! init surface values,
  ! will be updated during iteration according to surface restoring spec.
  do i=1, ToDim_nod2D     
     node=nod3D_below_nod2D(1,i)       
     Tsurf(i)=tracer(node,1)          
     Ssurf(i)=tracer(node,2)          
  end do


  ! ocean passive tracers

  if (use_passive_tracer) then
     if(mype==0) write(*,*) 'initialize passive tracers'
     call passive_tracer_init
     if(ptr_start_year<yearnew .or. &
          (ptr_start_year==yearnew .and. (daynew>1 .or. timenew>dt))) then
        if(mype==0) write(*,*) 'read passive tracer restart fields'
        call passive_tracer_input
     end if
  end if


  ! ocean age tracers

  if (use_age_tracer) then
     if(mype==0) write(*,*) 'initialize age tracers'
     call age_tracer_init
     if(age_tracer_start_year<yearnew .or. &
          (age_tracer_start_year==yearnew .and. (daynew>1 .or. timenew>dt))) then
        call age_tracer_input
     end if
  end if


  ! backup all tracers

  tracer0=tracer

end subroutine ocean_init
!
!----------------------------------------------------------------------------
!
subroutine ocean_init_back
  !A backup for ocean model test cases
  !reads the initial state or the restart file for the ocean
  use o_param
  use o_array
  use g_parfe
  use o_mesh
  use g_config
  use g_clock
  implicit none
  !
  integer                     :: i, eof_flag, num
  real(kind=8)                :: x, y, z, h
  real(kind=8)                :: dtdz, Tb
  real(kind=8)                :: Bu, N, r

  tracer=0.

  dtdz=0.
  tracer(:,1)=10.0
  dtdz=2.0/rho0/2.e-4/5000.
  tracer(:,1)=10.0+dtdz*coord_nod3D(3,:)

  do i=1,myDim_nod2d+eDim_nod2d         
     y=coord_nod2d(2,i)
     x=coord_nod2d(1,i)
     r=sqrt((y)**2+(x+250.e3/r_earth)**2)
     if(r<=25.e3/r_earth) then

        num=num_layers_below_nod2d(i)

        ! !!!!!!!!!!!!!!!!!!!!!! deleted
     end if
  end do

  !Bu=1.0
  !N=Bu*1.0e-4*25.0e3/4500.0
  !!vertical gradient of temperature
  !dtdz=(N**2)/g/2.e-4
  !!background T
  !tracer(:,1)=dtdz*coord_nod3D(3,:)

  do i=1, myDim_nod2D+eDim_nod2D    
     num=nod3D_below_nod2D(1,i)     
     Tsurf(i)=tracer(num,1)         
     Ssurf(i)=tracer(num,2)        
  end do
  tracer0=tracer

  if(mype==0) write(*,*) 'ocean initialization done'
end subroutine ocean_init_back
!
!----------------------------------------------------------------------------
!
subroutine ocean_array_setup
  ! Sets up the arrays needed by the ocean part
  use o_param
  use o_array
  use o_mixing_kpp_mod
  use o_mixing_pp_mod    
  use o_mixing_MY2p5_mod
  use o_mixing_tidal_mod
  use o_mesh
  use o_elements
  use g_config
  use g_parfe
  implicit none

  integer	:: size3D, size2D         

  size3d=ToDim_nod3D
  size2d=ToDim_nod2D

  ! density and pressure arrays
  allocate(density_insitu(size3D), density_ref(size3D)) 
  density_insitu=0.
  density_ref=0.
  allocate(bfsq_3D(size3D), dbsfc_3D(size3D))
  bfsq_3D=0.0
  dbsfc_3D=0.0
  if(grid_type/=2) then
     allocate(hpressure(size3D))
     hpressure=0.
  end if
  if(grid_type/=1) then
     allocate(PGF(2,max_num_layers-1,myDim_elem2D))
     PGF=0.
     call init_pressure_force
  end if

  ! Arrays used for ocean boundary forcing
  allocate(stress_x(size2d))
  allocate(stress_y(size2d)) 
  allocate(heat_flux(size2d)) 
  allocate(water_flux(size2d))
  allocate(Tsurf(size2d))
  allocate(Ssurf(size2d)) 
  allocate(ts_sfc_force(size2d, 2))
  allocate(uv_sfc_force(size2d, 2))
  allocate(uv_bott_force(size2d, 2))
  stress_x=0.
  stress_y=0.
  heat_flux=0.
  water_flux=0.
  ts_sfc_force=0.0
  uv_sfc_force=0.0
  uv_bott_force=0.0

  ! T, S fields, their increments and rhs
  allocate(tracer(size3d,num_tracer))     
  allocate(dtracer(size3d,num_tracer))   
  allocate(tracer_rhs(size3d,num_tracer)) 
  allocate(tracer0(size3d,num_tracer))
  tracer=0.0
  tracer0=0.0
  dtracer=0.0
  tracer_rhs=0.0

  ! u, v fields and their rhs 
#ifndef use_non_hydrostatic
  allocate(uf(2*size3D))             
  allocate(uf0(2*size3D))
  allocate(duf(2*size3D))           
  allocate(uv_rhs(2*size3D))
#else
  allocate(nhp(size3D), nhp0(size3D), nhp_rhs(size3D))   
  allocate(uf(3*size3D))          
  allocate(uf0(3*size3D))
  allocate(duf(3*size3D))           
  allocate(uv_rhs(3*size3D)) 
  nhp=0.
  nhp0=0.
#endif
  uf=0.
  uf0=0.
  duf=0.
  uv_rhs=0.

  !ssh
  allocate(ssh(size2D), ssh0(size2D), dssh(size2d), ssh_rhs(size2D)) 
  ssh=0.
  ssh0=0.
  dssh=0.
  uf=0.
  uf0=0.

  ! arrays for the AB2 coriolis case
  if(.not.use_cori_semi) then
     allocate(ucori(size3d), vcori(size3d))
     allocate(ucori_back(size3D), vcori_back(size3D))
     ucori=0.
     vcori=0.
     ucori_back=0.
     vcori_back=0.
  endif

  ! rhs of w equation and w-potential field
#ifndef use_non_hydrostatic
  allocate(wrhs(size3D),w(size3D))   
  wrhs=0.
  w=0.
#endif

  ! arrays for salt fluxes
  allocate(virtual_salt(size2d), relax_salt(size2d))
  virtual_salt=0.
  relax_salt=0.
#ifdef use_fullfreesurf
  allocate(real_salt_flux(size2d))
  real_salt_flux=0.
#endif  
  if(brine_rejection_param) then
     allocate(salt_brine_rejection(size2d))
     salt_brine_rejection=0.0
  end if

  ! Redi/GM
  allocate(Kh_relative(size3d))
  Kh_relative=0.0
  if (Redi_GM) then    
     if(nslope_version==1) then 
        allocate(neutral_slope(3,max_num_layers-1,myDim_elem2d))
        neutral_slope=0.0
     else
        allocate(neutral_slope_elem(3,myDim_elem3d))
        neutral_slope_elem=0.0
     end if
     allocate(BL_depth(size2d), BL_index(size2d))
     allocate(Sx_neutral_base(myDim_elem2d), Sy_neutral_base(myDim_elem2d))
     BL_depth=0.0
     Sx_neutral_base=0.0
     Sy_neutral_base=0.0
  end if

  ! vertical mixing 
  allocate(Av(ToDim_nod3d))
  Av=0.0
  if(trim(mix_scheme)=='KPP') then
     allocate(Kv(ToDim_nod3d,2))
     Kv=0.0
     call oce_mixing_kpp_init
  else
     allocate(Kv(ToDim_nod3d,1))
     Kv=0.0
  end if
  if(trim(mix_scheme)=='MY2p5') then
     call oce_mixing_MY2p5_init
  end if
  if(trim(mix_scheme)=='PP') then
     call oce_mixing_pp_init
  end if  
  if(tidal_mixing) call oce_mixing_tidal_init 

  ! initialize the fct scheme
#ifdef use_tracer_fct    
  call fct_init		    	   
#endif

  if (use_global_tides) then
     allocate(ssh_gp(size2d))
  end if
  if(mype==0) write(*,*) 'Ocean arrays have been set up'
end subroutine ocean_array_setup
!
!----------------------------------------------------------------------------
!
subroutine set_coriolis_param
  use o_param
  use o_array
  use o_mesh
  use o_elements
  use g_config
  use g_rotate_grid
  use g_parfe
  implicit none

  integer         :: i, j
  real(kind=8)    :: midlat, lon, lat, rlon, rlat

  ! nodal

  allocate(coriolis_param_nod2D(ToDim_nod2d))  
  coriolis_param_nod2D=0.0
  if (.not.rotated_grid) then
     coriolis_param_nod2D=2.0*omega*sin(coord_nod2D(2,:))
  else
     do i=1,ToDim_nod2d                        
        rlat=coord_nod2D(2,i)
        rlon=coord_nod2D(1,i)
        call r2g(lon, lat, rlon, rlat)
        coriolis_param_nod2D(i)=2.0*Omega*sin(lat)
     end do
  end if
  if(fplane) then
     coriolis_param_nod2D=f_fplane 
  end if
  if(betaplane) then    
     midlat=(maxval(coord_nod2d(2,:))+minval(coord_nod2d(2,:)))/2.0 
     coriolis_param_nod2D=f_fplane+beta_betaplane*(coord_nod2D(2,:)-midlat)*r_earth
  end if

  ! elementwise 

  allocate(coriolis_param_elem2D(myDim_elem2D))    
  coriolis_param_elem2D=0.0
  do i=1,myDim_elem2D                                            
     coriolis_param_elem2d(i)=sum(coriolis_param_nod2D(elem2D_nodes(:,i)))/3.0
  end do
end subroutine set_coriolis_param
!
!----------------------------------------------------------------------------
!
subroutine prepare_init_data
  ! this routine is only kept here for a backup, not used by the model
  ! read nc data and save to formatted data 
  use o_PARAM
  use g_config
  use g_parfe
  implicit none

#include "netcdf.inc" 

  integer			:: k, i, j
  integer,parameter             :: num_z=33
  integer			:: itime, latlen, lonlen
  integer			:: status, ncid, varid
  integer			:: lonid, latid
  integer			:: istart(3), icount(3)
  real(kind=8), allocatable	:: lon(:), lat(:)
  real(kind=8), allocatable	:: nc_temp(:,:,:), nc_salt(:,:,:)
  real(kind=8)                  :: dep(num_z)
  character                     :: cdep1*2, cdep2*3, cdep3*4
  character(15)			:: vari
  character(300)               	:: file  

  file=trim(ClimateDataPath)//'Winter_phc2.1_beta.dat.nc'

  ! open file
  status=nf_open(file, nf_nowrite, ncid)
  if (status.ne.nf_noerr)then
     print*,'ERROR: CANNOT READ init_data FILE CORRECTLY !!!!!'
     print*,'Error in opening netcdf file'//file
     call par_ex
     stop
  endif

  ! lat
  status=nf_inq_dimid(ncid, 'lat', latid)
  status=nf_inq_dimlen(ncid, latid, latlen)
  allocate(lat(latlen))
  status=nf_inq_varid(ncid, 'lat', varid)
  status=nf_get_vara_double(ncid,varid,1,latlen,lat)

  ! lon
  status=nf_inq_dimid(ncid, 'lon', lonid)
  status=nf_inq_dimlen(ncid, lonid, lonlen)
  allocate(lon(lonlen))
  status=nf_inq_varid(ncid, 'lon', varid)
  status=nf_get_vara_double(ncid,varid,1,lonlen,lon)

  ! depth
  dep=(/ 0., 10., 20., 30., 50., 75., 100., 125., 150., 200., &
       250., 300., 400., 500., 600., 700., 800., 900., 1000., 1100., &
       1200., 1300., 1400., 1500., 1750., 2000., 2500., 3000., 3500., &
       4000., 4500., 5000., 5500. /)

  ! data
  allocate(nc_temp(lonlen,latlen,33))
  allocate(nc_salt(lonlen,latlen,33))

  do k=1,33
     if(k==1) then
        vari='MAM_Temp_00m'
     elseif(dep(k)<100.) then
        write(cdep1,'(i2)') int(dep((k)))
        vari='MAM_Temp_'//cdep1//'m'
     elseif(dep(k)<1000.) then
        write(cdep2,'(i3)') int(dep((k)))
        vari='MAM_Temp_'//cdep2//'m'
     else
        write(cdep3,'(i4)') int(dep((k)))
        vari='MAM_Temp_'//cdep3//'m'
     endif

     status=nf_inq_varid(ncid, vari, varid)
     if (status.ne.nf_noerr)then
        write(*,*) 'error by getting varid for temp'
        call abort
     end if
     istart = (/1,1,itime/)
     icount= (/lonlen,latlen,1/)
     status=nf_get_vara_double(ncid,varid,istart,icount,nc_temp(:,:,k))
     if (status.ne.nf_noerr)then
        write(*,*) 'error when reading temp'
        call abort
     end if

     if(k==1) then
        vari='MAM_Salt_00m'
     elseif(dep(k)<100.) then
        write(cdep1,'(i2)') int(dep((k)))
        vari='MAM_Salt_'//cdep1//'m'
     elseif(dep(k)<1000.) then
        write(cdep2,'(i3)') int(dep((k)))
        vari='MAM_Salt_'//cdep2//'m'
     else
        write(cdep3,'(i4)') int(dep((k)))
        vari='MAM_Salt_'//cdep3//'m'
     endif

     status=nf_inq_varid(ncid, vari, varid)
     if (status.ne.nf_noerr)then
        write(*,*) 'error by getting varid for salt'
        call abort
     end if
     istart = (/1,1,itime/)
     icount= (/lonlen,latlen,1/)
     status=nf_get_vara_double(ncid,varid,istart,icount,nc_salt(:,:,k))
     if (status.ne.nf_noerr)then
        write(*,*) 'error when reading salt'
        call abort
     end if
  end do

  ! close file
  status=nf_close(ncid)

  ! change dep to negative values
  dep=-dep

!!$  ! save to formated output
!!$  open(36,file='Winter_phc2.1_beta.dat')
!!$!  write(36,*) lonlen, latlen, num_z 
!!$!  write(36,'(1f8.2)') lon
!!$!  write(36,'(1f8.2)') lat
!!$!  write(36,'(1f8.2)') dep
!!$  do i=1, lonlen
!!$     do j=1, latlen
!!$        write(36,'(1f10.5)') nc_temp(i,j,1:num_z)         
!!$     end do
!!$  end do
!!$  do i=1, lonlen
!!$     do j=1, latlen
!!$        write(36,'(1f10.5)') nc_salt(i,j,1:num_z)         
!!$     end do
!!$  end do
!!$  close(36) 


  open (1,file='PHC_t.dat')
  open (2,file='PHC_s.dat')
  write(1, *) -dep
  write(2, *) -dep
  do i=1, lonlen
     do j=1, latlen              
        write(1, *) lon(i), lat(j), (nc_temp(i,j, k), k=1, num_z)
        write(2, *) lon(i), lat(j), (nc_salt(i,j, k), k=1, num_z)
     end do
  end do
  close(1)
  close(2)  


  deallocate(lon, lat, nc_temp, nc_salt)
end subroutine prepare_init_data
!
!-----------------------------------------------------------------------------------
!
