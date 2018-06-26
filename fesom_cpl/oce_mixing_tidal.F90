module o_mixing_tidal_mod
  ! Note: currently only barotropic part is considered.
  ! 
  ! Compute vertical diffusivity and vertical viscosity
  ! deduced from barotropic and baroclinic tidal 
  ! dissipation. For the baroclinic dissipation, we follow
  ! Simmons etal, and for the barotropic dissipation we follow 
  ! Lee etal. Assume Prandtl number unity. 
  !
  ! Reference:
  !
  ! Simmons, Jayne, St. Laurent, and Weaver, 2004:
  ! Tidally driven mixing in a numerical model of the ocean 
  ! general circulation.  Ocean Modelling, vol. 6,
  ! pages 245-263.
  !
  ! Jayne and St. Laurent, 2001:
  ! Parameterizing tidal dissipation over rough topography.
  ! Geophysical Research Letters, vol. 28, pages 811-814.
  !
  ! Hyun-Chul Lee, A. Rosati, and M.J. Spelman, 2006: 
  ! Barotropic tidal mixing effects in a coupled climate model:
  ! ocean conditions in the northern Atlantic
  ! Ocean Modelling, vol 11, pages 464--477
  !
  ! Osborn, T.R., 1980: Estimates of the local rate of vertical diffusion 
  ! from dissipation measurements.  JPO, vol. 10, pages 83-89.
  !
  ! Munk and Anderson, 1948: Notes on a theory of the thermocline. 
  ! Journal of Marine Research, vol 3. pages 276-295.

  use o_param
  use o_array
  use o_mesh
  use o_elements
  use g_rotate_grid
  use g_parfe

  implicit none

  private

  public oce_mixing_tidal 
  public oce_mixing_tidal_init
  private vert_mix_drag 
  private compute_bvfreq

  real, dimension(:), allocatable :: tide_speed      ! tide speed (m/s) from barotropic tide model 
  real, dimension(:), allocatable :: bvfreq          ! buoyancy frequency (sec^-1) 

  real    :: munk_anderson_p          = 0.25     ! (dimensionless) from Munk and Anderson 
  real    :: munk_anderson_sigma      = 3.0      ! (dimensionless) from Munk and Anderson 

  real    :: von_karman               = 0.4
  real    :: bott_drag_cd             = 2.4e-3   ! (dimensionless) bottom drag from Lee etal

  real    :: background_diffusivity   = 0.1e-4   ! (m^2/sec)
  real    :: background_viscosity     = 0.1e-4   ! (m^2/sec)

  real    :: drhodz_min               = 1.e-10   ! (kg/m^4) minimum abs(drhodz) used to compute N^2 
  real    :: tide_speed_min           = 5.e-3    ! (m/s) below which set a mask=0 for drag mixing diffusivity 

  logical :: smooth_bf                = .true.   ! for smoothing N in vertical with 1-2-1
  integer :: num_121_smoothings       = 1        ! number of 1-2-1 smoothing 

contains


  !#######################################################################
  ! Initialization for the oce_mixing_tidal module.
  subroutine oce_mixing_tidal_init 

    integer					:: i, j, fid
    integer					:: num_lon_reg, num_lat_reg
    integer                                     :: p_flag
    real(kind=8)					:: lon, lat
    real(kind=8), allocatable, dimension(:)	:: lon_reg_4u, lat_reg_4u
    real(kind=8), allocatable, dimension(:)	:: lon_reg_4v, lat_reg_4v
    real(kind=8), allocatable, dimension(:)	:: xn2d, yn2d
    real, allocatable, dimension(:,:)		:: amp_reg
    real, allocatable, dimension(:)		:: tide_u_amp, tide_v_amp

    allocate(bvfreq(ToDim_nod3d))
    bvfreq = 0.0  

    allocate(tide_speed(ToDim_nod2D))
    tide_speed = default_tide_speed

    ! read tidal speed (m/s) from a tide model, such as the
    ! Global Inverse Solution TPX07.1 created by OSU.

    if(read_tide_speed) then 

       open(51,file=trim(TideForcingPath)//'lonlat_4U_'//trim(Tmix_tidalmodelname)//'.dat', status='old')
       read(51,*) num_lon_reg, num_lat_reg 
       allocate(lon_reg_4u(num_lon_reg), lat_reg_4u(num_lat_reg))
       do i=1,num_lon_reg
          read(51,*) lon_reg_4u(i)
       end do
       do i=1,num_lat_reg
          read(51,*) lat_reg_4u(i)
       end do
       close(51)
       open(52,file=trim(TideForcingPath)//'lonlat_4V_'//trim(Tmix_tidalmodelname)//'.dat', status='old')
       read(52,*) num_lon_reg, num_lat_reg 
       allocate(lon_reg_4v(num_lon_reg), lat_reg_4v(num_lat_reg))
       do i=1,num_lon_reg
          read(52,*) lon_reg_4v(i)
       end do
       do i=1,num_lat_reg
          read(52,*) lat_reg_4v(i)
       end do
       close(52)

       ! allocate arrays for global data on regular grids
       allocate(amp_reg(num_lon_reg, num_lat_reg)) 
       allocate(tide_u_amp(ToDim_nod2d))
       allocate(tide_v_amp(ToDim_nod2d))
       allocate(xn2d(ToDim_nod2d))
       allocate(yn2d(ToDim_nod2d))

       do i=1, ToDim_nod2d
          if(rotated_grid) then
             call r2g(lon, lat, coord_nod2d(1,i), coord_nod2d(2,i))
             xn2d(i)=lon/rad   ! change unit to degree
             yn2d(i)=lat/rad
          else
             xn2d(i)=coord_nod2d(1,i)/rad
             yn2d(i)=coord_nod2d(2,i)/rad
          end if
          ! change lon range to [0 360]
          if(xn2d(i)<0.) xn2d(i)=xn2d(i) + 360.0  
       end do

       fid=103
       open(fid,file=trim(TideForcingPath)//'U_'//Tmix_tidalconstituent//'_'//trim(Tmix_tidalmodelname)//'.dat', status='old')
       do i=1, num_lon_reg
          do j=1, num_lat_reg
             read(fid, *) amp_reg(i,j)         
          end do
       end do
       close(fid)
       p_flag=0
       call interp_2d_field(num_lon_reg, num_lat_reg, lon_reg_4u, lat_reg_4u, amp_reg, &
            ToDim_nod2d, xn2d, yn2d, tide_u_amp, p_flag)
       tide_u_amp=tide_u_amp/abs(coord_nod3d(3,bt_nds))  ! change transport (m^2/s) to depth mean velocity (m/s)

       open(fid,file=trim(TideForcingPath)//'V_'//Tmix_tidalconstituent//'_'//trim(Tmix_tidalmodelname)//'.dat', status='old')
       do i=1, num_lon_reg
          do j=1, num_lat_reg
             read(fid, *) amp_reg(i,j)         
          end do
       end do
       close(fid)
       p_flag=0
       call interp_2d_field(num_lon_reg, num_lat_reg, lon_reg_4v, lat_reg_4v, amp_reg, &
            ToDim_nod2d, xn2d, yn2d, tide_v_amp, p_flag)
       tide_v_amp=tide_v_amp/abs(coord_nod3d(3,bt_nds))  ! change transport (m^2/s) to depth mean velocity (m/s)

       tide_speed=sqrt(tide_u_amp*tide_u_amp + tide_v_amp*tide_v_amp)

       deallocate(amp_reg, tide_u_amp, tide_v_amp)
       deallocate(xn2d, yn2d)
       deallocate(lat_reg_4v,lon_reg_4v,lat_reg_4u,lon_reg_4u)

    else
       tide_speed=default_tide_speed
    endif

  end subroutine oce_mixing_tidal_init


  !#######################################################################
  !
  ! Compute vertical diffusivity and viscosity
  ! based on one or both of the following dissipation mechanisms:
  !
  ! 1. internal wave breaking as parameterized by Simmons etal.
  !
  ! 2. barotropic tides feeling the bottom drag, as parameterized by 
  !    Lee etal.  
  !
  ! Note: only the second part is coded. Part 1 need to be updated when required.

  subroutine oce_mixing_tidal(viscA, diffK)

    real, dimension(:,:), intent(inout) :: diffK
    real, dimension(:),   intent(inout) :: viscA

    call compute_bvfreq

    !if(use_wave_dissipation)  call vert_mix_wave(viscA, diffK)
    if(use_drag_dissipation)  call vert_mix_drag(viscA, diffK)

    diffK = diffK + background_diffusivity
    viscA = viscA + background_viscosity

  end subroutine oce_mixing_tidal


  !#####################################################################
  !
  ! Compute Brunt-Vaisala (buoyancy) frequency. 
  !
  subroutine compute_bvfreq

    real    :: dz, aux, bf_prev, tmp, fx, dens_up
    integer :: i, k, kn, mr, nodup, nodlo

    fx=g*rho0r
    do i=1,ToDim_nod2d
       kn=num_layers_below_nod2D(i)
       do k=1, kn
          nodup=nod3d_below_nod2D(k,i)
          nodlo=nod3D_below_nod2D(k+1,i)        
          call fcn_density(tracer(nodup,1),tracer(nodup,2),coord_nod3d(3,nodlo),dens_up)
          dz=coord_nod3d(3,nodup)-coord_nod3d(3,nodlo)
          aux=-(dens_up-density_insitu(nodlo)) / dz
	  bvfreq(nodup)=sqrt(fx*max(aux, drhodz_min))
       end do

       if(smooth_bf .and. kn>2) then
          do mr=1,num_121_smoothings
             bf_prev=0.25*bvfreq(nod3d_below_nod2d(1,i))
             do k=2,kn-1
                nodup=nod3d_below_nod2d(k,i)
                nodlo=nod3d_below_nod2d(k+1,i)
                tmp=bvfreq(nodup)
                bvfreq(nodup)=bf_prev+0.5*bvfreq(nodup)+0.25*bvfreq(nodlo)
                bf_prev=0.25*tmp
             end do
          end do
       end if
    end do

  end subroutine compute_bvfreq


  !#####################################################################
  !
  ! Computes tracer diffusivity based on the methods of Lee etal., 
  ! which consider the dissipation from barotropic tides
  ! rubbing against the ocean bottom. 
  !
  ! Assume a unit Prandtl number.
  !
  subroutine vert_mix_drag(viscA, diffK)

    integer :: i, k, kn, nodup, nodlo
    real    :: bottdep, speedr, fx, height, ri, diff_drag
    real, dimension(:,:), intent(inout) :: diffK
    real, dimension(:),   intent(inout) :: viscA

    fx=von_karman/sqrt(bott_drag_cd)

    do i=1,ToDim_nod2d
       if(tide_speed(i)<tide_speed_min) cycle
       kn=num_layers_below_nod2D(i)
       bottdep=-coord_nod3d(3,bt_nds(i))
       speedr=fx/tide_speed(i)

       do k=1, kn
          nodup=nod3d_below_nod2D(k,i)
          nodlo=nod3D_below_nod2D(k+1,i)
          height=bottdep + &
               0.5*(coord_nod3d(3,nodup)+coord_nod3d(3,nodlo))      
          ri = 2.0*(bvfreq(nodup)*height*speedr)**2

          ! compute drag induced diffusivity 
          ! (Lee etal equations (1), (2), and (3))
          diff_drag=max_drag_diffusivity &
               *(1.0 + munk_anderson_sigma*ri)**(-munk_anderson_p)
          diffK(nodup,:)=diffK(nodup,:)+diff_drag
          viscA(nodup)=viscA(nodup)+diff_drag
       enddo
    enddo

  end subroutine vert_mix_drag

end module o_mixing_tidal_mod
