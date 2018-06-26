!==========================================================================
!
subroutine init_tidal_opbnd
  ! initialize tides: read constituents (amplitude and phase)
  ! The global tidal data should be on a grid with longitude range 0-360 degree.
  ! Assumption: the number of columes/rows of the grids for global tidal data
  !             are the same for all constituents and types (U,V, z), though 
  !             the exact lon/lat values can be (are) different for different
  ! 	        types.
  use o_param
  use o_array
  use o_mesh
  use g_config
  use g_rotate_grid
  use g_parfe
  implicit none
  !
  integer					:: i, j, n, num, fid, nodmin,m(1)
  integer					:: num_lon_reg, num_lat_reg
  integer                                       :: p_flag
  character(2)    				:: cons_name
  real(kind=8)					:: lon, lat, dep
  real(kind=8)					:: x, y, x2, y2, d, dmin
  real(kind=8), allocatable, dimension(:)	:: lon_reg_4u, lat_reg_4u
  real(kind=8), allocatable, dimension(:)	:: lon_reg_4v, lat_reg_4v
  real(kind=8), allocatable, dimension(:)	:: lon_reg_4z, lat_reg_4z
  real(kind=8), allocatable, dimension(:,:)	:: amp_reg, pha_reg
  real(kind=8), allocatable, dimension(:)	:: lon_opbnd_n2d, lat_opbnd_n2d

  ! check
  num=len(trim(tidal_constituent))
  if(mod(num,2) /= 0 .or. num/2 /= nmbr_tidal_cons) then
     write(*,*) 'wrong specification of tidal constituents in O_module'
     call par_ex
     stop
  end if

  ! allocate arrays
  if(trim(tide_opbnd_type)=='Flather') then
     allocate(opbnd_u_tide(nmbr_opbnd_t2D))
     allocate(opbnd_v_tide(nmbr_opbnd_t2D))
     allocate(opbnd_z_tide(nmbr_opbnd_t2D))
     opbnd_u_tide=0.
     opbnd_v_tide=0.
     opbnd_z_tide=0.
  elseif(trim(tide_opbnd_type)=='ssh') then
     allocate(opbnd_z_tide(nmbr_opbnd_t2D))
     allocate(opbnd_z0_tide(nmbr_opbnd_t2D))
     opbnd_z_tide=0.  
     opbnd_z0_tide=0.
  else !'vel'
     allocate(opbnd_u_tide(nmbr_opbnd_t2D))
     allocate(opbnd_v_tide(nmbr_opbnd_t2D))
     opbnd_u_tide=0.
     opbnd_v_tide=0.
  end if

  ! allocate tidal amp and pha arrays
  allocate(tide_period_coeff(nmbr_tidal_cons))
  tide_period_coeff=0.
  if(trim(tide_opbnd_type)/='ssh') then
     allocate(tide_u_amp(nmbr_opbnd_t2d, nmbr_tidal_cons))
     allocate(tide_v_amp(nmbr_opbnd_t2d, nmbr_tidal_cons))
     allocate(tide_u_pha(nmbr_opbnd_t2d, nmbr_tidal_cons))
     allocate(tide_v_pha(nmbr_opbnd_t2d, nmbr_tidal_cons))
     tide_u_amp=0.
     tide_v_amp=0.
     tide_u_pha=0.
     tide_v_pha=0.
  end if
  if(trim(tide_opbnd_type)/='vel') then
     allocate(tide_z_amp(nmbr_opbnd_t2d, nmbr_tidal_cons))
     allocate(tide_z_pha(nmbr_opbnd_t2d, nmbr_tidal_cons))
     tide_z_amp=0.
     tide_z_pha=0.
  end if

  ! read the regular grid (for velocity)
  if(trim(tide_opbnd_type)/='ssh') then
     open(51,file=trim(TideForcingPath)//'lonlat_4U_'//trim(tidemodelname)//'.dat', status='old')
     read(51,*) num_lon_reg, num_lat_reg 
     allocate(lon_reg_4u(num_lon_reg), lat_reg_4u(num_lat_reg))
     do i=1,num_lon_reg
        read(51,*) lon_reg_4u(i)
     end do
     do i=1,num_lat_reg
        read(51,*) lat_reg_4u(i)
     end do
     close(51)
     open(52,file=trim(TideForcingPath)//'lonlat_4V_'//trim(tidemodelname)//'.dat', status='old')
     read(52,*) num_lon_reg, num_lat_reg 
     allocate(lon_reg_4v(num_lon_reg), lat_reg_4v(num_lat_reg))
     do i=1,num_lon_reg
        read(52,*) lon_reg_4v(i)
     end do
     do i=1,num_lat_reg
        read(52,*) lat_reg_4v(i)
     end do
     close(52)
  end if

  ! read the regular grid (for ssh)
  if(trim(tide_opbnd_type)/='vel') then
     open(103,file=trim(TideForcingPath)//'lonlat_4z_'//trim(tidemodelname)//'.dat', status='old')
     read(103,*) num_lon_reg, num_lat_reg 
     allocate(lon_reg_4z(num_lon_reg), lat_reg_4z(num_lat_reg))
     do i=1,num_lon_reg
        read(103,*) lon_reg_4z(i)
     end do
     do i=1,num_lat_reg
        read(103,*) lat_reg_4z(i)
     end do
     close(103)
  end if

  ! allocate arrays for global data on regular grids
  allocate(amp_reg(num_lon_reg, num_lat_reg), pha_reg(num_lon_reg, num_lat_reg)) 

  ! open-boundary nodes coordinates
  allocate(lon_opbnd_n2d(nmbr_opbnd_n2d), lat_opbnd_n2d(nmbr_opbnd_n2d))
  do i=1, nmbr_opbnd_n2d
     n=opbnd_n2d(i)
     if(rotated_grid) then
        call r2g(lon, lat, coord_nod2d(1,n), coord_nod2d(2,n))
        lon_opbnd_n2d(i)=lon/rad   ! change unit to degree
        lat_opbnd_n2d(i)=lat/rad
     else
        lon_opbnd_n2d(i)=coord_nod2d(1,n)/rad
        lat_opbnd_n2d(i)=coord_nod2d(2,n)/rad
     end if
     ! change lon range to [0 360]
     if(lon_opbnd_n2d(i)<0.) lon_opbnd_n2d(i)=lon_opbnd_n2d(i) + 360.0  
  end do

  ! read and interpolate each constituent (for U,V,z and their phase) 
  do n=1,nmbr_tidal_cons
     cons_name=tidal_constituent(2*n-1:2*n)
     call tide_period(cons_name, tide_period_coeff(n))

     if(trim(tide_opbnd_type)/='ssh') then
        fid=103+n
        open(fid,file=trim(TideForcingPath)//'U_'//cons_name//'_'//trim(tidemodelname)//'.dat', status='old')
        do i=1, num_lon_reg
           do j=1, num_lat_reg
              read(fid, *) amp_reg(i,j)         
           end do
        end do
        do i=1, num_lon_reg
           do j=1, num_lat_reg
              read(fid, *) pha_reg(i,j)         
           end do
        end do
        close(fid)

        p_flag=0
        call interp_2d_field(num_lon_reg, num_lat_reg, lon_reg_4u, lat_reg_4u, amp_reg, &
             nmbr_opbnd_n2D, lon_opbnd_n2d, lat_opbnd_n2d, tide_u_amp(1:nmbr_opbnd_n2d,n),p_flag)
        tide_u_amp(:,n)=tide_u_amp(:,n)/opbnd_dep  ! change transport (m^2/s) to depth mean velocity (m/s)

        p_flag=1
        call interp_2d_field(num_lon_reg, num_lat_reg, lon_reg_4u, lat_reg_4u, pha_reg, &
             nmbr_opbnd_n2D, lon_opbnd_n2d, lat_opbnd_n2d, tide_u_pha(1:nmbr_opbnd_n2d,n),p_flag)
        tide_u_pha(:,n)=tide_u_pha(:,n)*rad   ! change units of phase from degree to radian 


        open(fid,file=trim(TideForcingPath)//'V_'//cons_name//'_'//trim(tidemodelname)//'.dat', status='old')
        do i=1, num_lon_reg
           do j=1, num_lat_reg
              read(fid, *) amp_reg(i,j)         
           end do
        end do
        do i=1, num_lon_reg
           do j=1, num_lat_reg
              read(fid, *) pha_reg(i,j)         
           end do
        end do
        close(fid)

        p_flag=0
        call interp_2d_field(num_lon_reg, num_lat_reg, lon_reg_4v, lat_reg_4v, amp_reg, &
             nmbr_opbnd_n2D, lon_opbnd_n2d, lat_opbnd_n2d, tide_v_amp(1:nmbr_opbnd_n2d,n),p_flag)
        tide_v_amp(:,n)=tide_v_amp(:,n)/opbnd_dep  ! change transport (m^2/s) to depth mean velocity (m/s)

        p_flag=1
        call interp_2d_field(num_lon_reg, num_lat_reg, lon_reg_4v, lat_reg_4v, pha_reg, &
             nmbr_opbnd_n2D, lon_opbnd_n2d, lat_opbnd_n2d, tide_v_pha(1:nmbr_opbnd_n2d,n),p_flag)
        tide_v_pha(:,n)=tide_v_pha(:,n)*rad   ! change units of phase from degree to radian 

     end if


     if(trim(tide_opbnd_type)/='vel') then
        open(fid,file=trim(TideForcingPath)//'z_'//cons_name//'_'//trim(tidemodelname)//'.dat', status='old')
        do i=1, num_lon_reg
           do j=1, num_lat_reg
              read(fid, *) amp_reg(i,j)         
           end do
        end do
        do i=1, num_lon_reg
           do j=1, num_lat_reg
              read(fid, *) pha_reg(i,j)         
           end do
        end do
        close(fid)

        p_flag=0
        call interp_2d_field(num_lon_reg, num_lat_reg, lon_reg_4z, lat_reg_4z, amp_reg, &
             nmbr_opbnd_n2D, lon_opbnd_n2d, lat_opbnd_n2d, tide_z_amp(1:nmbr_opbnd_n2d,n),p_flag)

        p_flag=1
        call interp_2d_field(num_lon_reg, num_lat_reg, lon_reg_4z, lat_reg_4z, pha_reg, &
             nmbr_opbnd_n2D, lon_opbnd_n2d, lat_opbnd_n2d, tide_z_pha(1:nmbr_opbnd_n2d,n),p_flag)
        tide_z_pha(:,n)=tide_z_pha(:,n)*rad   ! change units of phase from degree to radian 
     end if
  end do

  ! amplifying the magnitude for process studies
  if(trim(tide_opbnd_type)/='ssh') then
     tide_u_amp=tide_u_amp*tide_amplify_coeff
     tide_v_amp=tide_v_amp*tide_amplify_coeff
  end if
  if(trim(tide_opbnd_type)/='vel') then
     tide_z_amp=tide_z_amp*tide_amplify_coeff
  end if


  ! deallocate temporary arrays
  deallocate(lat_opbnd_n2d, lon_opbnd_n2d)
  deallocate(pha_reg, amp_reg)
  if(trim(tide_opbnd_type)/='ssh') deallocate(lat_reg_4v,lon_reg_4v,lat_reg_4u,lon_reg_4u)
  if(trim(tide_opbnd_type)/='vel') deallocate(lat_reg_4z,lon_reg_4z)

  if(mype==0) write(*,*) 'Tidal constituents ', trim(tidal_constituent), ' have been loaded for opbnd'
end subroutine init_tidal_opbnd
!
!==========================================================================
!
subroutine update_tidal_opbnd
  use o_param
  use o_array
  use o_mesh
  use g_config
  implicit none
  integer		:: n
  real(kind=8)		:: aux
  !
  aux=2.0*pi*istep*dt
  if(trim(tide_opbnd_type)=='Flather') then
     opbnd_u_tide=0.
     opbnd_v_tide=0.
     opbnd_z_tide=0.
     do n=1, nmbr_tidal_cons
        opbnd_u_tide=opbnd_u_tide + tide_u_amp(:,n)*cos(aux/tide_period_coeff(n)-tide_u_pha(:,n))
        opbnd_v_tide=opbnd_v_tide + tide_v_amp(:,n)*cos(aux/tide_period_coeff(n)-tide_v_pha(:,n))
        opbnd_z_tide=opbnd_z_tide + tide_z_amp(:,n)*cos(aux/tide_period_coeff(n)-tide_z_pha(:,n))
     end do
  elseif(trim(tide_opbnd_type)=='ssh') then
     opbnd_z0_tide=opbnd_z_tide
     opbnd_z_tide=0.
     do n=1, nmbr_tidal_cons
        opbnd_z_tide=opbnd_z_tide + tide_z_amp(:,n)*cos(aux/tide_period_coeff(n)-tide_z_pha(:,n))
     end do
  else
     opbnd_u_tide=0.
     opbnd_v_tide=0.
     do n=1, nmbr_tidal_cons
        opbnd_u_tide=opbnd_u_tide + tide_u_amp(:,n)*cos(aux/tide_period_coeff(n)-tide_u_pha(:,n))
        opbnd_v_tide=opbnd_v_tide + tide_v_amp(:,n)*cos(aux/tide_period_coeff(n)-tide_v_pha(:,n))
     end do
  end if
  !
end subroutine update_tidal_opbnd
!
!==========================================================================
!
subroutine tide_period(cons_name, period)
  use g_parfe
  use g_config
  implicit none
  real(kind=8)		:: period
  real(kind=8)          :: sec_p_h
  character(2) 		:: cons_name
  !  
  sec_p_h=3600.0_8
  if(cons_name=='M2') then
     period=12.420601*sec_p_h
  elseif(cons_name=='S2') then
     period=12.0*sec_p_h
  elseif(cons_name=='N2') then
     period=12.658348*sec_p_h 
  elseif(cons_name=='K2') then
     period=11.967235*sec_p_h  
  elseif(cons_name=='K1') then
     period=23.93447*sec_p_h 
  elseif(cons_name=='O1') then
     period=25.819342*sec_p_h
  elseif(cons_name=='P1') then
     period=24.06589*sec_p_h
  elseif(cons_name=='Q1') then
     period=26.868357*sec_p_h
  elseif(cons_name=='S1') then
     period=24.0*sec_p_h
  elseif(cons_name=='Mf') then
     period=13.661*24.0*sec_p_h
  elseif(cons_name=='Mm') then
     period=27.555*24.0*sec_p_h
  else
     write(*,*) 'Error! One or more tidal constituents you selected are not specified in the code.'
     call par_ex
     stop
  endif
end subroutine tide_period
!
