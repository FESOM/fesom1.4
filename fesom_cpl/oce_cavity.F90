subroutine init_cavity_ts
  use o_mesh
  use o_param
  use o_array
  use g_config
  use g_rotate_grid
  use g_PARFE
  implicit none
  
  integer                     :: i, j, n, n2, n_bd, row, ind
  integer                     :: num_lat_reg, num_lon_reg, num_lay_reg
  integer, allocatable        :: ind_op(:), ind_op_glo(:), nearest_bd_nod(:)
  real(kind=8)                :: pp, pr, aux1,aux2, lon, lat, d, d_min
  real(kind=8)                :: x, y, temp_x, temp_y, temp_z
  real(kind=8), external      :: theta
  real(kind=8), allocatable   :: lon_reg(:), lat_reg(:), lay_reg(:)
  real(kind=8), allocatable   :: raw_data_T(:,:,:), raw_data_S(:,:,:)
  real(kind=8), allocatable   :: x_op(:), y_op(:), dep_op(:)
  real(kind=8), allocatable   :: x_op_glo(:), y_op_glo(:), dep_op_glo(:)

  ! define the cavity boundary line
  allocate(ind_op(nod2d), x_op(nod2d), y_op(nod2d), dep_op(nod2d))
  ind_op=0
  x_op=0.0
  y_op=0.0
  dep_op=0.0
  allocate(ind_op_glo(nod2d), x_op_glo(nod2d), y_op_glo(nod2d), dep_op_glo(nod2d))
  ind_op_glo=0
  x_op_glo=0.0
  y_op_glo=0.0
  dep_op_glo=0.0

  do row=1,myDim_nod2d  ! should not include eDim_nod2d 
     if(cavity_flag_nod2d(row)==0 .and. &
          any(cavity_flag_nod2d(nghbr_nod2d(row)%addresses)==1)) then 
        n=myList_nod2d(row)
        ind_op(n)=1
        x_op(n)=coord_nod2d(1,row)
        y_op(n)=coord_nod2d(2,row)
        if(rotated_grid) then
           call r2g(lon, lat, x_op(n), y_op(n))
           x_op(n)=lon
           y_op(n)=lat
        end if
        dep_op(n)=coord_nod3d(3,bt_nds(row))
     end if
  end do

  call MPI_Barrier(MPI_COMM_FESOM,MPIerr)
  call MPI_AllREDUCE(ind_op, ind_op_glo, &
       nod2d, MPI_INTEGER, MPI_SUM, &
       MPI_COMM_FESOM, MPIerr)
  call MPI_AllREDUCE(x_op, x_op_glo, &
       nod2d, MPI_DOUBLE_PRECISION, MPI_SUM, &
       MPI_COMM_FESOM, MPIerr)
  call MPI_AllREDUCE(y_op, y_op_glo, &
       nod2d, MPI_DOUBLE_PRECISION, MPI_SUM, &
       MPI_COMM_FESOM, MPIerr)
  call MPI_AllREDUCE(dep_op, dep_op_glo, &
       nod2d, MPI_DOUBLE_PRECISION, MPI_SUM, &
       MPI_COMM_FESOM, MPIerr)

  deallocate(ind_op, x_op, y_op, dep_op)

  allocate(nearest_bd_nod(myDim_nod2d+eDim_nod2d))
  do row=1,myDim_nod2d+eDim_nod2d
     if(cavity_flag_nod2d(row)==0) cycle
     
     x=coord_nod2d(1,row)
     y=coord_nod2d(2,row)
     if(rotated_grid) then
        call r2g(lon, lat, x, y)
        x=lon
        y=lat
     end if
          
     d_min=3000.0e3   !dist. in m
     do n=1,nod2d
        if(ind_op_glo(n)/=1) cycle
        call dist_on_earth(x, y, x_op_glo(n), y_op_glo(n), d)
        if(d<d_min) then
           ind=n
           d_min=d
        end if
     end do
     nearest_bd_nod(row)=ind
  end do

  ! open global T/S dataset files
  open(19,file=trim(ClimateDataPath)//trim(OceClimaDataName), status='old')
  ! read reg. grid
  read(19,*) num_lon_reg, num_lat_reg, num_lay_reg
  allocate(lon_reg(num_lon_reg))
  allocate(lat_reg(num_lat_reg))
  allocate(lay_reg(num_lay_reg))
  read(19,*) lon_reg
  read(19,*) lat_reg
  read(19,*) lay_reg
  allocate(raw_data_T(num_lon_reg,num_lat_reg,num_lay_reg))
  allocate(raw_data_S(num_lon_reg,num_lat_reg,num_lay_reg))
  ! read raw data: T
  do i=1, num_lon_reg
     do j=1, num_lat_reg
        read(19, *) raw_data_T(i,j,1:num_lay_reg)         
     end do
  end do
  ! read raw data: S
  do i=1, num_lon_reg
     do j=1, num_lat_reg
        read(19, *) raw_data_S(i,j,1:num_lay_reg)         
     end do
  end do
  close(19) 

  ! interpolate to the nodes under cavity
  do row=1, myDim_nod3d+eDim_nod3D  
     n2=nod2d_corresp_to_nod3d(row)
     if(cavity_flag_nod2d(n2)==0) cycle
     n_bd=nearest_bd_nod(n2)
     temp_z=coord_nod3d(3,row)
     if(temp_z<dep_op_glo(n_bd)) temp_z=dep_op_glo(n_bd)
     temp_x=x_op_glo(n_bd)/rad
     temp_y=y_op_glo(n_bd)/rad
     ! change lon range to [0 360]
     if(temp_x<0.) temp_x=temp_x + 360.0  
     
     call interp_3d_field(num_lon_reg, num_lat_reg, num_lay_reg, &
          lon_reg, lat_reg, lay_reg, &
          raw_data_T, 1, temp_x, temp_y, temp_z, aux1)

     call interp_3d_field(num_lon_reg, num_lat_reg, num_lay_reg, &
          lon_reg, lat_reg, lay_reg, &
          raw_data_S, 1, temp_x, temp_y, temp_z, aux2)
     tracer(row,2)=aux2

     ! Convert in situ temperature into potential temperature
     pr=0.0_8
     pp=abs(temp_z)
     tracer(row,1)=theta(aux2, aux1, pp, pr)
  end do    

  deallocate(ind_op_glo, x_op_glo, y_op_glo, dep_op_glo, nearest_bd_nod)
  deallocate(raw_data_T, raw_data_S, lay_reg, lat_reg, lon_reg)
end subroutine init_cavity_ts
!
!----------------------------------------------------------------------------------------
!
subroutine cavity_heat_water_fluxes
  ! compute the heat and freshwater fluxes under ice cavity
  use o_mesh
  use o_param
  use o_array
  use g_config
  use g_parfe
  implicit none

  integer        :: m, n, row
  real(kind=8)   :: gama, L
  real(kind=8)   :: c2, c3, c4, c5, c6
  real(kind=8)   :: t_i, s_i, p, t_fz

  ! parameter for computing heat and water fluxes
  gama = 1.0e-4     ! heat exchange velocity [m/s]
  L    = 334000.    ! water to ice latent heat [J/Kg], same as used by the ice model

  ! parameter for computing freezing temperature
  c3 = 1.710523e-3
  c4 = -2.154996e-4
  c5 = -0.0575
  c6 = -7.53e-4

  do n=1,myDim_nod2D+eDim_nod2D      
     if(cavity_flag_nod2d(n)==0) cycle   
     row=nod3d_below_nod2d(1,n)
     t_i = tracer(row,1)	
     s_i = tracer(row,2)
     t_fz = c3*(s_i**(3./2.)) + c4*(s_i**2) + c5*s_i + c6*abs(coord_nod3d(3,row))
     !
     heat_flux(n)=vcpw*gama*(t_i - t_fz)  !by Hunter2006 cpw=3974J/Kg is used
     water_flux(n) = -1.0*heat_flux(n)/(L*1000.0)  
  enddo

end subroutine cavity_heat_water_fluxes
!
!----------------------------------------------------------------------------------------
!
