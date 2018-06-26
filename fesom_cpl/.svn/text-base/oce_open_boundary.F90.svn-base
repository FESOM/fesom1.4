!example routines for initializing buffer zone restoring and 
!open boundary velocity restoring.
!----------------------------------------------------------------
subroutine init_restoring_bufferzone
  ! init tracer buffer zone for the Norht East Atlantic
  use o_mesh
  use o_param
  use o_array
  use g_config
  use g_rotate_grid
  use g_parfe
  implicit none
  
  integer                     :: row
  real(kind=8)                :: x, y, dx, dy, dd
  real(kind=8)                :: x1=-20., x2=0.
  real(kind=8)                :: y1=40., y2=55.

  allocate(tracer_restore_coeff(myDim_nod3d+eDim_nod3d))
  tracer_restore_coeff=0.0

  do row=1, myDim_nod3d+eDim_nod3d
     call r2g(x, y, coord_nod3d(1, row), coord_nod3d(2, row))
     x=x/rad
     y=y/rad
     if ((x>x1 .and. x<x2) .and. (y>y1 .and. y<y2)) then
	dx=min(abs(x-x1), abs(x-x2))
	dy=min(abs(y-y1), abs(y-y2))
	dd=min(dx, dy)
	tracer_restore_coeff(row)=min(.25*dd, 1.)*restore_ts_buff
     end if
  end do
  if(mype==0) write(*,*) 'restoring buffer zone ready: North East Atlantic'

end subroutine init_restoring_bufferzone
!
!--------------------------------------------------------------------
!
subroutine init_restoring_bufferzone_SO
  ! init tracer buffer zone for FO001 setup (a southern ocean setup)
  use o_mesh
  use o_param
  use o_array
  use g_config
  use g_rotate_grid
  use g_parfe
  implicit none
  
  integer                     :: row
  real(kind=8)                :: y, d, buffer_dist
  real(kind=8)                :: ymax_local, ymax_global

  buffer_dist=150.0e3  !in m, buffer zone scale

  ymax_local=maxval(coord_nod2d(2,:))

  call MPI_Barrier(MPI_COMM_FESOM,MPIerr)
  call MPI_AllREDUCE(ymax_local, ymax_global, &
       1, MPI_DOUBLE_PRECISION, MPI_MAX, &
       MPI_COMM_FESOM, MPIerr)

  allocate(tracer_restore_coeff(myDim_nod3d+eDim_nod3d))
  tracer_restore_coeff=0.0

  do row=1,myDim_nod3d+eDim_nod3d
     
     y=coord_nod3d(2,row)
     d=(ymax_global-y)*r_earth
     
     if(d<buffer_dist) then
        tracer_restore_coeff(row)=(buffer_dist-d)/buffer_dist*restore_ts_buff
     end if
  end do

  if(mype==0) write(*,*) 'restoring buffer zone ready: FO001'

end subroutine init_restoring_bufferzone_SO
!
!--------------------------------------------------------------------
!
subroutine init_restoring_bufferzone_FOt
  ! init tracer buffer zone for FOt setup
  use o_mesh
  use o_param
  use o_array
  use g_config
  use g_rotate_grid
  use g_parfe
  implicit none
  
  integer                     :: row, n
  real(kind=8)                :: x, y, d, d_min, buffer_dist

  integer, allocatable        :: ind_op(:), ind_op_glo(:)
  real(kind=8), allocatable   :: x_op(:), y_op(:)
  real(kind=8), allocatable   :: x_op_glo(:), y_op_glo(:)

  buffer_dist=60.0e3  !in m, buffer zone scale

  allocate(ind_op(nod2d), x_op(nod2d), y_op(nod2d))
  ind_op=0
  x_op=0.0
  y_op=0.0
  allocate(ind_op_glo(nod2d), x_op_glo(nod2d), y_op_glo(nod2d))
  ind_op_glo=0
  x_op_glo=0.0
  y_op_glo=0.0

  do row=1,myDim_nod2d
     if(index_nod3d(nod3d_below_nod2d(1,row))==12) then
        n=myList_nod2d(row)
        ind_op(n)=1
        x_op(n)=coord_nod2d(1,row)
        y_op(n)=coord_nod2d(2,row)
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

  deallocate(ind_op, x_op, y_op)

  allocate(tracer_restore_coeff(myDim_nod3d+eDim_nod3d))
  tracer_restore_coeff=0.0

  do row=1,myDim_nod3d+eDim_nod3d
     
     x=coord_nod3d(1,row)
     y=coord_nod3d(2,row)
     
     d_min=1000.0e3   !dist. in m
     do n=1,nod2d
        if(ind_op_glo(n)/=1) cycle
        call dist_on_earth(x, y, x_op_glo(n), y_op_glo(n), d)
        d_min=min(d_min, d)
     end do

     if(d_min<buffer_dist) then
        tracer_restore_coeff(row)=(buffer_dist-d_min)/buffer_dist*restore_ts_buff
     end if
  end do

  deallocate(ind_op_glo, x_op_glo, y_op_glo)

  if(mype==0) write(*,*) 'restoring buffer zone ready: FOt'

end subroutine init_restoring_bufferzone_FOt
!
!--------------------------------------------------------------------
!
subroutine init_restoring_bufferzone_Arctic
  ! init T/S buffer zone for Arctic setup
  ! T/S restored to PHC2.1 annual mean.
  use o_mesh
  use o_param
  use o_array
  use g_config
  use g_rotate_grid
  use g_parfe
  implicit none
  !
  integer                     :: i, j, n
  integer                     :: num_lat_reg, num_lon_reg, num_lay_reg
  real(kind=8)                :: pp, pr, tt, ss, lon, lat
  real(kind=8)                :: rest_bound, rest_range
  real(kind=8), external      :: theta
  real(kind=8), allocatable   :: lon_reg(:), lat_reg(:), lay_reg(:)
  real(kind=8), allocatable   :: raw_data(:,:,:)
  real(kind=8), allocatable   :: temp_x(:), temp_y(:)


  !1) T S fields for restoring: T_0, S_0

  ! open global T/S data files
  !open(19,file=trim(ClimateDataPath)//'Winter_phc2.1_beta_ts.out', status='old')
  open(19,file=trim(ClimateDataPath)//'annual_phc_ts.out', status='old')
  ! read reg. grid
  read(19,*) num_lon_reg, num_lat_reg, num_lay_reg
  allocate(lon_reg(num_lon_reg))
  allocate(lat_reg(num_lat_reg))
  allocate(lay_reg(num_lay_reg))
  read(19,*) lon_reg
  read(19,*) lat_reg
  read(19,*) lay_reg
  allocate(raw_data(num_lon_reg,num_lat_reg,num_lay_reg))

  ! model grid coordinates
  allocate(temp_x(myDim_nod3d+eDim_nod3D), temp_y(myDim_nod3d+eDim_nod3D))
  do n=1, myDim_nod3d+eDim_nod3D        
     if(rotated_grid) then
        call r2g(lon, lat, coord_nod3d(1,n), coord_nod3d(2,n))
        temp_x(n)=lon/rad   ! change unit to degree
        temp_y(n)=lat/rad
     else
        temp_x(n)=coord_nod3d(1,n)/rad   
        temp_y(n)=coord_nod3d(2,n)/rad
     end if
     ! change lon range to [0 360]
     if(temp_x(n)<0.) temp_x(n)=temp_x(n) + 360.0  
  end do

  ! read raw data and do interpolation
  do i=1, num_lon_reg
     do j=1, num_lat_reg
        read(19, *) raw_data(i,j,1:num_lay_reg)         
     end do
  end do
  call interp_3d_field(num_lon_reg, num_lat_reg, num_lay_reg, lon_reg, lat_reg, lay_reg, &
       raw_data, nod3d, temp_x, temp_y, coord_nod3d(3,:), tracer0(:,1))

  do i=1, num_lon_reg
     do j=1, num_lat_reg
        read(19, *) raw_data(i,j,1:num_lay_reg)         
     end do
  end do
  call interp_3d_field(num_lon_reg, num_lat_reg, num_lay_reg, lon_reg, lat_reg, lay_reg, &
       raw_data, nod3d, temp_x, temp_y, coord_nod3d(3,:), tracer0(:,2))

  close(19) 

  ! Convert in situ temperature into potential temperature
  pr=0.0_8
  do i=1,myDim_nod3d+eDim_nod3D                       
     tt=tracer0(i,1)
     ss=tracer0(i,2)
     pp=abs(coord_nod3D(3,i))
     tracer0(i,1)=theta(ss, tt, pp, pr)
  end do

  deallocate(temp_y, temp_x, raw_data, lay_reg, lat_reg, lon_reg)


  !2) where to apply restoring for T/S

  allocate(tracer_restore_coeff(myDim_nod3D+eDim_nod3D)) 
  tracer_restore_coeff=0.0
  ! east boundary
  rest_bound=maxval(coord_nod3d(1,:))
  rest_range=1.0*rad
  do i=1, myDim_nod3D+eDim_nod3D         
     if(coord_nod3D(1,i) >= rest_bound-rest_range) then
        tracer_restore_coeff(i)=1.0-abs(rest_bound - coord_nod3D(1,i))/rest_range
     end if
  end do
  ! west boundary
  rest_bound=minval(coord_nod3d(1,:))
  rest_range=3.0*rad
  do i=1, myDim_nod3D+eDim_nod3D        
     if(coord_nod3D(1,i) <= rest_bound+rest_range) then
        tracer_restore_coeff(i)=1.0-abs(rest_bound - coord_nod3D(1,i))/rest_range
     end if
  end do
  tracer_restore_coeff=tracer_restore_coeff*restore_ts_buff

  if(mype==0) write(*,*) 'restoring buffer zone ready'

end subroutine init_restoring_bufferzone_Arctic
!
!------------------------------------------------------------------------------
!
subroutine init_restoring_vel_Arctic
  ! init open boundary velocity for Arctic setup_Arctic
  ! In this version the barotropic velocity is specified (derived from
  ! streamfunction offline)
  use o_mesh
  use o_param
  use o_array
  use g_config
  use g_parfe
  implicit none
  !
  integer                     :: i, el2, elnodes2(3)
  integer                     :: q, row, mn(3)
  real(kind=8)                :: vol, aux
  real(kind=8)                :: tri_u(3), tri_v(3)
  real(kind=8)                :: uu, vv
  real(kind=8), allocatable   :: opbnd_u2d(:), opbnd_v2d(:)
  integer, allocatable        :: map_loc(:)
  ! read barotropic velocity on open boundaries

  if(nmbr_opbnd_tri>0) then   !=======================    
     allocate(opbnd_ssh_rhs(nmbr_opbnd_t2d))
     allocate(map_loc(nod2d))                          
     map_loc=0                                         
     do i=1, nmbr_opbnd_n2D                                
        map_loc(myList_nod2D(opbnd_n2d(i)))=i        
     end do
     allocate(opbnd_u2d(nmbr_opbnd_t2D))
     allocate(opbnd_v2d(nmbr_opbnd_t2D))
     opbnd_u2d=0.
     opbnd_v2d=0.
     open(20,file=trim(OpbndPath)//'Arc_ob_int_vel.out', status='old')
     do i=1,nod2d
        read(20,*) uu, vv         !depth integrated velocity
        if(map_loc(i).ne.0) then                        
           opbnd_u2d(map_loc(i))=uu                    
           opbnd_v2d(map_loc(i))=vv                    
        end if
     end do
     close(20) 
     ! depth averaged velocity
     opbnd_u2d=opbnd_u2d/opbnd_dep
     opbnd_v2d=opbnd_v2d/opbnd_dep

     ! prepare ssh_rhs contribution from open boundary
     ! and correct the result to make sure IN = OUT
     opbnd_ssh_rhs=0.0
     do el2=1, nmbr_opbnd_tri
        elnodes2=opbnd_tri(el2,1:3)
        elnodes2=nod2d_corresp_to_nod3d(elnodes2)
        mn=mapping_opbnd_n2d(elnodes2)
        vol=opbnd_nv(el2,4)/12.0
        tri_u=opbnd_u2d(mn)
        tri_v=opbnd_v2d(mn)
        tri_v=tri_u*opbnd_nv(el2,1)+tri_v*opbnd_nv(el2,2)
        aux=sum(tri_v)
        do q=1, 3
           row=elnodes2(q)
           opbnd_ssh_rhs(mn)=opbnd_ssh_rhs(mn) - (aux+tri_v(q))*vol 
        end do
     end do
     deallocate(opbnd_v2d, opbnd_u2d, map_loc)
  end if                        !=======================    

  ! correction
  aux=0.0
  row=0
  if(nmbr_opbnd_tri>0) then                         
     ! Summation should go over my nodes:            
     do i=1, nmbr_opbnd_t2D     ! all 2D OB nodes    
        q=opbnd_n2D(i)                               
        if(q<=myDim_nod2D) then                      
           aux=aux+opbnd_ssh_rhs(i)                   
           row=row+1                                 !! .... 
        end if
     end do
  end if

  vol=0.0
  q=0
  call MPI_Barrier(MPI_COMM_FESOM,MPIerr)
  call MPI_AllREDUCE( aux, vol, 1, MPI_DOUBLE_PRECISION,MPI_SUM, &
       MPI_COMM_FESOM, MPIerr)
  call MPI_AllREDUCE( row, q, 1, MPI_INTEGER,MPI_SUM, &
       MPI_COMM_FESOM, MPIerr)

  ! vol contains global sum of ssh_rhs, and q the global number 
  ! of open boundary nodes

  if(nmbr_opbnd_tri>0) then                         
     opbnd_ssh_rhs=opbnd_ssh_rhs-vol/real(q)
  end if

  if(mype==0) write(*,*) 'restoring open boundary velocity ready'

end subroutine init_restoring_vel_Arctic
