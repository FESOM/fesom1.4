!example routines for initializing and/or restoring the
!tracer fields in process studies (e.g., overflows)
!--------------------------------------------------------------------


subroutine init_ts_CavH1
  ! init T/S source for cavity test case of Hunter 1NN
  use o_mesh
  use o_param
  use o_array
  use g_config
  use g_rotate_grid
  use g_parfe
  implicit none

  tracer(:,1)=-1.9
  tracer(:,2)=34.4
  
  tracer0=tracer  
end subroutine init_ts_CavH1
!
!--------------------------------------------------------------------
!
subroutine init_tracers_FOtide
  ! init T/S source for FOtide
  use o_mesh
  use o_param
  use o_array
  use g_config
  use g_rotate_grid
  use g_parfe
  implicit none

  integer                     :: row, n, ind, fileID
  integer, allocatable        :: index_source(:), temp_arr2d(:)
  real(kind=8)                :: z
  real(kind=8)                :: zn(5), tn(5), sn(5)
  real(kind=8)                :: ts(4), ss(4), tint(4), sint(4)
  character*300               :: file_name

  data zn /0.0, -300.0, -550.0, -2000.0, -4000.0/
  data tn /-1.9, -1.75, 0.7, 0.0, -0.3/
  data sn /34.3, 34.42, 34.68, 34.66, 34.65/

  ts=(tn(1:4)-tn(2:5))/(zn(1:4)-zn(2:5))
  ss=(sn(1:4)-sn(2:5))/(zn(1:4)-zn(2:5))
  tint=tn(1:4)-ts*zn(1:4)
  sint=sn(1:4)-ss*zn(1:4)

  !background
  do row=1,myDim_nod3d+eDim_nod3d
     z=coord_nod3d(3,row)

     if(z>zn(2)) then
        tracer(row,1)=ts(1)*z+tint(1)
        tracer(row,2)=ss(1)*z+sint(1)        
     else if(z>zn(3)) then
        tracer(row,1)=ts(2)*z+tint(2)
        tracer(row,2)=ss(2)*z+sint(2)  
     else if(z>zn(4)) then
        tracer(row,1)=ts(3)*z+tint(3)
        tracer(row,2)=ss(3)*z+sint(3)    
     else if(z>zn(5)) then
        tracer(row,1)=ts(4)*z+tint(4)
        tracer(row,2)=ss(4)*z+sint(4)     
     else
        tracer(row,1)=tn(5)
        tracer(row,2)=sn(5)
     end if
     tracer(row,3)=0.0
  end do

  !source water (ISW)

  !read in info. specifying the initial source region 
  allocate(index_source(myDim_nod2D+eDim_nod2D))	 			

  file_name=trim(meshpath)//'index_source_region.out' 
  fileID=150
  open(fileID, file=file_name)

  allocate(temp_arr2d(nod2d))
  temp_arr2d=0
  do n=1, myDim_nod2D+eDim_nod2D
     temp_arr2d(myList_nod2D(n))=n
  end do

  do n=1,nod2D
     read(fileID,*) ind
     if (temp_arr2d(n)>0) then
        index_source(temp_arr2d(n))=ind
     end if
  end do
  close(fileID)

  ! set source values
  do n=1, myDim_nod3D+eDim_nod3D
     row=nod2d_corresp_to_nod3d(n)
     if(index_source(row)>0 .and. coord_nod3d(3,n)<=-400.0) then
        tracer(n,1)=-2.2
        tracer(n,2)=34.6
        tracer(n,3)=1.0
     end if
  end do

  tracer0=tracer  
 
  deallocate(index_source, temp_arr2d)

  if(mype==0) write(*,*) 'ambient and source water prescribed: FOtide'

end subroutine init_tracers_FOtide
!
!--------------------------------------------------------------------
!
subroutine restore_source_FOtide
  ! restore T/S source for FOtide setup
  use o_mesh
  use o_param
  use o_array
  use g_config
  use g_rotate_grid
  use g_parfe
  implicit none
 
  integer                     :: row
  real(kind=8)                :: x, y

  ! restoring for tracers

  do row=1,myDim_nod3d+eDim_nod3d
     x=coord_nod3d(1,row)
     y=coord_nod3d(2,row)
 
     if(y<=-76.0 .and. x>-42.0 .and. x<-31.0) then   
        tracer(row,:)=tracer0(row,:)
     end if

  end do
  
end subroutine restore_source_FOtide
!
!=============================================================================
!
subroutine init_source_RS1
  ! init T/S source for RS1 setup
  use o_mesh
  use o_param
  use o_array
  use g_config
  use g_rotate_grid
  use g_parfe
  implicit none

  integer                     :: i, k, n, n2, num_hor_reg, num_ver_reg, num_sn
  integer, allocatable        :: source_nodes(:)
  real(kind=8)                :: x, y, z, rx, ry
  real(kind=8), allocatable   :: ver_reg(:), hor_reg(:), raw_t(:,:), raw_s(:,:)


  !tracer fields for restoring

  !passive tracer
  tracer(:,3)=0.0

  !source profile
  open(39,file=trim(MeshPath)//'RS_init_special.out', status='old')
  read(39,*) num_hor_reg, num_ver_reg
  allocate(hor_reg(num_hor_reg),ver_reg(num_ver_reg))
  read(39,*) hor_reg  
  read(39,*) ver_reg
  allocate(raw_t(num_hor_reg,num_ver_reg))
  allocate(raw_s(num_hor_reg,num_ver_reg))
  do i=1,num_hor_reg
     read(39, *) raw_t(i,1:num_ver_reg)
  end do
  do i=1,num_hor_reg
     read(39, *) raw_s(i,1:num_ver_reg)
  end do
  close(39)

  !2d nodes in the source region
  open(40,file=trim(MeshPath)//'source_2d_nodes.out', status='old')
  read(40,*) num_sn
  allocate(source_nodes(num_sn))
  read(40,*) source_nodes 
  close(40)

  !interpolate
  do n2=1,myDim_nod2d+eDim_nod2d
     if(any(source_nodes==myList_nod2d(n2))) then
	if(rotated_grid) then
           rx=coord_nod2d(1,n2)
           ry=coord_nod2d(2,n2)
           call r2g(x, y, rx, ry)
           x=x/rad   ! degree
           y=y/rad
        else
           x=coord_nod2d(1,n2)
           y=coord_nod2d(2,n2)
        end if
        do k=1,num_layers_below_nod2d(n2)+1
           n=nod3d_below_nod2d(k,n2)
           z=coord_nod3d(3,n)
           call interp_vertical_section(num_hor_reg, num_ver_reg, hor_reg, ver_reg, raw_t, &
                1, y, z, tracer(n,1))
           call interp_vertical_section(num_hor_reg, num_ver_reg, hor_reg, ver_reg, raw_s, &
                1, y, z, tracer(n,2))
           if(tracer(n,2)>=34.62) tracer(n,3)=1.0
        end do
     end if
  end do

  tracer0=tracer  !!!

  deallocate(ver_reg, hor_reg, raw_t, raw_s, source_nodes)

  if(mype==0) write(*,*) 'source water prescribed'

end subroutine init_source_RS1
!
!--------------------------------------------------------------------
!
subroutine restore_source_RS1
  ! restore T/S source for RS1 setup
  use o_mesh
  use o_param
  use o_array
  use g_config
  use g_rotate_grid
  use g_parfe
  implicit none
 
  integer                     :: row
  real(kind=8)                :: x, y, rx, ry

  ! restoring for tracers

  do row=1,myDim_nod3d+eDim_nod3d
     
     if(rotated_grid) then
        rx=coord_nod3d(1,row)
        ry=coord_nod3d(2,row)
        call r2g(x, y, rx, ry)
        x=x/rad   ! degree
        y=y/rad
     else
        x=coord_nod3d(1,row)/rad
        y=coord_nod3d(2,row)/rad
     end if

     ! if((x>160. .and. x<170. .and. y<=-74.7) .or. (y<0.6897*x-193.3103 .and. y<-72.9)) then

     if(x>160. .and. x<171.2 .and. y<=-74.5) then   
        tracer(row,:)=tracer0(row,:)
     end if

  end do
  
end subroutine restore_source_RS1
