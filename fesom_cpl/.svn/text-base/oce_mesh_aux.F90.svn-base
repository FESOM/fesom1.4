subroutine mesh_cluster_setup
  use g_config
  implicit none

  call find_cluster_area
  call find_cluster_vol

  if (.not.r_restart) then
     call write_initial_mesh_diag
  end if

end subroutine mesh_cluster_setup
!
!---------------------------------------------------------------------------
!
subroutine find_cluster_area
  !
  ! Coded by Qiang Wang
  ! Reviewed by ??
  !-------------------------------------------------------------
  
  use o_MESH
  use o_ELEMENTS
  use g_PARFE
  use g_rotate_grid
  implicit none
  
  integer         :: i, elem, elnodes(3)
  real(kind=8)    :: inv3, vol
  real(kind=8)    :: x, y  
  integer         :: nn, ns  


  inv3=1.0/3.0_8

  allocate(cluster_area_2d(ToDim_nod2D))
  cluster_area_2d=0.0

  do elem=1,myDim_elem2d
     elnodes=elem2d_nodes(:,elem)
     vol=voltriangle(elem)*inv3
     cluster_area_2d(elnodes)=cluster_area_2d(elnodes)+vol
  end do

  !call com_2d(cluster_area_2d)

  vol=0.0
  do i=1,myDim_nod2d
     vol=vol+cluster_area_2d(i)
  end do

  call MPI_Barrier(MPI_COMM_FESOM,MPIerr)
  ocean_area=0.0
  call MPI_AllREDUCE(vol, ocean_area, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
       MPI_COMM_FESOM, MPIerr)
       

#if defined (__oasis) || defined (__uncplecham6)  
  nn=0
  ns=0  
  allocate(lump2d_north(myDim_nod2D), lump2d_south(myDim_nod2D))
  lump2d_north=0.
  lump2d_south=0.
  do i=1, myDim_nod2D
     call r2g(x, y, coord_nod2d(1, i), coord_nod2d(2, i))      
     if (y>0) then 
        nn=nn+1
        lump2d_north(i)=sum(voltriangle(nod_in_elem2d(i)%addresses))/3.
        lump2d_north(i)=abs(lump2d_north(i))
     else
        ns=ns+1     
        lump2d_south(i)=sum(voltriangle(nod_in_elem2d(i)%addresses))/3.
        lump2d_south(i)=abs(lump2d_south(i))
     end if	   
  end do   

  if (nn>0) allocate(ind_north(nn))
  if (ns>0) allocate(ind_south(ns))
  ns=0
  nn=0
  do i=1, myDim_nod2D
     call r2g(x, y, coord_nod2d(1, i), coord_nod2d(2, i))
     if (y>0) then
        nn=nn+1
	ind_north(nn)=i
     else
        ns=ns+1
	ind_south(ns)=i	
     end if	     
  end do     
#endif           
end subroutine find_cluster_area
!
!==========================================================================
!
subroutine find_cluster_vol
  !
  ! Coded by Qiang Wang
  ! Reviewed by ??
  !-------------------------------------------------------------
  
  use o_MESH
  use o_ELEMENTS
  use g_PARFE
  implicit none
  
  integer         :: i, elem, elnodes(4)
  real(kind=8)    :: inv4, vol

  inv4=1.0/4.0_8

  allocate(cluster_vol_3d(ToDim_nod3D))
  cluster_vol_3d=0.0

  do elem=1,myDim_elem3d
     elnodes=elem3d_nodes(:,elem)
     vol=voltetra(elem)*inv4
     cluster_vol_3d(elnodes)=cluster_vol_3d(elnodes)+vol
  end do

  !call com_3d(cluster_vol_3d)

  vol=0.0
  do i=1,myDim_nod3d
     vol=vol+cluster_vol_3d(i)
  end do

  call MPI_Barrier(MPI_COMM_FESOM,MPIerr)
  ocean_vol=0.0
  call MPI_AllREDUCE(vol, ocean_vol, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
       MPI_COMM_FESOM, MPIerr)
  
end subroutine find_cluster_vol
!
!==========================================================================
!
subroutine update_cluster_vol
  ! Update cluster volume in case of non-linear free surface.
  ! Based on the current version: only the surface nodes are moving,
  ! so we only need to update for the first two layers of nodes.
  !
  ! Coded by Qiang Wang, 10,02,2012
  ! Reviewed by ??
  !-------------------------------------------------------------
  
  use o_MESH
  use o_ELEMENTS
  use g_PARFE
  implicit none
  
  integer         :: i, k, nod 
  real(kind=8)    :: inv4

  inv4=1.0/4.0_8

  do i=1,myDim_nod2d
     do k=1,2
        nod=nod3d_below_nod2d(k,i)
        cluster_vol_3d(nod)=sum(voltetra(nod_in_elem3D(nod)%addresses))*inv4
     end do
  end do

end subroutine update_cluster_vol
