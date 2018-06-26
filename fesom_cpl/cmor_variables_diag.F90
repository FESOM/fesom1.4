module cmor_variables_diag
  implicit none
  public compute_diag, volo, opottemptend, pbo, soga, thetaoga, tos, sos, siarean, siareas, siextentn, siextents, sivoln, sivols
  private
  
  logical, save :: initialized = .false.
  real(kind=8), save :: volo
  real(kind=8), save, allocatable :: opottemptend(:)
  real(kind=8), save, allocatable :: pbo(:)
  real(kind=8), save, allocatable :: tos(:)
  real(kind=8), save, allocatable :: sos(:)
  real(kind=8), save :: soga, thetaoga
  real(kind=8), save :: siarean, siareas, siextentn, siextents, sivoln, sivols

contains

  subroutine init_diag
    use g_parfe, only : myDim_nod3D, myDim_nod2D
    use o_mesh, only : nod_in_elem3D
    use o_elements, only : voltetra, cluster_vol_3d
    use g_parfe, only : MPI_COMM_FESOM, MPI_DOUBLE_PRECISION, MPI_IN_PLACE, MPI_SUM, mype
    integer k, n3, el, ierr

    ! variable volo
    volo  = 0.
    do n3=1, myDim_nod3D
       volo = volo+cluster_vol_3d(n3)
    end do
    call MPI_AllREDUCE(MPI_IN_PLACE, volo, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FESOM, ierr)
  
    allocate(opottemptend(myDim_nod3D))
    allocate(pbo(myDim_nod2D))
    allocate(tos(myDim_nod2D))
    allocate(sos(myDim_nod2D))

    initialized = .true.
  end


  subroutine compute_diag
    use g_config, only : dt
    use o_array, only : dtracer, hpressure, ssh, tracer, density_insitu
    use g_parfe, only : myDim_nod3D, myDim_nod2D, ToDim_nod3D
    use o_param, only : vcpw, g, rho0, num_tracer
    use o_mesh, only : num_layers_below_nod2d, nod3d_below_nod2d, nod_in_elem3D, nod_in_elem2D, coord_nod2D, coord_nod3D
    use o_elements, only : voltetra, voltriangle, cluster_vol_3d, cluster_area_2d
    use g_parfe, only : MPI_COMM_FESOM, MPI_DOUBLE_PRECISION, MPI_IN_PLACE, MPI_SUM
    use i_array, only : a_ice, m_ice
    integer n2, k, n3, el, ierr

    real(kind=8)   :: z_up, z_lo, dens
    real(kind=8), allocatable, target, dimension(:,:), save :: previous_tracer
    integer        :: node_up, node_lo

    if(.not. allocated(previous_tracer)) then
      allocate(previous_tracer(ToDim_nod3D,num_tracer))
      previous_tracer = tracer
    end if

    if(.not. initialized) call init_diag
    
opottemptend=0.
pbo         =0.
do n2=1,myDim_nod2d
     pbo(n2)=ssh(n2)*rho0*g ! bottom pressure computation / SSH contribution

     node_up = nod3D_below_nod2D(1, n2)
     z_up    = coord_nod3d(3,node_up)

     do k=2, num_layers_below_nod2D(n2)+1
          node_lo = nod3D_below_nod2D(k, n2)
          z_lo    = coord_nod3d(3,node_lo)
          dens    = 0.5_8*(density_insitu(node_up)+density_insitu(node_lo))
          pbo(n2) = pbo(n2)+g*(z_up-z_lo)*dens !bottom pressure computation / steric contribution !

!The tendency diagnostics are expressed as rates  of  change  of heat  and  salt  content  in  grid  cells,  i.e. d(m*C_p*T)/dt, where
!T is temperature and m is  the mass per unit area of the grid cell and C_p the specific heat.
!In other words, m=rho*hight_of_the_gridcell and the prismatic approximation is assumed below:

          opottemptend(node_up) =opottemptend(node_up)+(tracer(node_up,1)-previous_tracer(node_up,1))/dt*vcpw*(z_up-z_lo)*0.5  ! dtracer is not "dtracer"! Qiang
          opottemptend(node_lo) =opottemptend(node_lo)+(tracer(node_lo,1)-previous_tracer(node_lo,1))/dt*vcpw*(z_up-z_lo)*0.5

          node_up = node_lo
          z_up    = z_lo
     end do
end do 
    
    ! variable soga, thetaoga
    soga = 0.
    thetaoga = 0.
    do n3=1, myDim_nod3D
       soga = soga+tracer(n3, 2)*cluster_vol_3d(n3)
       thetaoga = thetaoga+tracer(n3, 1)*cluster_vol_3d(n3)
    end do    
    call MPI_AllREDUCE(MPI_IN_PLACE, soga, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FESOM, ierr)
    call MPI_AllREDUCE(MPI_IN_PLACE, thetaoga, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FESOM, ierr)
    soga = soga/volo
    thetaoga = thetaoga/volo
  
    ! variable tos, sos
    do n2=1, myDim_nod2D
      n3 = nod3D_below_nod2D(1, n2)
      tos(n2) = tracer(n3,1)
      sos(n2) = tracer(n3,2)
    end do  

    ! variables siarean, siareas, siextentn, siextents, sivoln, sivols
    siarean = 0.
    siareas = 0.
    siextentn = 0.
    siextents = 0.
    sivoln = 0.
    sivols = 0.
    do n2=1, myDim_nod2D
       
       if (coord_nod2D(2, n2)>0.) then
          siarean=siarean+a_ice(n2)*cluster_area_2d(n2)
          if (a_ice(n2)>=.15) siextentn=siextentn+cluster_area_2d(n2)
          sivoln = sivoln + m_ice(n2)*cluster_area_2d(n2)
       else
          siareas=siareas+a_ice(n2)*cluster_area_2d(n2)
          if (a_ice(n2)>=.15) siextents=siextents+cluster_area_2d(n2)
          sivols = sivols + m_ice(n2)*cluster_area_2d(n2)
       end if
    end do  
    call MPI_AllREDUCE(MPI_IN_PLACE, siarean, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FESOM, ierr)
    call MPI_AllREDUCE(MPI_IN_PLACE, siareas, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FESOM, ierr)
    call MPI_AllREDUCE(MPI_IN_PLACE, siextentn, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FESOM, ierr)
    call MPI_AllREDUCE(MPI_IN_PLACE, siextents, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FESOM, ierr)
    siarean = siarean/1.e12
    siareas = siareas/1.e12
    siextentn = siextentn/1.e12
    siextents = siextents/1.e12

    call MPI_AllREDUCE(MPI_IN_PLACE, sivoln, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FESOM, ierr)
    call MPI_AllREDUCE(MPI_IN_PLACE, sivols, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_FESOM, ierr)
    sivoln = sivoln/1.e9
    sivols = sivols/1.e9
    
    previous_tracer = tracer
  end
end

