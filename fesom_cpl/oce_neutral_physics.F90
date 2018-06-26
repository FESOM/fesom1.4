subroutine prepare_neutral_physis
  use o_param
  use g_config
  implicit none

  call compute_neutral_slope

  call find_boundary_layer

  if(Kh_flow_depend .and. mod(istep, step_per_day)==1) then
  	call compute_neutral_diffusivity
  end if

end subroutine prepare_neutral_physis
!
!----------------------------------------------------------------------------
!
subroutine compute_neutral_diffusivity
  ! Calculate the flow dependent neutral diffusivity.
  ! Take the squared buoyancy frequency as the scaling factor. 
  ! A relative value is computed here; the diffusivity is the product
  ! of this relative value and a reference diffusivity (the diffusivity
  ! at the base of the surface boundary layer).
  use o_mesh
  use o_elements
  use o_param
  use o_array
  use g_config
  use g_parfe
  implicit none

  integer         :: n2, n3, n3up, k
  real(kind=8)    :: aux, bfsq_ref, zup, zlo
  real(kind=8)    :: bfsq_col(max_num_layers)

  do n2=1,ToDim_nod2D
     
     ! inside BL Kh is equal to the reference value
     Kh_relative(nod3d_below_nod2d(1:BL_index(n2),n2))=1.0
     if(BL_index(n2)==num_layers_below_nod2d(n2)+1) continue

     bfsq_col(1)=bfsq_3D(n2)
     n3up=n2
     do k=2, num_layers_below_nod2D(n2)+1
        n3=nod3d_below_nod2d(k,n2)
        if(k<num_layers_below_nod2d(n2)+1) then
           bfsq_col(k)=(bfsq_3D(n3up)+bfsq_3D(n3))*0.5
        else
           bfsq_col(k)=bfsq_3D(n3up)
        end if
        n3up=n3
     end do

     !interpolate to get the bfsq at the BL bottom
     zup=abs(coord_nod3d(3,nod3d_below_nod2d(BL_index(n2)-1,n2)))
     zlo=abs(coord_nod3d(3,nod3d_below_nod2d(BL_index(n2),n2)))
     bfsq_ref=bfsq_col(BL_index(n2))- &
          (bfsq_col(BL_index(n2))-bfsq_col(BL_index(n2)-1)) &
          *(zlo-BL_depth(n2))/(zlo-zup)
     
     ! if bfsq negative at the base, fine the first positive value below in the column
     if(bfsq_ref<=0.0) then
        do k=BL_index(n2), num_layers_below_nod2d(n2)+1
           if(bfsq_col(k)>0.0) then
              bfsq_ref=bfsq_col(k)
              exit
           end if
        end do
     end if

     if(bfsq_ref<=0.0) then
        do k=BL_index(n2)+1, num_layers_below_nod2d(n2)+1
           n3=nod3d_below_nod2d(k,n2)
           Kh_relative(n3)=1.0
        end do
     else
        do k=BL_index(n2)+1, num_layers_below_nod2D(n2)+1
           n3=nod3d_below_nod2d(k,n2)
           Kh_relative(n3)=max(0.2, min(bfsq_col(k)/bfsq_ref,1.0))
           ! limit from above and below
        end do
     end if

  end do

  !!call com_3d(Kh_relative)

end subroutine compute_neutral_diffusivity
!
!----------------------------------------------------------------------------
!
subroutine find_boundary_layer
  !Calculate the MLD as the BL depth
  !Define MLD following Larege etal 1997
  !MLD is the shallowest depth where the local buoyancy gradient
  !matches the maximum buoyancy gradient between the surface
  !and any discrete depth within the water column.
  use o_mesh
  use o_elements
  use o_param
  use o_array
  use g_config
  use g_parfe
  implicit none

  integer :: n2, n3, nlo, k, flag
  integer :: lay_m, lay_el(3), elem, elnodes_lo(3), elnodes2(3)
  real    :: z, dbgrad_max, aux, zm, bld_el
  real, allocatable :: smooth_array(:)

  ! find mixed layer depth and the neutral slope at the mld base
  do n2=1,myDim_nod2D
     dbgrad_max=0.0
     do k=2,num_layers_below_nod2d(n2)+1
        n3=nod3d_below_nod2d(k,n2)
        z=abs(coord_nod3d(3,n3))
        dbgrad_max=max(dbgrad_max, dbsfc_3D(n3)/z)  !1/s2
        !note: dbsfc_3d saves the value on a segment above a grid node
     end do

     flag=0
     do k=1,num_layers_below_nod2D(n2)
        n3=nod3d_below_nod2d(k,n2)
        nlo=nod3d_below_nod2d(k+1,n2)
        if(bfsq_3D(n3)>=dbgrad_max) then
!!$           Sx_neutral_base(n2)= &
!!$                sum(neutral_slope_elem(1,nod_in_elem3d(nlo)%addresses))/nod_in_elem3d(n3)%nmb
!!$           Sy_neutral_base(n2)= &
!!$                sum(neutral_slope_elem(2,nod_in_elem3d(nlo)%addresses))/nod_in_elem3d(n3)%nmb
           BL_depth(n2)=abs(coord_nod3d(3,nlo))
           flag=1
           exit
        end if
     end do
     if(flag==0) then
        BL_depth(n2)=abs(coord_nod3d(3,nod3d_below_nod2d(2,n2)))
     end if

!!$     ! limit S at the base to S_neutral_max
!!$     aux=sqrt(Sx_neutral_base(n2)*Sx_neutral_base(n2) + &
!!$          Sy_neutral_base(n2)*Sy_neutral_base(n2))
!!$     if(aux>S_neutral_max) then
!!$        aux=S_neutral_max/aux
!!$        Sx_neutral_base(n2)=Sx_neutral_base(n2)*aux
!!$        Sy_neutral_base(n2)=Sy_neutral_base(n2)*aux
!!$     end if

  end do

  call com_2D(BL_depth)

  !smooth mld
  allocate(smooth_array(ToDim_nod2d))
  do n2=1,myDim_nod2d
     smooth_array(n2)=sum(BL_depth(nghbr_nod2d(n2)%addresses))/nghbr_nod2d(n2)%nmb
  end do
  BL_depth=smooth_array
  call com_2D(BL_depth)
  deallocate(smooth_array)

  call com_2D(BL_depth)

  !index of the base of BL
  do n2=1,ToDim_nod2d
     do k=1,num_layers_below_nod2D(n2)+1
        n3=nod3d_below_nod2d(k,n2)
        if(abs(coord_nod3d(3,n3))>=BL_depth(n2) .or. &
             k==num_layers_below_nod2d(n2)+1) then
           BL_index(n2)=k
           exit
        end if
     end do
  end do

  do elem=1, myDim_elem2D                  
     if(grid_type_elem2d(elem)==1) cycle  !only apply GM to z-level grids
     elnodes2=elem2d_nodes(:,elem)
     bld_el=sum(BL_depth(elnodes2))/3.0

     lay_el=num_layers_below_nod2d(elnodes2)
     lay_m=minval(lay_el)

     do k=1,lay_m
        elnodes_lo=nod3d_below_nod2d(k+1,elnodes2)
        zm=abs(minval(coord_nod3d(3,elnodes_lo)))
        if(zm>=bld_el .or. k==lay_m) then
           Sx_neutral_base(elem)= neutral_slope(1,k,elem)
           Sy_neutral_base(elem)= neutral_slope(2,k,elem)
           exit
        end if
     end do

     aux=sqrt(Sx_neutral_base(elem)*Sx_neutral_base(elem) + &
          Sy_neutral_base(elem)*Sy_neutral_base(elem))
     if(aux>S_neutral_max) then
        aux=S_neutral_max/aux
        Sx_neutral_base(elem)=Sx_neutral_base(elem)*aux
        Sy_neutral_base(elem)=Sy_neutral_base(elem)*aux
     end if
  end do

end subroutine find_boundary_layer
!
!----------------------------------------------------------------------------
!
!!$subroutine find_boundary_layer
!!$  use o_mesh
!!$  use o_elements
!!$  use o_param
!!$  use o_array
!!$  use g_config
!!$  use g_parfe
!!$  implicit none
!!$
!!$  integer :: n2, n3, k, flag
!!$  real    :: sx, sy, aux
!!$
!!$  do n2=1,myDim_nod2D
!!$     flag=0
!!$     do k=1,num_layers_below_nod2D(n2)+1
!!$        n3=nod3d_below_nod2d(k,n2)
!!$
!!$        sx=sum(neutral_slope_elem(1,nod_in_elem3d(n3)%addresses))/nod_in_elem3d(n3)%nmb
!!$        sy=sum(neutral_slope_elem(2,nod_in_elem3d(n3)%addresses))/nod_in_elem3d(n3)%nmb
!!$
!!$        if(abs(sx)<S_neutral_max .and. abs(sy)<S_neutral_max) then
!!$           BL_index(n2)=k
!!$           BL_depth(n2)=abs(coord_nod3d(3,n3))
!!$           Sx_neutral_base(n2)=sx
!!$           Sy_neutral_base(n2)=sy
!!$           flag=1
!!$           exit
!!$        end if
!!$     end do
!!$     if(flag==0) then
!!$        BL_index(n2)=num_layers_below_nod2D(n2)+1
!!$        BL_depth(n2)=abs(coord_nod3d(3,n3))
!!$        aux=sqrt(sx*sx+sy*sy)+1.0e-10
!!$        Sx_neutral_base(n2)=sx/aux*S_neutral_max
!!$        Sy_neutral_base(n2)=sy/aux*S_neutral_max
!!$     end if
!!$  end do
!!$
!!$  call com_2D(BL_depth)
!!$  call com_2D(Sx_neutral_base)
!!$  call com_2D(Sy_neutral_base)
!!$
!!$  ssh=BL_depth
!!$  call output(1)
!!$  call par_ex
!!$  stop
!!$
!!$end subroutine find_boundary_layer
!
!----------------------------------------------------------------------------
!
subroutine compute_neutral_slope
  ! calculate neutral slopes 
  !
  ! OUTPUT:
  ! version 1: neutral_slope(3,layer,elem2d), over prism
  ! version 2: neutral_slope_elem(3, elem3d), elementwise neutral slope   
  ! 
  ! REFERENCE:
  !    McDougall, T.J. and  D.R. Jackett, 1988  
  !    On the helical nature of neutral surfaces
  !    Progress in Oceanography, vol 20, Pergamon, 153-183
  !    (about the definition of normal direction of neutral surface)
  !
  !    Griffies, S.M. et al. 1998
  !    Isoneutral Diffusion in a z-coordinate ocean model
  !    JPO, vol 28, 805-830
  !
  ! Coded by Qiang Wang, 25,11,2004
  ! Modified by Qiang Wang, 02,07,2010, add option of prismatic slope
  ! Reviewed by ??
  !-------------------------------------------------------------------

  use o_mesh
  use o_elements
  use o_param
  use o_array
  use g_config
  use g_parfe
  implicit none

  integer :: i, k, nod, elem, lay_m, elnodes2(3), elnodes(6), lay_el(3)
  integer :: elnodes_tetra(4)
  real    :: rho(6), zm, dz, dx2d(3), dy2d(3), denom, ano_1, ano_2
  real    :: rhograd_x(max_num_layers), rhograd_y(max_num_layers)
  real    :: rhograd_z(max_num_layers)
  real    :: t, s, p, sw_alpha, sw_beta
  real    ::  g_t_x, g_t_y, g_t_z, g_s_x, g_s_y, g_s_z
  real    :: inv2, inv3, inv4, inv6

  inv2=0.5
  inv3=1.0/3.0
  inv4=0.25
  inv6=1.0/6.0

  if(nslope_version==1) then
     do elem=1, myDim_elem2D                  
        if(grid_type_elem2d(elem)==1) cycle  !only apply GM to z-level grids
        elnodes2=elem2d_nodes(:,elem)
        dx2d=bafux_2d(:,elem)
        dy2d=bafuy_2d(:,elem)
        lay_el=num_layers_below_nod2d(elnodes2)
        lay_m=minval(lay_el)

        do k=1,lay_m
           elnodes(1:3)=nod3d_below_nod2d(k,elnodes2)
           elnodes(4:6)=nod3d_below_nod2d(k+1,elnodes2)
           zm=sum(coord_nod3d(3,elnodes))*inv6
           dz=1.0/(coord_nod3d(3,elnodes(1))-coord_nod3d(3,elnodes(4)))
           do i=1,6
              nod=elnodes(i)
              call fcn_density(tracer(nod,1),tracer(nod,2),zm,rho(i))
           end do
           rhograd_x(k)=(sum(dx2d*rho(1:3))+sum(dx2d*rho(4:6)))*inv2
           rhograd_y(k)=(sum(dy2d*rho(1:3))+sum(dy2d*rho(4:6)))*inv2
           rhograd_z(k)=(sum(rho(1:3))-sum(rho(4:6)))*inv3*dz
        end do

        do k=1,lay_m
           if(lay_m>1) then
              if(k==1) then
                 denom=(2.0*rhograd_z(1)+rhograd_z(2))*inv3
              elseif(k<lay_m) then
                 denom=(rhograd_z(k-1)+2.0*rhograd_z(k)+rhograd_z(k+1))*inv4
              else
                 denom=(2.0*rhograd_z(lay_m)+rhograd_z(lay_m-1))*inv3
              end if
           else
              denom=rhograd_z(1)
           end if
           if(denom<0.0) then
              denom=denom-1.0e-20
              neutral_slope(1,k,elem)=-rhograd_x(k)/denom
              neutral_slope(2,k,elem)=-rhograd_y(k)/denom
              neutral_slope(3,k,elem) = &
                   sqrt(neutral_slope(1,k,elem)**2 + neutral_slope(2,k,elem)**2)
           else
              neutral_slope(:,k,elem)=0.0
           end if

        end do
     end do

  else !version 2

     uv_rhs(1:ToDim_nod3d)=0.   ! to save memory we use uv_rhs as a temporary array

     do elem = 1, myDim_elem3d
        elnodes_tetra=elem3D_nodes(:,elem)

        ! alpha and beta
        t=inv4*sum(tracer(elnodes_tetra,1))
        s=inv4*sum(tracer(elnodes_tetra,2))
        p=abs(inv4*sum(coord_nod3D(3,elnodes_tetra))) 
        call sw_alpha_beta(t,s,p,sw_alpha,sw_beta)

        g_t_z = sum(bafuz_3d(:,elem)*tracer(elnodes_tetra,1))
        g_s_z = sum(bafuz_3d(:,elem)*tracer(elnodes_tetra,2))
        denom = -sw_alpha*g_t_z + sw_beta*g_s_z

        if(denom>=0.0) then
           uv_rhs(elnodes_tetra)=1.
           neutral_slope_elem(:,elem)=0.0
           cycle
        end if

        denom=denom-1.e-20

        g_t_x = sum(bafux_3d(:,elem)*tracer(elnodes_tetra,1))
        g_s_x = sum(bafux_3d(:,elem)*tracer(elnodes_tetra,2))
        g_t_y = sum(bafuy_3d(:,elem)*tracer(elnodes_tetra,1))
        g_s_y = sum(bafuy_3d(:,elem)*tracer(elnodes_tetra,2))  

        ano_1 = -sw_alpha*g_t_x + sw_beta*g_s_x
        ano_2 = -sw_alpha*g_t_y + sw_beta*g_s_y

        neutral_slope_elem(1,elem) = -ano_1/denom
        neutral_slope_elem(2,elem) = -ano_2/denom
        neutral_slope_elem(3,elem) = &
             sqrt(neutral_slope_elem(1,elem)**2 + neutral_slope_elem(2,elem)**2)
     end do

     ! in case of static instability
     do i=1,myDim_nod3D
        if(uv_rhs(i)>0.5) then
           neutral_slope_elem(:,nod_in_elem3d(i)%addresses)=0.0
        end if
     end do

  end if
end subroutine compute_neutral_slope
!
!----------------------------------------------------------------------------
!
subroutine sw_alpha_beta(t1, s1, p1, alpha, beta)
  !   A function to calculate the thermal expansion coefficient
  !   and saline contraction coefficient.
  !
  ! INPUT: t1 (c), s1 (psu), p1 (db)
  !
  ! OUTPUT:
  !    alpha = Thermal expansion coeff (alpha) [degree_C.^-1]
  !    beta  = Saline contraction coeff (beta) [psu.^-1]
  !
  ! REFERENCE:
  !    McDougall, T.J. 1987.  Neutral Surfaces
  !    Journal of Physical Oceanography, vol 17, 1950-1964,
  !
  ! CHECK VALUE:
  !    sw_beta=0.72088e-3 psu^-1 @ S=40.0psu, ptmp=10.0C (ITS-90), p=4000db
  !    a_over_b=0.34765 psu*C^-1 @ S=40.0psu, ptmp=10.0C, p=4000db
  !
  ! Coded by Qiang Wang, 25,11,2004
  ! Reviewed by ??
  !-----------------------------------------------------------------

  implicit none

  real :: t1, s1, p1, alpha, beta
  real :: t1_2,t1_3,t1_4,p1_2,p1_3,s35,s35_2 
  real :: a_over_b    

  t1_2 = t1*t1
  t1_3 = t1_2*t1
  t1_4 = t1_3*t1
  p1_2 = p1*p1
  p1_3 = p1_2*p1
  s35 = s1-35.0_8
  s35_2 = s35*s35

  ! calculate beta
  beta = 0.785567e-3 - 0.301985e-5*t1 &
       + 0.555579e-7*t1_2 - 0.415613e-9*t1_3 &
       + s35*(-0.356603e-6 + 0.788212e-8*t1 &
       + 0.408195e-10*p1 - 0.602281e-15*p1_2) &
       + s35_2*(0.515032e-8) & 
       + p1*(-0.121555e-7 + 0.192867e-9*t1 - 0.213127e-11*t1_2) &
       + p1_2*(0.176621e-12 - 0.175379e-14*t1) &
       + p1_3*(0.121551e-17)

  ! calaculate the thermal expansion / saline contraction ratio
  a_over_b = 0.665157e-1 + 0.170907e-1*t1 &
       - 0.203814e-3*t1_2 + 0.298357e-5*t1_3 &
       - 0.255019e-7*t1_4 &
       + s35*(0.378110e-2 - 0.846960e-4*t1 &
       - 0.164759e-6*p1 - 0.251520e-11*p1_2) &
       + s35_2*(-0.678662e-5) &
       + p1*(0.380374e-4 - 0.933746e-6*t1 + 0.791325e-8*t1_2) &
       + p1_2*t1_2*(0.512857e-12) &
       - p1_3*(0.302285e-13)

  ! calculate alpha
  alpha = a_over_b*beta

end subroutine sw_alpha_beta
!
