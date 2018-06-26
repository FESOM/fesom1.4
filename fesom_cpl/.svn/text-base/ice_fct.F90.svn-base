! The FCT ice advection scheme

subroutine fct_ice_init
  ! fct init.
  !
  ! Coded by Sergey Danilov
  ! Reviewed by Qiang Wang
  !----------------------------------------------------------

  use o_mesh
  use i_array
  use g_parfe
  implicit none
  integer        :: n2

  n2=myDim_nod2D+eDim_nod2D   

  allocate(m_icel(n2), a_icel(n2), m_snowl(n2))  ! low-order solutions
  allocate(icefluxes(myDim_elem2D,3))
  allocate(icepplus(n2), icepminus(n2))

  m_icel=0.0
  a_icel=0.0 
  m_snowl=0.0

end subroutine fct_ice_init
!
!----------------------------------------------------------------------------
!
subroutine fct_ice_solve
  ! driving routine of the ice fct scheme
  !
  ! Coded by Sergey Danilov
  ! Reviewed by Qiang Wang
  !----------------------------------------------------------
  
  use i_array
  use i_solver
  use i_dyn_parms
  use g_PARFE
  implicit none

  if(lump_ice_matrix) then
     call ice_solve    
  else                            
     call solveIce(solve_m_ice)       
     call solveIce(solve_a_ice)       
     call solveIce(solve_m_snow)
     call com_2D(dm_ice)
     call com_2D(da_ice)
     call com_2D(dm_snow)      
  end if

  call ice_solve_low_order

  call ice_fem_fct(1)    ! m_ice   
  call ice_fem_fct(2)    ! a_ice
  call ice_fem_fct(3)    ! m_snow

  call com_2D(m_ice)
  call com_2D(a_ice)
  call com_2D(m_snow)
end subroutine fct_ice_solve
!
!----------------------------------------------------------------------------
!
subroutine ice_solve_low_order
  ! Low-order solution
  ! It is assumed that m_ice, a_ice and m_snow from the previous time step 
  ! are known at 1:myDim_nod2D+eDim_nod2D.
  ! One adds diffusive contribution to the rhs. It is realized as
  ! difference between the consistent and lumped mass matrices
  ! acting on the field from the previous time step. The mass matrix on the 
  ! lhs is replaced with lumped one. 
  !
  ! Coded by Sergey Danilov
  ! Reviewed by Qiang Wang
  !----------------------------------------------------------  
    
  use o_mesh
  use i_array
  use i_dyn_parms
  use g_parfe
  implicit none

  integer      :: row, clo, clo2, cn, location(100)
  real(kind=8) :: gamma

  gamma=ice_gamma_fct     
  do row=1,myDim_nod2D             
     clo=icestiff%rowptr(row)-icestiff%rowptr(1)+1  
     clo2=icestiff%rowptr(row+1)-icestiff%rowptr(1) 
     cn=clo2-clo+1
     location(1:cn)=nghbr_nod2D(row)%addresses      
     m_icel(row)=(rhs_m(row)+gamma*sum(icestiff%values(clo:clo2)* &
          m_ice(location(1:cn))))/ice_lump(row) + &
          (1.-gamma)*m_ice(row)
     a_icel(row)=(rhs_a(row)+gamma*sum(icestiff%values(clo:clo2)* &
          a_ice(location(1:cn))))/ice_lump(row) + &
          (1.-gamma)*a_ice(row)
     m_snowl(row)=(rhs_ms(row)+gamma*sum(icestiff%values(clo:clo2)* &
          m_snow(location(1:cn))))/ice_lump(row) + &
          (1.-gamma)*m_snow(row)
  end do

  call com_2D(m_icel)
  call com_2D(a_icel)
  call com_2D(m_snowl)
end subroutine ice_solve_low_order
!
!----------------------------------------------------------------------------
!
subroutine ice_fem_fct(tr_array_id)
  ! Flux corrected transport algorithm for ice advection
  !
  ! It is based on Loehner et al. (Finite-element flux-corrected 
  ! transport (FEM-FCT) for the Euler and Navier-Stokes equation, 
  ! Int. J. Numer. Meth. Fluids, 7 (1987), 1093--1109) as described by Kuzmin and
  ! Turek. (kuzmin@math.uni-dortmund.de) 
  !
  ! Coded by Sergey Danilov
  ! Reviewed by Qiang Wang
  !----------------------------------------------------------  
  
  use o_mesh
  use o_elements
  use i_array
  use i_dyn_parms
  use g_parfe
  use g_config
  implicit none

  integer        :: tr_array_id
  integer        :: n, q, elem, elnodes(3), row
  real(kind=8), allocatable, dimension(:) :: tmax, tmin 
  real(kind=8)   :: vol, flux, ae, gamma, inv12, icoef(3,3)

  inv12=1.0_8/12.0_8

  gamma=ice_gamma_fct         

  !==========================
  ! Compute elemental antidiffusive fluxes to nodes
  !==========================
  ! This is the most unpleasant part --- 
  ! it takes memory and time. For every element 
  ! we need its antidiffusive contribution to 
  ! each of its 3 nodes

  allocate(tmax(myDim_nod2D), tmin(myDim_nod2D))

  ! Auxiliary elemental operator (mass matrix- lumped mass matrix)
  icoef=1.0
  do n=1,3   ! three upper nodes
     icoef(n,n)=-2.0
  end do
  icoef=icoef*inv12

  do elem=1, myDim_elem2D     
     elnodes=elem2D_nodes(:,elem)
     vol=voltriangle(elem)*dt_inv

     if (tr_array_id==1) then
        do q=1,3       
           icefluxes(elem,q)=-sum(icoef(:,q)*(gamma*m_ice(elnodes) + & 
                dm_ice(elnodes)))*vol/ice_lump(elnodes(q))
        end do
     end if

     if (tr_array_id==2) then
        do q=1,3       
           icefluxes(elem,q)=-sum(icoef(:,q)*(gamma*a_ice(elnodes) + &  
                da_ice(elnodes)))*vol/ice_lump(elnodes(q))
        end do
     end if

     if (tr_array_id==3) then
        do q=1,3       
           icefluxes(elem,q)=-sum(icoef(:,q)*(gamma*m_snow(elnodes) + &  
                dm_snow(elnodes)))*vol/ice_lump(elnodes(q))
        end do
     end if
  end do

  !==========================   
  ! Screening the low-order solution
  !==========================
  ! TO BE ADDED IF FOUND NECESSARY
  ! Screening means comparing low-order solutions with the
  ! solution on the previous time step and using whichever 
  ! is greater/smaller in computations of max/min below

  !==========================
  ! Cluster min/max
  !==========================
  if (tr_array_id==1) then
     do row=1, myDim_nod2D       
        n=nghbr_nod2D(row)%nmb
        tmax(row)=max(maxval(m_icel(nghbr_nod2D(row)%addresses(1:n))), maxval(m_ice(nghbr_nod2D(row)%addresses(1:n))))
        tmin(row)=min(minval(m_icel(nghbr_nod2D(row)%addresses(1:n))), minval(m_ice(nghbr_nod2D(row)%addresses(1:n))))    
        ! Admissible increments
        tmax(row)=tmax(row)-m_icel(row)                             
        tmin(row)=tmin(row)-m_icel(row)                             
     end do
  end if

  if (tr_array_id==2) then
     do row=1, myDim_nod2D          
        n=nghbr_nod2D(row)%nmb
        tmax(row)=max(maxval(a_icel(nghbr_nod2D(row)%addresses(1:n))),  maxval(a_ice(nghbr_nod2D(row)%addresses(1:n)))) 
        tmin(row)=min(minval(a_icel(nghbr_nod2D(row)%addresses(1:n))),  minval(a_ice(nghbr_nod2D(row)%addresses(1:n))))   
        ! Admissible increments
        tmax(row)=tmax(row)-a_icel(row)                            
        tmin(row)=tmin(row)-a_icel(row)                           
     end do
  end if

  if (tr_array_id==3) then
     do row=1, myDim_nod2D       
        n=nghbr_nod2D(row)%nmb
        tmax(row)=max(maxval(m_snowl(nghbr_nod2D(row)%addresses(1:n))), maxval(m_snow(nghbr_nod2D(row)%addresses(1:n))))  
        tmin(row)=min(minval(m_snowl(nghbr_nod2D(row)%addresses(1:n))), minval(m_snow(nghbr_nod2D(row)%addresses(1:n))))   
        ! Admissible increments
        tmax(row)=tmax(row)-m_snowl(row)                         
        tmin(row)=tmin(row)-m_snowl(row)                            
     end do
  end if


  !=========================
  ! Sums of positive/negative fluxes to node row
  !=========================

  icepplus=0.
  icepminus=0.
  do elem=1, myDim_elem2D          
     elnodes=elem2D_nodes(:,elem)
     do q=1,3
        n=elnodes(q) 
        flux=icefluxes(elem,q)     
        if (flux>0) then
           icepplus(n)=icepplus(n)+flux
        else
           icepminus(n)=icepminus(n)+flux	  
        end if
     end do
  end do

  !========================
  ! The least upper bound for the correction factors
  !========================
  do n=1,myDim_nod2D              
     flux=icepplus(n)
     if (abs(flux)>0) then
        icepplus(n)=min(1.0,tmax(n)/flux)     
     else
        icepplus(n)=0.
     end if

     flux=icepminus(n)
     if (abs(flux)>0) then
        icepminus(n)=min(1.0,tmin(n)/flux)      
     else
        icepminus(n)=0.
     end if
  end do
  ! pminus and pplus are to be known to neighbouting PE
  call com_2D(icepminus)
  call com_2D(icepplus) 

  !========================	 
  ! Limiting
  !========================	 
  do elem=1, myDim_elem2D                                          
     elnodes=elem2D_nodes(:,elem)
     ae=1.0
     do q=1,3
        n=elnodes(q)  
        flux=icefluxes(elem,q)     
        if(flux>=0.) ae=min(ae,icepplus(n))
        if(flux<0.) ae=min(ae,icepminus(n))
     end do
     icefluxes(elem,:)=ae*icefluxes(elem,:) 
     !if (ae.le.0.0) write (*,*) 'ae is too large', ae 
  end do


  !==========================
  ! Update the solution 
  !==========================
  if(tr_array_id==1) then
     do n=1,myDim_nod2D          
        m_ice(n)=m_icel(n)
     end do
     do elem=1, myDim_elem2D      
        elnodes=elem2D_nodes(:,elem)
        do q=1,3
           n=elnodes(q)  
           m_ice(n)=m_ice(n)+icefluxes(elem,q)   
        end do
     end do
  end if

  if(tr_array_id==2) then
     do n=1,myDim_nod2D           
        a_ice(n)=a_icel(n)
     end do
     do elem=1, myDim_elem2D      
        elnodes=elem2D_nodes(:,elem)
        do q=1,3
           n=elnodes(q)  
           a_ice(n)=a_ice(n)+icefluxes(elem,q)    
        end do
     end do
  end if

  if(tr_array_id==3) then
     do n=1,myDim_nod2D         
        m_snow(n)=m_snowl(n)
     end do
     do elem=1, myDim_elem2D     
        elnodes=elem2D_nodes(:,elem)
        do q=1,3
           n=elnodes(q)  
           m_snow(n)=m_snow(n)+icefluxes(elem,q)    
        end do
     end do
  end if

  deallocate(tmin, tmax)
end subroutine ice_fem_fct
