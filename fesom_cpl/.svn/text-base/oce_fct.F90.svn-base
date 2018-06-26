! This file collect subroutines implementing FE-FCT
! advection scheme by Loehner et al.
! There is a tunable paremeter gamma_fct in ts_solve_low_order and fem_fct.
! Increasing it leads to positivity preserving solution.
!----------------------------------------------------------------------------

subroutine fct_init
  ! Init. the fct scheme
  !
  ! Coded by Sergey Danilov
  ! Reviewed by Qiang Wang
  !--------------------------------------------------------
  
  use o_param
  use o_mesh
  use o_array
  use g_config
  use g_parfe
  !
  allocate(tral(myDim_nod3D+eDim_nod3D,num_tracer))       
  allocate(trafluxes(4,myDim_elem3D))
  allocate(pplus(myDim_nod3D+eDim_nod3D), pminus(myDim_nod3D+eDim_nod3D)) 
  tral=0.0

end subroutine fct_init
!
!----------------------------------------------------------------------------
!
subroutine fct_tracer_solve
  ! The driving routine for the fct scheme
  !
  ! Coded by Sergey Danilov
  ! Reviewed by Qiang Wang
  ! Adapted to the tetrahedra code by Qiang Wang
  !--------------------------------------------------------
  
  use o_solver
  use o_param
  use o_array
  use g_config
  use g_parfe

  integer     :: i

  ! high order solution
  if(lump_ts_matrix) then
     call tracer_solve
  else
     do i=1,num_tracer
        call solve(solve_tra+i-1)
     end do
  endif 

  ! low order solution
  call tra_solve_low_order

  ! fct
  do i=1,num_tracer
     call fem_fct(i)
  end do
  do i=1,num_tracer
     call com_3d(tracer(:,i))
  end do

end subroutine fct_tracer_solve
!
!----------------------------------------------------------------------------
!
subroutine tra_solve_low_order
  ! Low-order solution
  ! One adds diffusive contribution to the rhs. It is realized as
  ! difference between the consistent and lumped mass matrices
  ! acting on the field from the previous time step.   
  !
  ! Coded by Sergey Danilov
  ! Reviewed by Qiang Wang
  !--------------------------------------------------------

  use o_matrices
  use o_mesh
  use o_array
  use o_param
  use g_config
  use g_parfe
  implicit none
  
  integer      ::  j, row, clo, clo2, cn, location(100)
  
  do row=1,myDim_nod3D               
     clo=tsstiff%rowptr(row)-tsstiff%rowptr(1)+1  
     clo2=tsstiff%rowptr(row+1)-tsstiff%rowptr(1) 
     cn=clo2-clo+1
     location(1:cn)=nghbr_nod3D(row)%addresses    
     do j=1,num_tracer
        tral(row,j)=(tracer_rhs(row,j)+gamma_fct*sum(tsstiff%values(clo:clo2)* &
             tracer(location(1:cn),j)))/ts_lump(row) + &
             (1.-gamma_fct)*tracer(row,j)
     end do
   end do
  
   do j=1,num_tracer
      call com_3D(tral(:,j))     ! solution must be known to neighbours
   end do
end subroutine tra_solve_low_order
!
!----------------------------------------------------------------------------
!
subroutine fem_fct(tra_id)
  ! Flux corrected transport algorithm for tracer advection
  !
  ! It is based on Loehner et al. (Finite-element flux-corrected 
  ! transport (FEM-FCT) for the Euler and Navier-Stokes equation, 
  ! Int. J. Numer. Meth. Fluids, 7 (1987), 1093--1109) as described by Kuzmin and
  ! Turek. (kuzmin@math.uni-dortmund.de) 
  !
  ! Steps:
  !
  ! Construct a low-order solution by adding an artificial diffusion 
  ! to the rhs and lumping the mass matrix
  !
  ! Compute a high-order solution (run ts_solve)
  ! 
  ! Limit antidiffusive fluxes 
  ! 
  ! Update the low-order solution to the high order but using 
  ! the limited fluxes
  
  ! Coded by Sergey Danilov
  ! Reviewed by Qiang Wang
  ! Adapted to the tetrahedra code by Qiang Wang
  !--------------------------------------------------------
  
  use o_matrices
  use o_mesh
  use o_elements
  use o_array
  use o_param
  use g_config
  use g_parfe
  implicit none

  integer        :: tra_id
  integer        :: i, n, q, row
  integer        :: elem, elnodes(4)
  real(kind=8)   :: flux, ae, vol, aux(4), icoef(4,4), inv20
  real(kind=8), allocatable :: tmax(:), tmin(:) 

  allocate(tmax(myDim_nod3D), tmin(myDim_nod3D))
  inv20=1.0_8/20.0_8

  !==========================
  ! Compute elemental antidiffusive fluxes to nodes
  !==========================
  ! This is the most unpleasant part: it takes memory and time. 
  ! For every element we need its antidiffusive contribution to 
  ! each of its 4 nodes
  !
  ! Auxiliary elemental operator (lumped mass matrix - mass matrix)
  icoef=-inv20
  do n=1,4   
     icoef(n,n)=3.0*inv20
  end do
  ! antidiffusive fluxes 
  do elem=1, myDim_elem3D            
     elnodes=elem3D_nodes(:,elem)
     vol=voltetra(elem)*dt_inv   
     aux=gamma_fct*tracer(elnodes,tra_id) + dtracer(elnodes,tra_id)
     do q=1,4       
        trafluxes(q,elem)=sum(icoef(:,q)*aux)*vol/ts_lump(elnodes(q))
     end do
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
  do row=1, myDim_nod3D               
     n=nghbr_nod3D(row)%nmb
     !-----------------------
     !if(geolat(row)<45.0*rad .or. (geolat(row)<76.5*rad .and. geolon(row)>-80.0*rad .and. geolon(row)<20.0*rad)) then
     !  tmax(row)=maxval(tral(nghbr_nod3D(row)%addresses(1:n),tra_id))
     !  tmin(row)=minval(tral(nghbr_nod3D(row)%addresses(1:n),tra_id))
     !else
       tmax(row)=max(maxval(tral(nghbr_nod3D(row)%addresses(1:n),tra_id)),  maxval(tracer(nghbr_nod3D(row)%addresses(1:n),tra_id)))
       tmin(row)=min(minval(tral(nghbr_nod3D(row)%addresses(1:n),tra_id)),  minval(tracer(nghbr_nod3D(row)%addresses(1:n),tra_id)))
     !end if
     !-----------------------
     ! Admissible increments
     tmax(row)=tmax(row)-tral(row,tra_id)    
     tmin(row)=tmin(row)-tral(row,tra_id)   
  end do

  !=========================
  ! Sums of positive/negative fluxes to node row
  !=========================
  pplus=0.
  pminus=0.
  do elem=1, myDim_elem3D           
     elnodes=elem3D_nodes(:,elem)
     do q=1,4
        n=elnodes(q) 
        flux=trafluxes(q,elem)    
        if (flux>0.) then
           pplus(n)=pplus(n)+flux
        else
           pminus(n)=pminus(n)+flux	  
        end if
     end do
  end do

  !========================
  ! The least upper bound for the correction factors
  !========================
  do n=1,myDim_nod3D                
     flux=pplus(n)
     if (abs(flux)>0.) then    !!!!
        pplus(n)=min(1.0,tmax(n)/flux)  
     else
        pplus(n)=0.
     end if
     flux=pminus(n)
     if (abs(flux)>0.) then    !!!! 
        pminus(n)=min(1.0,tmin(n)/flux)  
     else
        pminus(n)=0.
     end if
  end do
  ! pminus and pplus are to be known to neighbouring PE
  call com_3D(pplus)
  call com_3D(pminus)	   

  !========================	 
  ! Limiting
  !========================	 
  do elem=1, myDim_elem3D           
     elnodes=elem3D_nodes(:,elem)
     ae=1.0
     do q=1,4
        n=elnodes(q)  
        flux=trafluxes(q,elem)            
        if(flux>=0.) ae=min(ae,pplus(n))
        if(flux<0.) ae=min(ae,pminus(n))
     end do
     trafluxes(:,elem)=ae*trafluxes(:,elem) 
  end do

  !==========================
  ! Update the solution
  !==========================
  do n=1,myDim_nod3D                 
     tracer(n,tra_id)=tral(n,tra_id)
  end do
  do elem=1, myDim_elem3D             
     elnodes=elem3D_nodes(:,elem)
     do q=1,4
        n=elnodes(q)  
        tracer(n,tra_id)=tracer(n,tra_id)+trafluxes(q,elem) 
     end do
  end do

  deallocate(tmin, tmax)
end subroutine fem_fct
!
!----------------------------------------------------------------------------
!
