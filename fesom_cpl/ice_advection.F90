subroutine ice_matrices_setup
  use g_parfe
  implicit none

  call icestiff_matrix
  call icestiff_fill_tg  ! fill mass matrix and get lumped mass matrix 

end subroutine ice_matrices_setup
!
!-------------------------------------------------------------------------
!
subroutine icestiff_matrix
  ! Sets the structure of the stiffness matrix for ice advection
  use i_array
  use o_MESH
  use o_ELEMENTS
  use o_PARAM
  use g_PARFE
  use g_config
  implicit none
  !
  integer                           :: k
  integer, allocatable              :: pnza(:), rpnza(:) 
  integer                           :: nini, nend
  !
  ! a)
  icestiff%dim=nod2D
  allocate(icestiff%rowptr(myDim_nod2D+1))         
  icestiff%rowptr(1)=1
  do k=1,myDim_nod2D                               
     icestiff%rowptr(k+1)=icestiff%rowptr(k)+nghbr_nod2D(k)%nmb
  end do
  icestiff%nza=icestiff%rowptr(myDim_nod2D+1)-1    
  !

  ! b)
  allocate(icestiff%colind(icestiff%nza))
  allocate(icestiff%values(icestiff%nza))
  icestiff%values=0.0
  !
  ! =================================            
  ! Exchange between PE to know their 
  ! numbers of non-zero entries (nza):
  ! =================================     
  allocate(pnza(npes),rpnza(npes))
  pnza(1:npes)=0
  pnza(mype+1)=icestiff%nza
  call MPI_Barrier(MPI_COMM_FESOM,MPIerr)
  call MPI_AllREDUCE( pnza, rpnza, &
       npes, MPI_INTEGER,MPI_SUM, &
       MPI_COMM_FESOM, MPIerr)

  if (mype==0) then
     k=0
  else
     k=sum(rpnza(1:mype))  ! This is offset for mype 
  end if
  icestiff%rowptr=icestiff%rowptr+k       

  ! =================================
  ! replace local nza with a global one
  ! =================================
  icestiff%nza=sum(rpnza(1:npes))             
  deallocate(rpnza,pnza)

  !c)
  do k=1,myDim_nod2D                                 
     nini=icestiff%rowptr(k)-icestiff%rowptr(1)+1    
     nend=icestiff%rowptr(k+1)-icestiff%rowptr(1)    
     icestiff%colind(nini:nend)= &
          nghbr_nod2D(k)%addresses
  end do

  ! ==================================
  ! addresses are now in local numbering. We need to make them contiguous
  do k=1,icestiff%rowptr(myDim_nod2D+1)-icestiff%rowptr(1)    
     icestiff%colind(k)=myCList_nod2D(icestiff%colind(k))
  end do

end subroutine icestiff_matrix
!
!----------------------------------------------------------------
!
subroutine icestiff_fill_tg
  ! Fill in ice stiffness(mass) matrix for TG/FCT scheme
  use o_MESH
  use o_ELEMENTS
  use o_PARAM
  use i_array
  use g_parfe
  use g_config
  implicit none
  !
  integer                            :: i, q, row, col, ipos, offset
  integer                            :: elem, elnodes(3)
  real(kind=8)                       :: vol, inv12
  !
  inv12=1.0_8/12.0_8
  icestiff%values=0.0
  allocate(ice_lump(myDim_nod2d+eDim_nod2D))              
  ice_lump=0.0_8
  !
  do elem=1,myDim_elem2D                                 

     elnodes=elem2D_nodes(:,elem)
     vol=voltriangle(elem)*dt_inv*inv12
     do i=1,3             ! all rows into which the element elem could contribute
        row=elnodes(i)
	if(row>myDim_nod2D) cycle                         
        do q=1,nghbr_nod2D(row)%nmb
           col_pos(nghbr_nod2D(row)%addresses(q))=q
        end do
        offset=icestiff%rowptr(row)-icestiff%rowptr(1)    
        do q=1,3          ! all columns 
           col=elnodes(q)
           ipos = offset+col_pos(col)
	   icestiff%values(ipos)=icestiff%values(ipos)+vol
           if (q==i) then
              icestiff%values(ipos)=icestiff%values(ipos)+vol
           end if
        end do
        ice_lump(row)=ice_lump(row) + vol*4.0

     end do
  end do
  call com_2D(ice_lump)

end subroutine icestiff_fill_tg
!
!------------------------------------------------------------------------
!
#ifndef use_ice_gls
subroutine ice_rhs_tg
  use o_mesh
  use o_elements
  use o_param
  use i_dyn_parms
  use i_array
  use g_PARFE
  use g_config
  implicit none 
  !
  integer                      :: m, i, row, elnodes(3), elem, edgnod(2)
  real(kind=8)                 :: dx(3), dy(3), vol, aux, adv
  real(kind=8)                 :: uvel(3), vvel(3), usum, vsum, um, vm
  real(kind=8)                 :: m_el(3), a_el(3), ms_el(3)
  real(kind=8)                 :: mdiffx, mdiffy, adiffx, adiffy 
  real(kind=8)                 :: msdiffx, msdiffy, m_sum, a_sum, ms_sum
  real(kind=8)                 :: um_sum, vm_sum, ua_sum, va_sum
  real(kind=8)                 :: ums_sum, vms_sum, div
  real(kind=8)                 :: inv2, inv3, inv6, inv12
  real(kind=8)                 :: nv(2), uedg(2), vedg(2)
  real(kind=8)                 :: m_edg(2), a_edg(2), ms_edg(2)
  real(kind=8)                 :: aux1, aux2, aux3

  inv2=0.5_8
  inv3=1.0_8/3.0_8
  inv12=1.0_8/12.0_8
  do row=1, myDim_nod2d          

     rhs_m(row)=0.
     rhs_a(row)=0.
     rhs_ms(row)=0.
  enddo
  !
  do elem=1,myDim_elem2d         
     elnodes=elem2D_nodes(:,elem)
     dx=bafux_2d(:, elem)
     dy=bafuy_2d(:, elem)
     vol=voltriangle(elem)
     uvel=u_ice(elnodes)
     vvel=v_ice(elnodes)
     usum=sum(uvel)
     vsum=sum(vvel)
     um=inv3*usum
     vm=inv3*vsum	
     m_el=m_ice(elnodes)
     a_el=a_ice(elnodes)
     ms_el=m_snow(elnodes)
     mdiffx=sum(dx*m_el)
     mdiffy=sum(dy*m_el)
     adiffx=sum(dx*a_el)
     adiffy=sum(dy*a_el)
     msdiffx=sum(dx*ms_el)
     msdiffy=sum(dy*ms_el)
     m_sum=sum(m_el)
     a_sum=sum(a_el)
     ms_sum=sum(ms_el)
     um_sum=sum(m_el*uvel)
     vm_sum=sum(m_el*vvel)
     ua_sum=sum(a_el*uvel)
     va_sum=sum(a_el*vvel)
     ums_sum=sum(ms_el*uvel)
     vms_sum=sum(ms_el*vvel)
     div=sum(dx*uvel+dy*vvel)
     !
     !assembling over nodes
     do i=1,3 
        row=elnodes(i)

        !diffusion
        rhs_m(row)=rhs_m(row)-Kh_ice*(mdiffx*dx(i)+mdiffy*dy(i))*vol
        rhs_a(row)=rhs_a(row)-Kh_ice*(adiffx*dx(i)+adiffy*dy(i))*vol
        rhs_ms(row)=rhs_ms(row)-Kh_ice*(msdiffx*dx(i)+msdiffy*dy(i))*vol

        !advection
        adv=(dx(i)*(usum*m_sum+um_sum)+dy(i)*(vsum*m_sum+vm_sum))*inv12
        rhs_m(row)=rhs_m(row)+adv*vol
        adv=(dx(i)*(usum*a_sum+ua_sum)+dy(i)*(vsum*a_sum+va_sum))*inv12
        rhs_a(row)=rhs_a(row)+adv*vol
        adv=(dx(i)*(usum*ms_sum+ums_sum)+dy(i)*(vsum*ms_sum+vms_sum))*inv12
        rhs_ms(row)=rhs_ms(row)+adv*vol

        !TG stabilization
        aux=(um*dx(i)+vm*dy(i))*dt*inv2*vol
        adv=um*mdiffx+vm*mdiffy+div*m_sum*inv3
        rhs_m(row)=rhs_m(row)-aux*adv
        adv=um*adiffx+vm*adiffy+div*a_sum*inv3
        rhs_a(row)=rhs_a(row)-aux*adv
        adv=um*msdiffx+vm*msdiffy+div*ms_sum*inv3
        rhs_ms(row)=rhs_ms(row)-aux*adv
     end do
  end do
  !
#ifdef use_opbnd_restoring 
  inv6=1.0_8/6.0_8
  do m=1,nmbr_opbnd_edg
     edgnod=opbnd_edg(m,1:2)
     nv=opbnd_edg_nv(m,1:2)
     vol=opbnd_edg_nv(m,3)
     uedg=u_ice(edgnod)
     vedg=v_ice(edgnod)
     vedg=uedg*nv(1)+vedg*nv(2)
     if(vedg(1)<0.) vedg(1)=0.     ! assume no ice is entering the model domain
     if(vedg(2)<0.) vedg(2)=0.

     m_edg=m_ice(edgnod)
     a_edg=a_ice(edgnod)
     ms_edg=m_snow(edgnod)

     aux1=(vedg(1)*m_edg(2)+vedg(2)*m_edg(1)+sum(vedg*m_edg))*inv12
     aux2=(vedg(1)*a_edg(2)+vedg(2)*a_edg(1)+sum(vedg*a_edg))*inv12
     aux3=(vedg(1)*ms_edg(2)+vedg(2)*ms_edg(1)+sum(vedg*ms_edg))*inv12

     do i=1,2
        row=edgnod(i)
        rhs_m(row)=rhs_m(row)-(aux1+vedg(i)*m_edg(i)*inv6)*vol
        rhs_a(row)=rhs_a(row)-(aux2+vedg(i)*a_edg(i)*inv6)*vol
        rhs_ms(row)=rhs_ms(row)-(aux3+vedg(i)*ms_edg(i)*inv6)*vol
     end do
  end do
#endif

end subroutine ice_rhs_tg
#endif
!
!=========================================================================
!
#ifdef use_ice_gls
subroutine icestiff_fill_gls
  use o_PARAM
  use o_MESH
  use o_ELEMENTS
  use i_dyn_parms
  use i_array
  use g_PARFE
  use g_config
  implicit none

  integer                           :: m, i, q, ipos, offset, edgnod(2)
  integer                           :: col, col2, row, elem, elnodes(3)
  real(kind=8)                      :: dx(3), dy(3), entries(3), entries_t(3)
  real(kind=8)                      :: v, uel(3), vel(3), usum, vsum, um, vm
  real(kind=8)                      :: div, aux, hh, v_nabla_psi, adv, diff
  real(kind=8)                      :: ept, ept_invdt, inv12, inv6, inv3, inv2
  real(kind=8)                      :: nv(2), uedg(2), vedg(2)
  real(kind=8)                      :: m_edg(2), a_edg(2), ms_edg(2)
  real(kind=8)                      :: aux1, aux2, aux3

  inv12=1.0_8/12.0_8
  inv3=1.0_8/3.0_8
  inv2=0.5_8
  do row=1,myDim_nod2d             
     rhs_m(row)=0.
     rhs_a(row)=0.
     rhs_ms(row)=0.
     col=icestiff%rowptr(row)-icestiff%rowptr(1)+1   
     col2=icestiff%rowptr(row+1)-icestiff%rowptr(1)  
     icestiff%values(col:col2)=0.0
  enddo

  do elem=1, myDim_elem2d           
     elnodes=elem2D_nodes(:,elem)
     v=voltriangle(elem)
     dx=bafux_2D(:,elem)
     dy=bafuy_2D(:,elem)
     uel=u_ice(elnodes)
     vel=v_ice(elnodes)
     usum=sum(uel)
     vsum=sum(vel)
     um=usum*inv3
     vm=vsum*inv3
     div=sum(uel*dx + vel*dy)
     aux=v*dt_inv*inv12   
     hh=sqrt(v)
     ept=0.1*dt_inv + (abs(um)+abs(vm))/hh + 2.0*Kh_ice/(hh*hh)
     ept=v/ept
     ept_invdt=ept*dt_inv

     do i=1,3             ! all rows into which the element elem could contribute
        row=elnodes(i)
        if(row>myDim_nod2D) cycle     
        do q=1,nghbr_nod2D(row)%nmb
           col_pos(nghbr_nod2D(row)%addresses(q))=q
        end do

        v_nabla_psi=um*dx(i)+vm*dy(i) 

        entries=0.0
        do q=1,3          ! all columns 

           diff=Kh_ice*(dx(i)*dx(q)+dy(i)*dy(q))                ! pure numerical things
           adv=-(dx(i)*(usum+uel(q))+dy(i)*(vsum+vel(q)))*inv12        
	   entries(q)=(adv+diff)*v                              ! advection + numerical diff

	   adv= um*dx(q) + vm*dy(q) + div*inv3 
           entries(q)=entries(q)+ ept*v_nabla_psi*adv	        ! stabilization

           col=elnodes(q)
           rhs_m(row)=rhs_m(row)-entries(q)*m_ice(col)		
           rhs_a(row)=rhs_a(row)-entries(q)*a_ice(col)		
           rhs_ms(row)=rhs_ms(row)-entries(q)*m_snow(col)       ! fill rhs

           entries_t(q)=aux                             	! mass matrix

           entries_t(q)=entries_t(q)+ept_invdt*v_nabla_psi*inv3	! stabilization, t          

        end do

        entries_t(i)=entries_t(i)+aux                    	! completes mass matrix 
        entries=inv2*entries+entries_t       

        ! put the entries to the appropriate place
        offset=icestiff%rowptr(row)-icestiff%rowptr(1)    
        do q=1,3        ! all columns 
           col=elnodes(q)
           ipos = offset+col_pos(col)
           icestiff%values(ipos)=icestiff%values(ipos)+entries(q)
        end do
     end do
  end do

#ifdef use_opbnd_restoring 
  inv6=1.0_8/6.0_8
  do m=1,nmbr_opbnd_edg
     edgnod=opbnd_edg(m,1:2)
     nv=opbnd_edg_nv(m,1:2)
     v=opbnd_edg_nv(m,3)
     uedg=u_ice(edgnod)
     vedg=v_ice(edgnod)
     vedg=uedg*nv(1)+vedg*nv(2)
     if(vedg(1)<0.) vedg(1)=0.     ! assume no ice coming into the model domain
     if(vedg(2)<0.) vedg(2)=0.

     m_edg=m_ice(edgnod)
     a_edg=a_ice(edgnod)
     ms_edg=m_snow(edgnod)

     aux1=(vedg(1)*m_edg(2)+vedg(2)*m_edg(1)+sum(vedg*m_edg))*inv12
     aux2=(vedg(1)*a_edg(2)+vedg(2)*a_edg(1)+sum(vedg*a_edg))*inv12
     aux3=(vedg(1)*ms_edg(2)+vedg(2)*ms_edg(1)+sum(vedg*ms_edg))*inv12

     do i=1,2
        row=edgnod(i)
        rhs_m(row)=rhs_m(row)-(aux1+vedg(i)*m_edg(i)*inv6)*v
        rhs_a(row)=rhs_a(row)-(aux2+vedg(i)*a_edg(i)*inv6)*v
        rhs_ms(row)=rhs_ms(row)-(aux3+vedg(i)*ms_edg(i)*inv6)*v
     end do
  end do
#endif

end subroutine icestiff_fill_gls
#endif
!
!----------------------------------------------------------------------------
!
subroutine ice_solve
  use o_mesh
  use i_array
  use o_array
  use g_parfe
  use i_dyn_parms
  use g_config
  implicit none
  !
  integer                                 :: n, row, clo, clo2, cn, location(100)
  real(kind=8)                            :: rhs_new
  real(kind=8), allocatable               :: auxarray(:,:)

  allocate(auxarray(3,myDim_nod2d))

  !the first approximation
  do row=1,myDim_nod2D                  
     dm_ice(row)=rhs_m(row)/ice_lump(row)
     da_ice(row)=rhs_a(row)/ice_lump(row)
     dm_snow(row)=rhs_ms(row)/ice_lump(row)
  end do

  call com_2D(dm_ice)
  call com_2D(da_ice)
  call com_2D(dm_snow)

  !iterate 
  do n=1,num_iter_solve_ice-1
     do row=1,myDim_nod2D               
        clo=icestiff%rowptr(row)-icestiff%rowptr(1)+1  
        clo2=icestiff%rowptr(row+1)-icestiff%rowptr(1) 
        cn=clo2-clo+1
        location(1:cn)=nghbr_nod2D(row)%addresses      
        rhs_new=rhs_m(row) - sum(icestiff%values(clo:clo2)*dm_ice(location(1:cn)))
        auxarray(1,row)=dm_ice(row)+rhs_new/ice_lump(row)    
        rhs_new=rhs_a(row) - sum(icestiff%values(clo:clo2)*da_ice(location(1:cn)))
        auxarray(2,row)=da_ice(row)+rhs_new/ice_lump(row)     
        rhs_new=rhs_ms(row) - sum(icestiff%values(clo:clo2)*dm_snow(location(1:cn)))
        auxarray(3,row)=dm_snow(row)+rhs_new/ice_lump(row)  
     end do

     do row=1,myDim_nod2D             
        dm_ice(row)=auxarray(1,row)      
        da_ice(row)=auxarray(2,row)      
	dm_snow(row)=auxarray(3,row)    
     end do
     call com_2D(dm_ice)
     call com_2D(da_ice)
     call com_2D(dm_snow)
  end do

  deallocate(auxarray)
end subroutine ice_solve
