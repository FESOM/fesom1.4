! Set up and initialized ocean matrices

subroutine ocean_matrices_setup
  !driving routine for setting up matrices
  
  use o_param
  use o_mesh
  use g_parfe
  implicit none

  call set_coriolis_param      ! set f before assembling ssh matrix

  call sshstiff_matrix

#ifdef use_non_hydrostatic
  call nhpstiff_matrix  
#endif

  call tsstiff_construct

#ifndef use_tracer_gls
  call tsstiff_fill_tg
#endif

  call uvstiff_matrix

#ifndef use_non_hydrostatic
  call build_wpot_matrix
#endif

  if(mype==0) write(*,*) 'Ocean matrices have been set up'

end subroutine ocean_matrices_setup
!
!----------------------------------------------------------------------
!
subroutine uvstiff_matrix
  ! Stiffness matrix for 3D velocities (u and v entries are not coupled).
  ! mass matrix does not contain dt_inv 
  !
  ! Coded by Sergey Danilov and Qiang Wang
  ! Reviewed by ??
  !--------------------------------------------------------
  
  use o_param
  use o_array
  use o_mesh
  use o_elements
  use o_matrices
  use g_config
  use g_parfe
  implicit none

  integer                            :: k
  integer                            :: row, col, elnodes(4)
  integer                            :: i, q, elem, ind
  integer                            :: ipos, offset, is, ie
  real(kind=8)                       :: vol, inv20
  integer, allocatable              :: pnza(:), rpnza(:) 
  integer                           :: nini, nend

  ! a)
  uvstiff%dim=nod3D
  allocate(uvstiff%rowptr(myDim_nod3D+1))        
  uvstiff%rowptr(1)=1
  do k=1,myDim_nod3D                            
     uvstiff%rowptr(k+1)=uvstiff%rowptr(k)+nghbr_nod3D(k)%nmb
  end do
  uvstiff%nza=uvstiff%rowptr(myDim_nod3D+1)-1    
  ! b)
  allocate(uvstiff%colind(uvstiff%nza))
  allocate(uvstiff%values(uvstiff%nza))
  uvstiff%values=0.0

  ! =================================            
  ! Exchange between PE to know their 
  ! numbers of non-zero entries (nza):
  ! =================================     
  allocate(pnza(npes), rpnza(npes))
  pnza(1:npes)=0
  pnza(mype+1)=uvstiff%nza
  call MPI_Barrier(MPI_COMM_FESOM,MPIerr)
  call MPI_AllREDUCE( pnza, rpnza, &
       npes, MPI_INTEGER,MPI_SUM, &
       MPI_COMM_FESOM, MPIerr)

  if (mype==0) then
     i=0
  else
     i=sum(rpnza(1:mype))  ! This is offset for mype 
  end if
  uvstiff%rowptr=uvstiff%rowptr+i       

  ! =================================
  ! replace local nza with a global one
  ! =================================
  uvstiff%nza=sum(rpnza(1:npes))                  
  deallocate(rpnza, pnza)

  ! c)
  do k=1,myDim_nod3D                                 
     nini=uvstiff%rowptr(k)-uvstiff%rowptr(1)+1      
     nend=uvstiff%rowptr(k+1)-uvstiff%rowptr(1)      
     uvstiff%colind(nini:nend)= &                    
          nghbr_nod3D(k)%addresses
  end do

  ! ==================================
  ! addresses are now in local numbering. We need to make them contiguous
  do k=1,uvstiff%rowptr(myDim_nod3D+1)-uvstiff%rowptr(1)   
     uvstiff%colind(k)=myCList_nod3D(uvstiff%colind(k))          
  end do
  ! ==================================

  ! d) fill in
  allocate(uv_lump(myDim_nod3D+eDim_nod3D))           
  uv_lump=0.
  inv20=1.0_8/20.0_8

  do elem=1,myDim_elem3d                         
     elnodes=elem3D_nodes(:,elem)
     vol=voltetra(elem)*inv20
     do i=1,4             ! all rows into which the element elem could contribute
        row=elnodes(i)
	if(row>myDim_nod3D) cycle                   
        do q=1,nghbr_nod3D(row)%nmb
           col_pos(nghbr_nod3D(row)%addresses(q))=q
        end do
        offset=uvstiff%rowptr(row)-uvstiff%rowptr(1)  
        do q=1,4          ! all columns 
           col=elnodes(q)
           ipos = offset+col_pos(col)
	   uvstiff%values(ipos)=uvstiff%values(ipos)+vol
           if (q==i) then
              uvstiff%values(ipos)=uvstiff%values(ipos)+vol
           end if
        end do
        uv_lump(row)=uv_lump(row) + vol*5.0
     end do
  end do

  !boundary condition (no slip bc.)
  do i=1,myDim_nod3d                              
     ind = index_nod3D(i)
     if (ind==11 .or. ind==21 .or. ind==31) then
        is=uvstiff%rowptr(i)-uvstiff%rowptr(1)+1  
        ie=uvstiff%rowptr(i+1)-uvstiff%rowptr(1)   
        where(uvstiff%colind(is:ie)==myCList_nod3D(i))   
           uvstiff%values(is:ie)=1.0
        elsewhere
           uvstiff%values(is:ie)=0.0
        end where
     end if
  end do

  call com_3D(uv_lump)  

end subroutine uvstiff_matrix
!
!========================================================================
!
subroutine sshstiff_matrix
  ! Stiffness matrix for ssh (surface pressure)
  !
  ! Coded by Sergey Danilov and Qiang Wang
  ! Reviewed by ??
  !--------------------------------------------------------
  
  use o_PARAM
  use o_array
  use o_MESH
  use o_ELEMENTS
  use o_MATRICES
  use g_config
  use G_parfe
  implicit none
  integer                           :: i, k, q, ipos, is, ie, row, col
  integer                           :: elem, elem2, elnodes2(3), mn(3)
  real(kind=8)                      :: vol, dx(3), dy(3), val, aux, tri_v(3)
  real(kind=8)                      :: cori_p, dparam,  beta, gamma
  integer, allocatable              :: pnza(:), rpnza(:) 
  integer                           :: nini, nend

  ! a)
  sshstiff%dim=nod2D
  allocate(sshstiff%rowptr(myDim_nod2D+1))         
  sshstiff%rowptr(1)=1              ! This will be updated as
  ! contiguous numbering is required 
  do k=1,myDim_nod2D                                
     sshstiff%rowptr(k+1)=sshstiff%rowptr(k)+nghbr_nod2D(k)%nmb
  end do
  sshstiff%nza=sshstiff%rowptr(myDim_nod2D+1)-1      

  ! b)
  allocate(sshstiff%colind(sshstiff%nza))
  allocate(sshstiff%values(sshstiff%nza))
  sshstiff%values=0.0

  ! =================================                  
  ! Now we need to exchange between PE to know their 
  ! numbers of non-zero entries (nza):
  ! =================================     
  allocate(pnza(npes),rpnza(npes))
  pnza=0
  pnza(mype+1)=sshstiff%nza
  call MPI_Barrier(MPI_COMM_FESOM,MPIerr)
  call MPI_AllREDUCE( pnza, rpnza, &
       npes, MPI_INTEGER,MPI_SUM, &
       MPI_COMM_FESOM, MPIerr)

  if (mype==0) then
     i=0
  else
     i=sum(rpnza(1:mype))  ! This is offset for mype 
  end if
  sshstiff%rowptr=sshstiff%rowptr+i       

  ! =================================
  ! replace local nza with a global one
  ! =================================
  sshstiff%nza=sum(rpnza(1:npes))                    
  deallocate(rpnza, pnza)

  ! c) FIND colind 
  do k=1,myDim_nod2D
     nini=sshstiff%rowptr(k)-sshstiff%rowptr(1)+1     
     nend=sshstiff%rowptr(k+1)-sshstiff%rowptr(1)     
     sshstiff%colind(nini:nend)= &
          nghbr_nod2D(k)%addresses
  end do
  ! ==================================
  ! colind are now in local numbering. We need to make them contiguous
  do k=1,sshstiff%rowptr(myDim_nod2D+1)-sshstiff%rowptr(1)  
     sshstiff%colind(k)=myCList_nod2D(sshstiff%colind(k))   
  end do

  ! ==================================
  ! d) fill in

  ! g*dt*\int_{-H}^0(\nabla\phi_i * \nabla\phi_j) d\Omega
  do elem=1, myDim_elem3d                         
     elem2 = elem2D_corresp_to_elem3D(elem)
     elnodes2=elem2d_nodes(:,elem2)
     vol=voltetra(elem)
     dx=bafux_2d(:,elem2)
     dy=bafuy_2d(:,elem2)
     if(use_cori_semi) then
        !coriolis parameter
        cori_p=coriolis_param_elem2d(elem2)
        !terms for semi-implicit schemes
        dparam=dt_inv**2 + alpha_trapez**2*cori_p**2
        beta=dt_inv/dparam
        gamma=cori_p*alpha_trapez/dparam
     endif
#ifdef use_semiimplicit_scheme
     vol=vol*theta_ssh*theta_vel
#endif
     do i=1,3
        row=elnodes2(i)
	if (row>myDim_nod2D) cycle                  
        do q=1, nghbr_nod2D(row)%nmb
           col_pos(nghbr_nod2D(row)%addresses(q)) = q
        enddo
        do q=1,3
           col=elnodes2(q)
           if(use_cori_semi) then
              val=(dx(i)*dx(q)+dy(i)*dy(q))*vol*g*beta
              val=val+(dx(i)*dy(q)-dy(i)*dx(q))*vol*g*gamma
           else
              val = (dx(i)*dx(q) + dy(i)*dy(q))*vol*g*dt
           endif
           ipos=sshstiff%rowptr(row)+col_pos(col)-sshstiff%rowptr(1) 
           sshstiff%values(ipos)=sshstiff%values(ipos)+val
        end do
     end do
  end do

  ! Mass matrix part
  do i=1, myDim_nod2d                                  
     is=sshstiff%rowptr(i)-sshstiff%rowptr(1)+1       
     ie=sshstiff%rowptr(i+1)-sshstiff%rowptr(1)       
     do q=1, nod_in_elem2D(i)%nmb
        elem2=nod_in_elem2D(i)%addresses(q)
        vol=voltriangle(elem2)
        elnodes2=elem2D_nodes(:, elem2)
        do k=1, 3
           col=elnodes2(k)
           aux=1./12.
           if (col==i) aux=1./6.
           row=myCList_nod2D(col)
           where(sshstiff%colind(is:ie)==row)     
              sshstiff%values(is:ie)=sshstiff%values(is:ie)+aux*vol*dt_inv
           end where
        end do
     end do
  end do

  ! open boundary
#ifdef use_opbnd_tide
  if(nmbr_opbnd_tri==0) return                           
  if(trim(tide_opbnd_type)=='ssh') then
     do i=1,nmbr_opbnd_t2d
        row=opbnd_n2d(i)
        is=sshstiff%rowptr(row)-sshstiff%rowptr(1)+1     
        ie=sshstiff%rowptr(row+1)-sshstiff%rowptr(1)     
	q=myCList_nod2D(row)
        where(sshstiff%colind(is:ie)==q)                 
           sshstiff%values(is:ie)=1.0
        elsewhere
           sshstiff%values(is:ie)=0.0
        end where
     end do
  elseif(trim(tide_opbnd_type)=='Flather') then
     do elem2=1, nmbr_opbnd_tri
        elnodes2=opbnd_tri(elem2,1:3)
        elnodes2=nod2d_corresp_to_nod3d(elnodes2)    
        mn=mapping_opbnd_n2d(elnodes2)
        vol=opbnd_nv(elem2,4)
        tri_v=sqrt(g/opbnd_dep(mn)) 
        do i=1,3        
           row=elnodes2(i)
           if(row>myDim_nod2d) cycle
           is=sshstiff%rowptr(row)-sshstiff%rowptr(1)+1        
           ie=sshstiff%rowptr(row+1)-sshstiff%rowptr(1)         
           do k=1, 3
              col=elnodes2(k)
              aux=1./12.
              if (row==col) aux=1./6.
	      q=myCList_nod2D(col)
              where(sshstiff%colind(is:ie)==q)               
                 sshstiff%values(is:ie)=sshstiff%values(is:ie)+aux*vol*tri_v(k)
              end where
           end do
        end do
     end do
  end if
#endif

end subroutine sshstiff_matrix
!
!============================================================================
!
#ifdef use_non_hydrostatic
subroutine nhpstiff_matrix
  ! Stiffness matrix for non-hydrostatic pressure
  ! This is only a test version!  
  !
  ! Coded by Qiang Wang
  ! Reviewed by ??
  !--------------------------------------------------------
  
  use o_PARAM
  use o_array
  use o_MESH
  use o_ELEMENTS
  use o_MATRICES
  use g_config
  use g_parfe
  implicit none
  integer                           :: i, k, ipos, q
  integer                           :: row, col, elem, elnodes(4), elem2
  real(kind=8)                      :: vol, dx(4), dy(4), dz(4), val
  real(kind=8)                      :: cori_p, dparam, beta, gamma
  integer, allocatable              :: pnza(:), rpnza(:) 
  integer                           :: nini, nend
  !
  ! a)
  nhpstiff%dim=nod3D
  allocate(nhpstiff%rowptr(myDim_nod3D+1))          
  nhpstiff%rowptr(1)=1
  do k=1,myDim_nod3D                                 
     nhpstiff%rowptr(k+1)=nhpstiff%rowptr(k)+nghbr_nod3D(k)%nmb
  end do
  nhpstiff%nza=nhpstiff%rowptr(myDim_nod3D+1)-1      
  ! b)
  allocate(nhpstiff%colind(nhpstiff%nza))
  allocate(nhpstiff%values(nhpstiff%nza))
  nhpstiff%values=0.0
  ! =================================                 
  ! Exchange between PE to know their 
  ! numbers of non-zero entries (nza):
  ! =================================     
  allocate(pnza(npes), rpnza(npes))
  pnza(1:npes)=0
  pnza(mype+1)=nhpstiff%nza
  call MPI_Barrier(MPI_COMM_FESOM,MPIerr)
  call MPI_AllREDUCE( pnza, rpnza, &
       npes, MPI_INTEGER,MPI_SUM, &
       MPI_COMM_FESOM, MPIerr)

  if (mype==0) then
     i=0
  else
     i=sum(rpnza(1:mype))  ! This is offset for mype 
  end if
  nhpstiff%rowptr=nhpstiff%rowptr+i       

  ! =================================
  ! replace local nza with a global one
  ! =================================
  nhpstiff%nza=sum(rpnza(1:npes))
  deallocate(rpnza, pnza)                           

  ! c)
  do k=1,myDim_nod3D                                
     nini=nhpstiff%rowptr(k)-nhpstiff%rowptr(1)+1     
     nend=nhpstiff%rowptr(k+1)-nhpstiff%rowptr(1)     
     nhpstiff%colind(nini:nend)= &                    
          nghbr_nod3D(k)%addresses
  end do
  ! ==================================
  ! addresses are now in local numbering. We need to make them contiguous
  do k=1,nhpstiff%rowptr(myDim_nod3D+1)-nhpstiff%rowptr(1) 
     nhpstiff%colind(k)=myCList_nod3D(nhpstiff%colind(k))        
  end do
  ! ==================================


  ! d) fill in 
  do elem=1,myDim_elem3d                  
     elnodes=elem3d_nodes(:,elem)
     vol=voltetra(elem)
     dx=bafux_3d(:,elem)
     dy=bafuy_3d(:,elem)
     dz=bafuz_3d(:,elem)
     if(use_cori_semi) then
        elem2=elem2d_corresp_to_elem3d(elem)
        !coriolis parameter
        cori_p=coriolis_param_elem2d(elem2)
        !terms for semi-implicit schemes
        dparam=dt_inv**2 + alpha_trapez**2*cori_p**2
        beta=dt_inv/dparam
        gamma=cori_p*alpha_trapez/dparam
     endif
     do i=1,4
        row=elnodes(i)
	if (row>myDim_nod3D) cycle                        
        do q=1, nghbr_nod3D(row)%nmb
           col_pos(nghbr_nod3D(row)%addresses(q)) = q
        enddo
        do q=1,4
           col=elnodes(q)
           if(use_cori_semi) then
              val=(dx(i)*dx(q)+dy(i)*dy(q))*beta
              val=val+(dx(i)*dy(q)-dy(i)*dx(q))*gamma
              val=(val+dz(i)*dz(q)*dt)*vol
           else
              val=(dx(i)*dx(q)+dy(i)*dy(q)+dz(i)*dz(q))*vol*dt
           endif
           ipos=nhpstiff%rowptr(row)+col_pos(col)-nhpstiff%rowptr(1)  
           nhpstiff%values(ipos)=nhpstiff%values(ipos)+val
        end do
     end do
  end do
end subroutine nhpstiff_matrix
#endif
!
!============================================================================
!
#ifndef use_non_hydrostatic
subroutine build_wpot_matrix
  ! w potential matrix
  ! Taking into account the property of vertical dirivatives in tetrahedra:
  ! only a three-diagonal matrix is required.
  !
  ! Coded by Qiang Wang and Sergey Danilov
  ! Reviewed by ??
  !--------------------------------------------------------
  
  use o_param
  use o_mesh
  use o_elements
  use o_matrices
  use g_config
  use g_parfe
  implicit none

  integer            :: n2, nlay, k, i, q, p
  integer            :: row, elem, elnodes(4)
  integer            :: nodup, nodlo
  real(kind=8)       :: vol
  real(kind=8)       :: a(max_num_layers), b(max_num_layers)
  real(kind=8)       :: c(max_num_layers)

  allocate(wpot_matrix(3,max_num_layers,myDim_nod2d))
  wpot_matrix=0.0

  do n2=1,myDim_nod2d
     nlay=num_layers_below_nod2d(n2)+1

     ! assemble three diagonal matrix
     do k=1,nlay
        a(k)=0.
        b(k)=0.
        c(k)=0.

        row=nod3d_below_nod2d(k,n2)
        if(k>1) nodup=nod3d_below_nod2d(k-1,n2)
        if(k<nlay) nodlo=nod3d_below_nod2d(k+1,n2)
        do i=1,nod_in_elem3D(row)%nmb
           elem=nod_in_elem3D(row)%addresses(i)
           elnodes=elem3D_nodes(:,elem)  
           vol=voltetra(elem)

           do q=1,4
              if(elnodes(q)==row) then
                 p=q
                 exit
              end if
           end do

           !first entry
           if(k>1) then
              do q=1,4
                 if(elnodes(q)==nodup) then
                    a(k)=a(k) + bafuz_3d(p,elem) * bafuz_3d(q,elem) * vol
                    exit
                 end if
              end do
           end if

           !second entry
           b(k)=b(k) + bafuz_3d(p,elem) * bafuz_3d(p,elem) * vol

           !third entry
           if(k<nlay) then
              do q=1,4
                 if(elnodes(q)==nodlo) then
                    c(k)=c(k) + bafuz_3d(p,elem) * bafuz_3d(q,elem)  * vol
                    exit
                 end if
              end do
           end if

        end do  !i
     end do  !k

     wpot_matrix(1,1:nlay,n2)=a(1:nlay)      
     wpot_matrix(2,1:nlay,n2)=b(1:nlay)    
     wpot_matrix(3,1:nlay,n2)=c(1:nlay)    

  end do  !m

end subroutine build_wpot_matrix
#endif
!
!=========================================================================
!
subroutine tsstiff_construct
  !T/S stiffness matrix
  !
  ! Coded by Sergey Danilov and Qiang Wang
  ! Reviewed by ??
  !--------------------------------------------------------
  
  use o_MATRICES
  use o_MESH
  use o_ELEMENTS
  use o_PARAM
  use g_config
  use g_PARFE
  implicit none
  !
  integer                           :: k
  integer, allocatable              :: pnza(:), rpnza(:) 
  integer                           :: nini, nend

  ! a)
  tsstiff%dim=nod3D
  allocate(tsstiff%rowptr(myDim_nod3D+1))     
  tsstiff%rowptr(1)=1
  do k=1,myDim_nod3D                          
     tsstiff%rowptr(k+1)=tsstiff%rowptr(k)+nghbr_nod3D(k)%nmb
  end do
  tsstiff%nza=tsstiff%rowptr(myDim_nod3D+1)-1   

  ! b)
  allocate(tsstiff%colind(tsstiff%nza))
  allocate(tsstiff%values(tsstiff%nza))
  tsstiff%values=0.0
  ! =================================           
  ! Exchange between PE to know their 
  ! numbers of non-zero entries (nza):
  ! =================================     
  allocate(pnza(npes),rpnza(npes))
  pnza(1:npes)=0
  pnza(mype+1)=tsstiff%nza

  call MPI_Barrier(MPI_COMM_FESOM,MPIerr)
  call MPI_AllREDUCE( pnza, rpnza, &
       npes, MPI_INTEGER,MPI_SUM, &
       MPI_COMM_FESOM, MPIerr)

  if (mype==0) then
     k=0
  else
     k=sum(rpnza(1:mype))  ! This is offset for mype 
  end if
  tsstiff%rowptr=tsstiff%rowptr+k       

  ! =================================
  ! replace local nza with a global one
  ! =================================
  tsstiff%nza=sum(rpnza(1:npes))              
  deallocate(rpnza,pnza)
  ! c)
  ! ===== FIND irowind ======
  do k=1,myDim_nod3D                            
     nini=tsstiff%rowptr(k)-tsstiff%rowptr(1)+1 
     nend=tsstiff%rowptr(k+1)-tsstiff%rowptr(1) 
     tsstiff%colind(nini:nend)= &                
          nghbr_nod3D(k)%addresses
  end do
  ! ==================================
  ! addresses are now in local numbering. We need to make them contiguous
  do k=1,tsstiff%rowptr(myDim_nod3D+1)-tsstiff%rowptr(1) 
     tsstiff%colind(k)=myCList_nod3D(tsstiff%colind(k))        
  end do
  ! ==================================

end subroutine tsstiff_construct
!
!=========================================================================
!
#ifdef use_tracer_gls
subroutine tsstiff_fill_gls
  ! Fill in T/S stiffness matrix for GLS scheme
  ! 1) Fills in  advective, diffusive and time contributions
  ! 2) Computes RHSs for temperature and salinity
  !
  ! Coded by Sergey Danilov and Qiang Wang
  ! Reviewed by ??
  !--------------------------------------------------------
  
  use o_MATRICES
  use o_MESH
  use o_ELEMENTS
  use o_PARAM
  use o_array
  use g_config
  use g_forcing_arrays
  use g_PARFE
  implicit none 

  integer                       :: k, i, q, elem, elem2, elem_type
  integer                       :: elnodes(4), is, ie, lay
  integer                       :: row, col, col2, ipos, offset, ntr
  real(kind=8)                  :: uvel(4), vvel(4), usum, vsum
  real(kind=8)                  :: u2d, v2d, um, vm, wm
  real(kind=8)                  :: dx(4), dy(4), dz(4), tra_elem(4,num_tracer)
  real(kind=8)                  :: auxf, advection, diffusion
  real(kind=8)                  :: vnabla_t, ept, ept_dt_inv
  real(kind=8)                  :: v, vtri, hh, hv, entries(4), entries_t(4)
  real(kind=8)                  :: Kh, Kv_el, r_coeff
  real(kind=8)                  :: inv2, inv4, inv5, inv20
  real(kind=8)                  :: fc, aux1, aux2, dif(4)
  real(kind=8)                  :: rotate_coe, temp_coe, temp_coe2
  real(kind=8)                  :: swr_conv
  !
  !variables used for Redi/GM scheme 
  real(kind=8)                  :: K_GM, Kh_diag, Kh_other
  real(kind=8)		        :: S(3), fcn1,fcn2, lambd, Z_mean, depth_scale
  real(kind=8)        		:: S_d, c_speed
  data S_d, c_speed /1.0e-3, 2.0/
  !
  real(kind=8)                  :: dparam, beta, gamma
#ifndef use_non_hydrostatic
  integer                       :: elnodes2(3)
  real(kind=8)                  :: vc(3)
#else
  integer                       :: elnodes23(4)
  real(kind=8)                  :: vc(4), wsum, wvel(4), w2d
#endif
#ifdef use_fullfreesurf
  integer                       :: n_el
  real(kind=8)                  :: wmove(4), wmove_sum
  real(kind=8)                  :: v_new, auxf_new, entries_rhs(4)
  logical                       :: flag_move=.false.
#endif

  inv2=0.5_8
  inv4=0.25_8
  inv5=0.2_8
  inv20=0.05_8
  
  tracer_rhs=0.0
  tsstiff%values(:)=0.0
 
  do elem=1, myDim_elem3d     
     elnodes=elem3D_nodes(:,elem)
     elem2 = elem2D_corresp_to_elem3D(elem)
     elem_type=grid_type_elem2d(elem2)
#ifndef use_non_hydrostatic
     elnodes2=elem2D_nodes(:,elem2)
#else
     elnodes23=nod2d_corresp_to_nod3d(elnodes)   
#endif
     fc=coriolis_param_elem2d(elem2)

     ! elementwise lateral(neutral) and vertical(dianeutral) diffusivity 
     Kh=Kh0
     if(scale_mixing_h) Kh=Kh*(voltriangle(elem2)/scalevol)**(1.0/real(scale_mixing_type))
     Kv_el=Kv0+sum(Kv(elnodes,1))/4.0

     !-------------------------------------------------------------------------------
     if (Redi_GM .and. elem_type==0) then 

        !GM diffusivity
        K_GM  = Kh*ratio_K_GM   
           
        !neutral slope
        if(nslope_version==1) then
           lay=elem3d_layer(elem)
           S = neutral_slope(:,lay,elem2)   ! S(1:3): Sx,Sy and |S|
        else
           S = neutral_slope_elem(:,elem) 
        end if

        !prepare for tapering
        !define 2 functions fcn1 and fcn2, which are required for tapering
        fcn1=1.0_8
        fcn2=1.0_8

        ! fcn1, hyperbolic tangent, used for steep slope region
        if(ODM95) fcn1 = 0.5_8*(1.0_8 + tanh((S_neutral_max - S(3))/S_d))

        !we need to check if the element is near the surface 
        !If yes, then we need the tapering function fcn2, a sine function of depth.
        if(LDD97) then
           !the first baroclinic Rossby radius
           lambd = c_speed/abs(fc)

           !limit lambda [following Large et al(1997)] to handle singularity
           !at the equator
           if (lambd < 15000.) lambd = 15000.
           if (lambd > 100000.) lambd = 100000.

           !critical depth, above which sine tapering is necessary.
           depth_scale = 100.0 + lambd*S(3)
           !here 100.0m depth is assumed to be turbulent and does not need Redi/GM
           !a well defined turbulent layer should be updated later

           !the mean depth of this element
           Z_mean = abs(sum(coord_nod3D(3,elnodes)))/4.0

           !if in the surface layer we add tapering function f2
           if (Z_mean < depth_scale)  then
              fcn2 = 0.5*(1.0 + sin(pi*Z_mean/depth_scale - pi/2.0))
           end if
        end if
      
        ! apply tapering
        ! For steep slope region:
        ! a) no taper applied to diagonal piece of horizontal neutral operator
        ! b) hyperbolic tangent(exponential) taper applied to off-diagonal piece of
        !    horizontal operator and to diagonal and off-diagonal piece of vertical
        !    neutral diffusion operator. a)+b) means we transfer the tracer diffusion
        !    to a horizontal-vertical manner in regions of steep neutral slopes.
        ! c) Exponential taper applied to GM operator.
        ! For surface layer with small slope:
        ! a) sine taper applied to both neutral operator and GM operator, except the
        !    diagonal piece of the horizontal diffusion.
        ! In one word, here we use ldd97, but always keep the diagonal part of the 
        ! horizontal diffusion following the suggestion of Griffies (2004).
        
        ! diffusion part:
        Kh_diag = Kh
        Kh_other = Kh*fcn1*fcn2
        ! skewion part:
        K_GM = K_GM*fcn1*fcn2
        
     end if  	!Redi_GM:  tapered neutral diffusivity computed
     !------------------------------------------------------------------------------

     vtri=voltriangle(elem2)
     v=voltetra(elem)
     hh=sqrt(vtri)
     hv=v/vtri     

     dx=bafux_3D(:,elem)
     dy=bafuy_3D(:,elem)
     dz=bafuz_3D(:,elem)

     uvel=uf0(elnodes)         !u*
     vvel=uf0(elnodes+myDim_nod3d+eDim_nod3D)   !v*      
     usum=sum(uvel)
     vsum=sum(vvel)

     tra_elem=tracer(elnodes,:)

     if(use_cori_semi) then
        dparam=dt_inv**2 + alpha_trapez**2*fc**2
        beta=dt_inv/dparam
        gamma=fc*alpha_trapez/dparam
     endif

     ! shortwave penetration
#ifdef use_sw_pene
     swr_conv=sum(dz*sw_3d(elnodes))*v
#endif

#ifndef use_non_hydrostatic
#ifdef use_semiimplicit_scheme
     vc=theta_ssh*ssh(elnodes2)-(gamma_stab-1.0+theta_ssh)*ssh0(elnodes2)	
#else
     vc=ssh(elnodes2)-gamma_stab*ssh0(elnodes2)
#endif
     aux1=sum(vc*bafux_2d(:,elem2))*g
     aux2=sum(vc*bafuy_2d(:,elem2))*g
     if(use_cori_semi) then
        u2d=beta*aux1+gamma*aux2
        v2d=beta*aux2-gamma*aux1
     else
        u2d=aux1*dt
        v2d=aux2*dt	
     endif
     um=inv4*usum-u2d
     vm=inv4*vsum-v2d	
     wm=sum(dz*w(elnodes))
#else
     !non-hydrostatic case
     wvel=uf0(elnodes+2*(myDim_nod3D+eDim_nod3D))          
     wsum=sum(wvel)
     vc=g*ssh(elnodes23)+nhp(elnodes) - &
          gamma_stab*g*ssh0(elnodes23)-gamma_stab_nh*nhp0(elnodes)
     aux1=sum(vc*dx)
     aux2=sum(vc*dy)
     if(use_cori_semi) then  
        u2d=beta*aux1+gamma*aux2
        v2d=beta*aux2-gamma*aux1
     else
        u2d=aux1*dt
        v2d=aux2*dt	
     endif
     w2d=sum(vc*dz)*dt
     um=inv4*usum-u2d
     vm=inv4*vsum-v2d
     wm=inv4*wsum-w2d
#endif

#ifdef use_fullfreesurf
     n_el=map_elem(elem)
     flag_move=.false.
     v_new=v
     if(n_el/=0) then
        flag_move=.true.
        v_new=voltetra_new(n_el)
        wmove=0.
        do i=1, 4
           row=elnodes(i)
           if(row<=nod2d) wmove(i)=-(ssh(row)-ssh0(row))*dt_inv
        end do
        wmove_sum=sum(wmove)  
     end if
#endif

     ept=2.0*Kh/(hh*hh)+Kv_el/(hv*hv)+abs(um)/hh+abs(vm)/hh+abs(wm)/hv + 0.01*dt_inv 
#ifdef use_fullfreesurf
     ept=ept+abs(wmove_sum)/hv 
#endif
     ept=v/ept 
     ept_dt_inv=ept*dt_inv
     auxf=v*dt_inv*inv20
#ifdef use_fullfreesurf
     auxf_new=v_new*dt_inv*inv20
#endif

     !for along sigma mixing
     if (elem_type==1) then
        lay=elem3d_layer(elem)
        S(1)=grid_slope(1,lay,elem2)
        S(2)=grid_slope(2,lay,elem2)
        aux1=S(1)**2+S(2)**2
        S(3)=sqrt(aux1)
        rotate_coe=1.0/(1.0+aux1) 
     end if


     entries=0.
     do i=1,4             ! all rows into which the element elem could contribute
        row=elnodes(i)
        if(row>myDim_nod3D) cycle                      
        do q=1,nghbr_nod3D(row)%nmb
           col_pos(nghbr_nod3D(row)%addresses(q))=q
        end do
        vnabla_t=um*dx(i)+vm*dy(i)+wm*dz(i) 

#ifndef use_fullfreesurf  
	! linear freesurface
        do q=1,4          ! all columns 
           if(elem_type==1) then !sigma diff.
              diffusion= &
                   (dx(i)*dx(q)*(1.0+S(2)*S(2))+dy(i)*dy(q)* &
                   (1.0+S(1)*S(1)))*Kh*rotate_coe + &
                   dz(i)*dz(q)*(Kv_el+Kh*S(3)*S(3)*rotate_coe) + &  
                   (S(1)*dx(i)+S(2)*dy(i))*dz(q)*Kh*rotate_coe - &  
                   (dx(i)*dy(q)+dy(i)*dx(q))*S(1)*S(2)*Kh*rotate_coe + &
                   (S(1)*dx(q)+S(2)*dy(q))*dz(i)*Kh*rotate_coe
           else
	      if (Redi_GM) then
                 diffusion=Kh_diag*(dx(i)*dx(q)+dy(i)*dy(q)) &
                      + (Kh_other-K_GM)*(S(1)*dx(i)+S(2)*dy(i))*dz(q) &
                      + (Kv_el+Kh_other*S(3)*S(3))*dz(i)*dz(q) &    
                      + (Kh_other+K_GM)*(S(1)*dx(q)+S(2)*dy(q))*dz(i) 
              else
                 diffusion=Kh*(dx(i)*dx(q)+dy(i)*dy(q))+Kv_el*dz(i)*dz(q)
	      end if
           end if
#ifndef use_non_hydrostatic
           advection=((usum+uvel(i))*inv5-u2d)*dx(q) + &
                ((vsum+vvel(i))*inv5-v2d)*dy(q) + & 
                wm*dz(q) 
#else
           advection=((usum+uvel(i))*inv5-u2d)*dx(q) + &
                ((vsum+vvel(i))*inv5-v2d)*dy(q) + & 
                ((wsum+wvel(i))*inv5-w2d)*dz(q) 
#endif  
	   entries(q)=(diffusion+advection*inv4)*v              ! diff and adv

	   advection= um*dx(q) + vm*dy(q) + wm*dz(q) 
           entries(q)=entries(q)+ ept*vnabla_t*advection       	! stabilization, upwind 

           col=elnodes(q)
           do ntr=1,num_tracer
              tracer_rhs(row,ntr)=tracer_rhs(row,ntr)-entries(q)*tra_elem(q,ntr)   ! fill tracer_rhs
           end do

           entries_t(q)=auxf                            	! mass matrix
           entries_t(q)=entries_t(q)+inv4*vnabla_t*ept_dt_inv 	! stabilization, t
        end do

        entries_t(i)=entries_t(i)+auxf                    	! completes mass matrix 
        entries=inv2*entries+entries_t
        ! put the entries to the appropriate place
        offset=tsstiff%rowptr(row)-tsstiff%rowptr(1)     
        do q=1,4        ! all columns 
           col=elnodes(q)
           ipos = offset+col_pos(col)
           tsstiff%values(ipos)=tsstiff%values(ipos)+entries(q)
        end do

#else   
	! full free surface, mass matrix updated
        do q=1,4        ! all columns 
           if(elem_type==1) then !sigma diff.
              diffusion= &
                   (dx(i)*dx(q)*(1.0+S(2)*S(2))+dy(i)*dy(q)* &
                   (1.0+S(1)*S(1)))*Kh*rotate_coe + &
                   dz(i)*dz(q)*(Kv_el+Kh*S(3)*S(3)*rotate_coe) + &  
                   (S(1)*dx(i)+S(2)*dy(i))*dz(q)*Kh*rotate_coe - &  
                   (dx(i)*dy(q)+dy(i)*dx(q))*S(1)*S(2)*Kh*rotate_coe + &
                   (S(1)*dx(q)+S(2)*dy(q))*dz(i)*Kh*rotate_coe
           else
	      if (Redi_GM) then
                 diffusion=Kh_diag*(dx(i)*dx(q)+dy(i)*dy(q)) &
                      + (Kh_other-K_GM)*(S(1)*dx(i)+S(2)*dy(i))*dz(q) &
                      + (Kv_el+Kh_other*S(3)*S(3))*dz(i)*dz(q) &    
                      + (Kh_other+K_GM)*(S(1)*dx(q)+S(2)*dy(q))*dz(i) 
	      else
                 diffusion=Kh*(dx(i)*dx(q)+dy(i)*dy(q))+Kv_el*dz(i)*dz(q)
	      end if
           end if

#ifndef use_non_hydrostatic
           advection=(dx(i)*(usum+uvel(q))+dy(i)*(vsum+vvel(q)))*inv20 + &
		dz(i)*wm*inv4 - (dx(i)*u2d+dy(i)*v2d)*inv4 
           if(flag_move) advection=advection+dz(i)*(wmove_sum+wmove(q))*inv20
#else
           advection=(dx(i)*(usum+uvel(q))+dy(i)*(vsum+vvel(q))+ &
		dz(i)*(wsum+wvel(q)))*inv20 - &
		(dx(i)*u2d+dy(i)*v2d+dz(i)*w2d)*inv4 
           if(flag_move) advection=advection+dz(i)*(wmove_sum+wmove(q))*inv20
#endif  
	   entries(q)=(diffusion-advection)*v               	! diff and adv

	   advection= um*dx(q) + vm*dy(q) + wm*dz(q)
           if(flag_move) advection=advection+wmove_sum*inv4*dz(q)
           entries(q)=entries(q)+ ept*vnabla_t*advection     	! stabilization, upwind 

           entries_t(q) = auxf_new                            
           entries_rhs(q) = entries(q)+(auxf_new-auxf)     
           if(i==q) then
              entries_t(q)=entries_t(q)+auxf_new              	! mass matrix
              entries_rhs(q)=entries_rhs(q)+(auxf_new-auxf)  	! rhs entries
           endif

           col=elnodes(q)
           do ntr=1,num_tracer
              tracer_rhs(row,ntr)=tracer_rhs(row,ntr)-entries_rhs(q)*tra_elem(q,ntr)   ! fill tracer_rhs
           end do

           entries_t(q)=entries_t(q)+inv4*vnabla_t*ept_dt_inv 	! stabilization, t
        end do

        entries=inv2*entries+entries_t
        ! put the entries to the appropriate place
        offset=tsstiff%rowptr(row)-tsstiff%rowptr(1)         
        do q=1,4        ! all columns 
           col=elnodes(q)
           ipos = offset+col_pos(col)
           tsstiff%values(ipos)=tsstiff%values(ipos)+entries(q)
        end do

#endif

        ! other rhs terms

        ! in case of considering shortwave penetration into the ocean
#ifdef use_sw_pene
        tracer_rhs(row,1)=tracer_rhs(row,1)+swr_conv*inv4 
#endif

        if(buffer_zone) then
           aux1=tracer_restore_coeff(row)*v*inv20
           do ntr=1,num_tracer
              dif=tracer0(elnodes,ntr)-tra_elem(:,ntr)
              tracer_rhs(row,ntr)=tracer_rhs(row,ntr)+aux1*(sum(dif)+dif(i))
           end do
        end if

     end do  !rows to which this element can contribute 
  end do  ! 3d element

end subroutine tsstiff_fill_gls
#endif
!
!--------------------------------------------------------------------------
!
#ifndef use_tracer_gls
subroutine tsstiff_fill_tg
  ! Fill in T/S stiffness matrix for TG scheme.
  ! mass matrix contains dt_inv
  !
  ! Coded by Sergey Danilov and Qiang Wang
  ! Reviewed by ??
  !--------------------------------------------------------
  
  use o_MATRICES
  use o_MESH
  use o_ELEMENTS
  use o_PARAM
  use o_array
  use g_config
  use g_parfe
  implicit none
  !
  integer                            :: m, i, q, row, col, ipos, offset
  integer                            :: elem, elnodes(4)
  real(kind=8)                       :: vol, inv20
  !
  inv20=0.05_8
  tsstiff%values=0.0
  allocate(ts_lump(myDim_nod3d+eDim_nod3D))      
  ts_lump=0.
  !
  !
  do elem=1,myDim_elem3d                                                 
     elnodes=elem3D_nodes(:,elem)
     vol=voltetra(elem)*dt_inv*inv20    ! contains dt_inv
     do i=1,4             ! all rows into which the element elem could contribute
        row=elnodes(i)
	if(row>myDim_nod3D) cycle               
        do q=1,nghbr_nod3D(row)%nmb
           col_pos(nghbr_nod3D(row)%addresses(q))=q
        end do
        offset=tsstiff%rowptr(row)-tsstiff%rowptr(1)   
        do q=1,4          ! all columns 
           col=elnodes(q)
           ipos = offset+col_pos(col)
	   tsstiff%values(ipos)=tsstiff%values(ipos)+vol
           if (q==i) then
              tsstiff%values(ipos)=tsstiff%values(ipos)+vol
           end if
        end do
        ts_lump(row)=ts_lump(row) + vol*5.0
     end do
  end do

  call com_3d(ts_lump)   !required for the FCT scheme
  
end subroutine tsstiff_fill_tg
#endif
!
!--------------------------------------------------------------------------
!
#ifdef use_fullfreesurf
subroutine update_matrices
  ! Update matrices due to moving surface
  !
  ! Coded by Qiang Wang
  ! Reviewed by ??
  !--------------------------------------------------------
  use o_param
  use o_array
  use o_ELEMENTS
  use o_MESH
  use o_matrices
  use g_config
  use g_parfe
  implicit none
  !
  integer          :: m, el, n_el, el2, elnodes(4), elnodes2(3), n2
  integer          :: i, q, p, row, k, col, ind, offset, ipos
  integer          :: nodup, nodlo
  real(kind=8)     :: vol, vol_new, dx2(3), dy2(3), dz(4), dz_new(4)
  real(kind=8)     :: aux, aux2, inv20, a, b, c
  real(kind=8)     :: cori_p, dparam,  beta, gamma
#ifdef use_non_hydrostatic
  real(kind=8)     :: dx(4), dy(4), dx_new(4), dy_new(4)
#endif

  inv20=0.05_8

  do el=1,myDim_elem3d                           
     n_el=map_elem(el)
     if(n_el==0) cycle 
     elnodes=elem3D_nodes(:,el)
     el2=elem2D_corresp_to_elem3D(el)
     elnodes2=elem2d_nodes(:,el2)
     vol=voltetra(el)
     vol_new=voltetra_new(n_el)
     dx2=bafux_2d(:,el2)
     dy2=bafuy_2d(:,el2)
     dz=bafuz_3d(:,el)
     dz_new=bafuz_3d_new(:,n_el)
#ifdef use_non_hydrostatic
     dx=bafux_3d(:,el)
     dy=bafuy_3d(:,el)
     dx_new=bafux_3d_new(:,n_el)
     dy_new=bafuy_3d_new(:,n_el)
#endif
     !
     !
     ! update uv matrix
     aux=(vol_new-vol)*inv20
     do i=1,4             ! all rows into which the element elem could contribute
        row=elnodes(i)
	if(row>myDim_nod3D) cycle                      
        ind=index_nod3d(row)
        if(mod(ind,10)==1) cycle
        do q=1,nghbr_nod3D(row)%nmb
           col_pos(nghbr_nod3D(row)%addresses(q))=q
        end do
        offset=uvstiff%rowptr(row)-uvstiff%rowptr(1)   
        do q=1,4          ! all columns 
           col=elnodes(q)
           ipos = offset+col_pos(col)
	   uvstiff%values(ipos)=uvstiff%values(ipos)+aux
           if (q==i) then
              uvstiff%values(ipos)=uvstiff%values(ipos)+aux
           end if
        end do
     end do
     !
     !
     if(biharmonic_visc .or. lump_uv_matrix) then
        aux=(vol_new-vol)*inv20
        do i=1,4  
           row=elnodes(i)
	   if(row>myDim_nod3D) cycle                  
           ind=index_nod3d(row)
           if(mod(ind,10)==1) cycle
           uv_lump(row)=uv_lump(row) + aux*5.0
        end do
     endif
     !
     !
     ! update ssh matrix
     if(use_cori_semi) then
        !coriolis parameter
        cori_p=coriolis_param_elem2d(el2)
        !terms for semi-implicit schemes
        dparam=dt_inv**2 + alpha_trapez**2*cori_p**2
        beta=dt_inv/dparam
        gamma=cori_p*alpha_trapez/dparam
     endif
     do i=1,3
        row=elnodes2(i)
	if(row>myDim_nod2D) cycle                         
        do q=1, nghbr_nod2D(row)%nmb
           col_pos(nghbr_nod2D(row)%addresses(q)) = q
        enddo
        do q=1,3
           col=elnodes2(q)
           if(use_cori_semi) then
              aux=(dx2(i)*dx2(q)+dy2(i)*dy2(q))*(vol_new-vol)*g*beta
              aux=aux+(dx2(i)*dy2(q)-dy2(i)*dx2(q))*(vol_new-vol)*g*gamma
           else
              aux = (dx2(i)*dx2(q) + dy2(i)*dy2(q))*(vol_new-vol)*g*dt
           endif
           ipos=sshstiff%rowptr(row)+col_pos(col)-sshstiff%rowptr(1)   
           sshstiff%values(ipos)=sshstiff%values(ipos)+aux
        end do
     end do
     !
     !
#ifdef use_non_hydrostatic
     ! update nhp matrix
     do i=1,4
        row=elnodes(i)
	if(row>myDim_nod3D) cycle                          
        do q=1, nghbr_nod3D(row)%nmb
           col_pos(nghbr_nod3D(row)%addresses(q)) = q
        enddo
        do q=1,4
           col=elnodes(q)
           if(use_cori_semi) then
              aux=(dx_new(i)*dx_new(q)+dy_new(i)*dy_new(q))*beta
              aux=aux+(dx_new(i)*dy_new(q)-dy_new(i)*dx_new(q))*gamma
              aux=(aux+dz_new(i)*dz_new(q)*dt)*vol_new
              aux2=(dx(i)*dx(q)+dy(i)*dy(q))*beta
              aux2=aux2+(dx(i)*dy(q)-dy(i)*dx(q))*gamma
              aux=aux-(aux2+dz(i)*dz(q)*dt)*vol
           else
              aux=(dx_new(i)*dx_new(q)+dy_new(i)*dy_new(q)+ &
                   dz_new(i)*dz_new(q))*vol_new*dt
              aux=aux-(dx(i)*dx(q)+dy(i)*dy(q)+dz(i)*dz(q))*vol*dt
           endif
           ipos=nhpstiff%rowptr(row)+col_pos(col)-nhpstiff%rowptr(1)  
           nhpstiff%values(ipos)=nhpstiff%values(ipos)+aux
        end do
     end do
#endif
     !
     !
#ifndef use_tracer_gls
     ! update ts matrix 
     aux=(vol_new-vol)*dt_inv*inv20
     do i=1,4             ! all rows into which the element el could contribute
        row=elnodes(i)
        if(row>myDim_nod3D) cycle                             
        do q=1,nghbr_nod3D(row)%nmb
           col_pos(nghbr_nod3D(row)%addresses(q))=q
        end do
        offset=tsstiff%rowptr(row)-tsstiff%rowptr(1)            
        do q=1,4          ! all columns 
           col=elnodes(q)
           ipos = offset+col_pos(col)
           tsstiff%values(ipos)=tsstiff%values(ipos)+aux
           if (q==i) then
              tsstiff%values(ipos)=tsstiff%values(ipos)+aux
           end if
        end do
        ts_lump(row)=ts_lump(row) + aux*5.0
     end do
#endif
  end do
  !
  !
#ifndef use_non_hydrostatic
  ! update w potential matrix
  
  do n2=1,myDim_nod2d
     do k=1,2   !only the first layer is moved
        a=0.
        b=0.
        c=0.
        row=nod3d_below_nod2d(k,n2)
        if(k==2) nodup=nod3d_below_nod2d(k-1,n2)
        if(k==1) nodlo=nod3d_below_nod2d(k+1,n2)

        do i=1,nod_in_elem3D(row)%nmb
           el=nod_in_elem3D(row)%addresses(i)
           elnodes=elem3D_nodes(:,el)  
           n_el=map_elem(el)
           if(n_el==0) cycle
           vol=voltetra(el)
           vol_new=voltetra_new(n_el)
           dz=bafuz_3d(:,el)
           dz_new=bafuz_3d_new(:,n_el)

           do q=1,4
              if(elnodes(q)==row) then
                 p=q
                 exit
              end if
           end do         

           !first entry
           if(k>1) then
              do q=1,4
                 if(elnodes(q)==nodup) then
                    a=a + dz_new(p)*dz_new(q)*vol_new - dz(p)*dz(q)*vol
                    exit
                 end if
              end do
           end if
           
           !second entry
           b=b + dz_new(p)*dz_new(p)*vol_new - dz(p)*dz(p)*vol
        
           !third entry
           if(k<2) then
              do q=1,4
                 if(elnodes(q)==nodlo) then
                    c=c + dz_new(p)*dz_new(q)*vol_new - dz(p)*dz(q)*vol
                    exit
                 end if
              end do
           end if
        end do  !i

        wpot_matrix(1,k,n2)= wpot_matrix(1,k,n2)+a
        wpot_matrix(2,k,n2)= wpot_matrix(2,k,n2)+b
        wpot_matrix(3,k,n2)= wpot_matrix(3,k,n2)+c

     end do  !k
  end do  !n2
#endif

end subroutine update_matrices
#endif
!
!------------------------------------------------------------------
!
