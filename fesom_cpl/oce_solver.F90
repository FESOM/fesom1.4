! Solver for the ocean model


#ifdef PETSC

subroutine solve(ident)

  ! All stiffness matrices are in the transposed form ?
  ! PIL_PMVALS --- each processor takes only owned rows
  ! PIL_MVALS  --- matrix should be present at each processor
  ! RHS  --- only owned rows are taken as default

  use g_PARFE
  use g_config
  use o_PARAM
  use o_array
  use o_MATRICES
  use o_MESH
  use o_solver
  implicit none

#include "petscf.h"
  !include "petscf.h"
  !include "pilutf.h"
  !include "hypref.h"

  integer,intent(IN)    :: ident
  integer               :: Pmode, n3 
  integer               :: maxiter, restart, lutype, fillin
  integer               :: COMM, myrows
  real(kind=8)          :: droptol, soltol, rinfo(20,20)
  save rinfo

  COMM=MPI_COMM_FESOM

  maxiter=2000
  restart=15
  fillin=3
  lutype=2
  droptol=1.e-6
  soltol=1.e-6
  n3=myDim_nod3D+eDim_nod3D
  call MPI_Barrier(MPI_COMM_FESOM, MPIERR)

  if (mype==0) WRITE (*,*) 'Solver: PETSc'

#ifdef use_fullfreesurf     
  ! matrix is updated

  if (ident==solve_u) then  ! du* or du
     Pmode = PET_BLOCKP+PET_STRUCT + PET_PMVALS + PET_SOLVE + PET_BICGSTAB + PET_RCM &
          + PET_QUIET 
     if (.not. iteruv_first)   Pmode = Pmode - PET_STRUCT  + PET_NEWPC 
     call PETSC_S(Pmode, &
          ident, uvstiff%dim, uvstiff%nza,  myrows, &
          maxiter, & 
          restart, &  
          fillin,  &  
          droptol, & 
          soltol,  &   
          part3D, uvstiff%rowptr, uvstiff%colind, uvstiff%values, &
          uv_rhs(1:n3), duf(1:n3), rinfo(:,1), COMM) 
  endif

  if (ident==solve_v) then  ! dv* or dv
     Pmode = PET_BLOCKP+PET_SOLVE + PET_BICGSTAB + PET_QUIET 
     call PETSC_S(Pmode, &
          solve_u, uvstiff%dim, uvstiff%nza,  myrows, &
          maxiter, & 
          restart, &  
          fillin,  &  
          droptol, & 
          soltol,  &   
          part3D, uvstiff%rowptr, uvstiff%colind, uvstiff%values, &
          uv_rhs(1+n3:n3*2), duf(1+n3:n3*2), rinfo(:,2), COMM) 
  endif

#ifdef use_non_hydrostatic
  if (ident==solve_w) then  ! dw* or dw for non-hydrostatic case
     Pmode = PET_BLOCKP+PET_SOLVE + PET_BICGSTAB + PET_QUIET
     call PETSC_S(Pmode, &
          solve_u, uvstiff%dim, uvstiff%nza,  myrows, &
          maxiter, & 
          restart, &  
          fillin,  &  
          droptol, & 
          soltol,  &   
          part3D, uvstiff%rowptr, uvstiff%colind, uvstiff%values, &
          uv_rhs(1+n3*2:n3*3), duf(1+n3*2:n3*3), rinfo(:,3), COMM) 
  endif
#endif

  if (ident==solve_ssh) then  !dssh
     Pmode = PET_BLOCKP+PET_STRUCT + PET_PMVALS + PET_SOLVE + PET_BICGSTAB + PET_ILU + PET_PCBJ + PET_QUIET
     !if (.not. iter_first)   Pmode = Pmode - PET_STRUCT
     if (.not. iter_first)   Pmode = Pmode - PET_STRUCT  + PET_NEWPC 
     call PETSC_S(Pmode,ident,sshstiff%dim, sshstiff%nza, myrows, &
          maxiter, & 
          restart, &  
          fillin,  &  
          droptol*1.0e-2, &  
          soltol*1.0e-2,  &   ! Sol. Residuum reduction 
          part2D, sshstiff%rowptr, sshstiff%colind, sshstiff%values, &
          ssh_rhs, dssh, rinfo(:,4), COMM) 
  endif

#ifdef use_non_hydrostatic  
  if (ident==solve_nhp) then
     Pmode = PET_BLOCKP+PET_STRUCT + PET_PMVALS + PET_SOLVE + PET_BICGSTAB + PET_ILU + PET_PCBJ +PET_QUIET 
     !if (.not. iter_first)   Pmode = Pmode - PET_STRUCT  + PET_NEWPC
     if (.not. iter_first)   Pmode = Pmode - PET_STRUCT  
     call PETSC_S(Pmode,ident,nhpstiff%dim, nhpstiff%nza, myrows, &
          maxiter, & 
          restart, &  
          fillin,  &  
          droptol*10., &  
          soltol*10.,  &  
          part3D, nhpstiff%rowptr, nhpstiff%colind, nhpstiff%values, &
          nhp_rhs, nhp, rinfo(:,5), COMM) 
  endif
#endif

  if(ident>=solve_tra) then
     if(ident==solve_tra) then  ! dTF
        Pmode = PET_BLOCKP+PET_STRUCT + PET_PMVALS + PET_SOLVE + PET_BICGSTAB + PET_RCM &
             + PET_QUIET
        if (.not. iter_first)  Pmode = Pmode - PET_STRUCT  + PET_NEWPC 
        call PETSC_S(Pmode, &
             ident,tsstiff%dim, tsstiff%nza,  myrows, &
             maxiter, & 
             restart, &  
             fillin,  &  
             droptol, &  ! droptol 
             soltol, &   ! Sol. Residuum reduction 
             part3D, tsstiff%rowptr, tsstiff%colind, tsstiff%values, &
             tracer_rhs(:,1), dtracer(:,1), rinfo(:,6), COMM) 
     else  !other tracers
        Pmode = PET_BLOCKP+PET_SOLVE + PET_BICGSTAB + PET_QUIET  
        call PETSC_S(Pmode , &
             solve_tra, &
             tsstiff%dim, tsstiff%nza,  myrows, &
             maxiter, & 
             restart, &  
             fillin,  &  
             droptol, &  
             soltol, &   ! Sol. Residuum reduction 
             part3D, tsstiff%rowptr, tsstiff%colind, tsstiff%values, &
             tracer_rhs(:,ident-solve_tra+1), dtracer(:,ident-solve_tra+1), rinfo(:,7), COMM)
     endif
  endif


  !-------------------------------------------------------------------------------------
#else   

  ! linear surface

  if (ident==solve_u) then  ! du* or du
     Pmode =PET_BLOCKP+PET_SOLVE + PET_BICGSTAB + PET_QUIET
     if (iteruv_first) Pmode=Pmode+PET_STRUCT+PET_PMVALS + PET_RCM
     call PETSC_S(Pmode, &
          ident, uvstiff%dim, uvstiff%nza,  myrows, &
          maxiter, & 
          restart, &  
          fillin,  &  
          droptol, & 
          soltol,  &   
          part3D, uvstiff%rowptr, uvstiff%colind, uvstiff%values, &
          uv_rhs(1:n3), duf(1:n3), rinfo(:,1), COMM) 
  endif

  if (ident==solve_v) then  ! dv* or dv
     Pmode = PET_BLOCKP+PET_SOLVE + PET_BICGSTAB + PET_QUIET 
     call PETSC_S(Pmode, &
          solve_u, uvstiff%dim, uvstiff%nza,  myrows, &
          maxiter, & 
          restart, &  
          fillin,  &  
          droptol, & 
          soltol,  &   
          part3D, uvstiff%rowptr, uvstiff%colind, uvstiff%values, &
          uv_rhs(1+n3:n3*2), duf(1+n3:n3*2), rinfo(:,2), COMM) 
  endif

#ifdef use_non_hydrostatic
  if (ident==solve_w) then  ! dw* or dw for non-hydrostatic case
     Pmode = PET_BLOCKP+PET_SOLVE + PET_BICGSTAB + PET_QUIET 
     call PETSC_S(Pmode, &
          solve_u, uvstiff%dim, uvstiff%nza,  myrows, &
          maxiter, & 
          restart, &  
          fillin,  &  
          droptol, & 
          soltol,  &   
          part3D, uvstiff%rowptr, uvstiff%colind, uvstiff%values, &
          uv_rhs(1+n3*2:n3*3), duf(1+n3*2:n3*3), rinfo(:,3), COMM) 
  endif
#endif

  if (ident==solve_ssh) then  ! dssh
     Pmode =PET_BLOCKP+PET_SOLVE+ PET_BICGSTAB+PET_QUIET !+PET_REPORT
     if (iter_first) Pmode=Pmode+PET_STRUCT+PET_PMVALS + PET_ILU + PET_PCBJ 
     call PETSC_S(Pmode,ident,sshstiff%dim, sshstiff%nza, myrows, &
          maxiter, & 
          restart, &  
          fillin,  &  
          droptol*1.0e-2, &  
          soltol*1.0e-2,  &   ! Sol. Residuum reduction 
          part2D, sshstiff%rowptr, sshstiff%colind, sshstiff%values, &
          ssh_rhs, dssh, rinfo(:,4), COMM) 
  endif

#ifdef use_non_hydrostatic  
  if (ident==solve_nhp) then
     Pmode =PET_BLOCKP+PET_SOLVE+ PET_BICGSTAB !+PET_QUIET 
     if (iter_first) Pmode=Pmode+PET_STRUCT+PET_PMVALS+PET_PCBJ+PET_ILU
     call PETSC_S(Pmode,ident,nhpstiff%dim, nhpstiff%nza, myrows, &
          maxiter, & 
          restart, &  
          fillin,  &  
          droptol*10., &  
          soltol*10.,  &  
          part3D, nhpstiff%rowptr, nhpstiff%colind, nhpstiff%values, &
          nhp_rhs, nhp, rinfo(:,5), COMM) 
  endif
#endif

  if(ident>=solve_tra) then
     if(ident==solve_tra) then  ! dTF
#ifdef use_tracer_gls
        Pmode = PET_BLOCKP+PET_STRUCT + PET_PMVALS + PET_SOLVE + PET_BICGSTAB + PET_RCM + PET_QUIET
        if (.not. iter_first)  Pmode = Pmode - PET_STRUCT  + PET_NEWPC
#else
        Pmode =PET_BLOCKP+PET_SOLVE + PET_BICGSTAB + PET_QUIET
        if (iter_first) Pmode=Pmode+PET_STRUCT+PET_PMVALS + PET_RCM
#endif
        call PETSC_S(Pmode, &
             ident,tsstiff%dim, tsstiff%nza,  myrows, &
             maxiter, & 
             restart, &  
             fillin,  &  
             droptol, &  
             soltol, &   
             part3D, tsstiff%rowptr, tsstiff%colind, tsstiff%values, &
             tracer_rhs(:,1), dtracer(:,1), rinfo(:,2), COMM) 

     else  !other tracers

        Pmode =PET_BLOCKP+ PET_SOLVE + PET_BICGSTAB + PET_QUIET  
        call PETSC_S(Pmode , &
             solve_tra, &
             tsstiff%dim, tsstiff%nza,  myrows, &
             maxiter, & 
             restart, &  
             fillin,  &  
             droptol, &  
             soltol, &  
             part3D, tsstiff%rowptr, tsstiff%colind, tsstiff%values, &
             tracer_rhs(:,ident-solve_tra+1), dtracer(:,ident-solve_tra+1), rinfo(:,2), COMM)
     endif
  endif
#endif

end subroutine solve	

#elif defined(PARMS)

subroutine solver_init(ident, stiff, stype, pctype, pcilutype, &
            lutype, fillin, droptol, maxiter, restart, soltol, reuse, comm)

  use g_PARFE
  use g_config
  use o_PARAM
  use o_array
  use o_MATRICES
  use o_MESH
  use o_solver
  use o_data_types

  implicit none
  
  integer, intent(in) :: ident, stype, pctype, pcilutype, lutype
  integer, intent(in) :: fillin, maxiter, restart, reuse, comm
  real(kind=8), intent(in) :: droptol, soltol
  type(sparse_matrix), intent(in) :: stiff

  if (mype==0) WRITE (*,*) 'Solver: pARMS'

if(stiff%dim == nod2D) then
  call psolver_init(ident, stype, pctype, pcilutype, lutype, &
              fillin, droptol, maxiter, restart, soltol, &
              part2D-1, stiff%rowptr(:)-stiff%rowptr(1), &
	      stiff%colind-1, stiff%values, reuse, comm)
else
  call psolver_init(ident, stype, pctype, pcilutype, lutype, &
              fillin, droptol, maxiter, restart, soltol, &
              part3D-1, stiff%rowptr(:)-stiff%rowptr(1), &
	      stiff%colind-1, stiff%values, reuse, comm)
end if
end subroutine solver_init

subroutine solve(ident)

  ! All stiffness matrices are in the transposed form ?
  ! PIL_PMVALS --- each processor takes only owned rows
  ! PIL_MVALS  --- matrix should be present at each processor
  ! RHS  --- only owned rows are taken as default

  use g_PARFE
  use g_config
  use o_PARAM
  use o_array
  use o_MATRICES
  use o_MESH
  use o_solver
 
  implicit none
  
#include "fparms.h"

  integer,intent(IN)    :: ident
  integer               :: n3, comm
  integer               :: maxiter, restart, lutype, fillin
  real(kind=8)          :: droptol, soltol
  COMM=MPI_COMM_FESOM

  maxiter=2000
  restart=15
  fillin=3
  lutype=2
  droptol=1.e-6
  soltol=1.e-6

  n3=myDim_nod3D+eDim_nod3D
  call MPI_Barrier(COMM, MPIERR)

#ifdef use_fullfreesurf     
  ! matrix is updated
 
  if(ident == solve_u) then
    if(iteruv_first) then
       call solver_init(ident, uvstiff, SOLBICGS, PCBJ, PCILUK, &
            lutype, fillin, droptol, maxiter, restart, soltol, 1, comm)
      call psolve(ident, uv_rhs(1:n3), 0., duf(1:n3), 0)			 
    else
      call psolve(ident, uv_rhs(1:n3), uvstiff%values, duf(1:n3), 1)			 
    end if	       

 end if
  if(ident == solve_v) then
     call psolve(solve_u, uv_rhs(1+n3:2*n3), 0., duf(1+n3:2*n3),0)			 
  end if

#ifdef use_non_hydrostatic
  if (ident==solve_w) then  ! dw* or dw for non-hydrostatic case
     call psolve(solve_u, uv_rhs(1+2*n3:3*n3), 0., duf(1+2*n3:3*n3),0) 
  end if	 
#endif

  if (ident==solve_ssh) then  !dssh
    if(iter_first) then
      call solver_init(ident, sshstiff, SOLBICGS, PCBJ, PCILUK, &
           lutype, fillin, droptol*1.e-2, maxiter, restart, soltol*1.e-2, 1, comm)
      call psolve(ident, ssh_rhs, 0., dssh, 0)
    else    
      call psolve(ident, ssh_rhs, sshstiff%values, dssh, 1) 
    end if	 
  endif

#ifdef use_non_hydrostatic  
  if (ident==solve_nhp) then
    if(iter_first) then
      call solver_init(ident, nhpstiff, SOLBICGS, PCBJ, PCILUK, &
           lutype, fillin, droptol*10., maxiter, restart, soltol*10., 1, comm)
      call psolve(ident, nhp_rhs, 0., nhp, 0)	
    else
      call psolve(ident, nhp_rhs, nhpstiff%values, nhp, 1)	
    end if  
  end if   
#endif

  if(ident>=solve_tra) then
    if(ident == solve_tra) then
      if(iter_first) then
        call solver_init(ident, tsstiff, SOLBICGS, PCBJ, PCILUK, &
             lutype, fillin, droptol, maxiter, restart, soltol, 1, comm)
	     call psolve(ident, tracer_rhs(:,1), 0., dtracer(:,1), 0)
      else	
	call psolve(ident, tracer_rhs(:,1), tsstiff%values, dtracer(:,1), 1)
      end if		   
   else
      call psolve(solve_tra, tracer_rhs(:,ident-solve_tra+1),&
                  0., dtracer(:,ident-solve_tra+1, 0)
    end if	
  endif


  !-------------------------------------------------------------------------------------
#else   
  ! linear surface

  if (ident==solve_u) then  ! du* or du

    if(iteruv_first) then
      call solver_init(ident, uvstiff, SOLBICGS, PCBJ, PCILUK, &
      lutype, fillin, droptol, maxiter, restart, soltol, 0, comm)
    end if 
    call psolve(ident, uv_rhs(1:n3), 0., duf(1:n3), 0)			        
  end if
  if(ident == solve_v) then
    call psolve(solve_u, uv_rhs(1+n3:2*n3), 0., duf(1+n3:2*n3), 0)			 
  end if

#ifdef use_non_hydrostatic
  if (ident==solve_w) then  ! dw* or dw for non-hydrostatic case
    call psolve(solve_u, uv_rhs(1+2*n3:3*n3), 0., duf(1+2*n3:3*n3), 0) 
  end if
#endif
  if (ident==solve_ssh) then  ! dssh
    if(iter_first) then
      call solver_init(ident, sshstiff, SOLBICGS, PCBJ, PCILUK, &
		       lutype, fillin, droptol*1.e-2, maxiter, restart, &
		         soltol*1.e-2, 0, comm)
    end if
    call psolve(ident, ssh_rhs, 0., dssh, 0)
  endif
#ifdef use_non_hydrostatic  
  if (ident==solve_nhp) then
    if(iter_first) then
      call solver_init(ident, nhpstiff, SOLBICGS, PCBJ, PCILUK, &
                       lutype, fillin, droptol*10., maxiter, restart,
  soltol*10., 0, comm)
    end if	
    call psolve(ident, nhp_rhs, 0., nhp, 0)	
  end if   
#endif

  if(ident>=solve_tra) then
    if(ident==solve_tra) then  ! dTF


#ifdef use_tracer_gls

        if(iter_first) then
          call solver_init(ident, tsstiff, SOLBICGS, PCBJ, PCILUK, &
               lutype, fillin, droptol, maxiter, restart, soltol, 1, comm)
	  call psolve(ident, tracer_rhs(:,1), 0., dtracer(:,1),0)	
        else
          call psolve(ident, tracer_rhs(:,1), tsstiff%values, dtracer(:,1),1)
	end if	
#else 
      if(iter_first) then
         call solver_init(ident, tsstiff, SOLBICGS, PCBJ, PCILUK, &
               lutype, fillin, droptol, maxiter, restart, soltol, 0, comm)
      end if 
      call psolve(ident, tracer_rhs(:,1), 0., dtracer(:,1),0)	
       
#endif
    else 
     call psolve(solve_tra, tracer_rhs(:,ident-solve_tra+1),&
                  0., dtracer(:,ident-solve_tra+1), 0)
    end if	
  endif

#endif

end subroutine solve

#endif
