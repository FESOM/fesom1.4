! solver for solving ice advection equations

subroutine solveIce(ident)
  use g_PARFE
  use o_mesh
  use o_solver
  USE i_array
  use i_solver
  implicit none

#ifdef PETSC
#include "petscf.h"
  !include "pilutf.h"
  !include "hypref.h"
#endif

#ifdef PARMS
#include "fparms.h"
#endif

  integer        :: ident
  INTEGER        :: Pmode, restart, maxiter, lutype, fillin
  INTEGER        :: COMM, myrows
  REAL(kind=8)   :: droptol, soltol, rinfo(20,20)
  save rinfo

  COMM=MPI_COMM_FESOM

  maxiter=2000
  restart=15
  fillin=2
  lutype=2
  droptol=1.e-6
  soltol=1.e-6

  call MPI_Barrier(COMM, MPIERR)

#ifdef PETSC

  if (ident==solve_m_ice) then ! m_ice
#ifdef use_ice_gls
     Pmode=PET_BLOCKP+PET_SOLVE+PET_BICGSTAB+PET_PMVALS+PET_RCM + PET_QUIET !+PET_REPORT
     if (iter_first) Pmode=Pmode+PET_STRUCT 
#else
     Pmode=PET_BLOCKP+PET_SOLVE+PET_BICGSTAB + PET_QUIET !+PET_REPORT
     if (iter_first) Pmode=Pmode+PET_STRUCT+PET_PMVALS+PET_RCM
#endif
     call PETSC_S(Pmode,ident,icestiff%dim, icestiff%nza,  myrows, &
          maxiter, & 
          restart, &  
          fillin,  &  
          droptol, &  
          soltol,  &  
          part2D, icestiff%rowptr, icestiff%colind, icestiff%values, &
          rhs_m, dm_ice, rinfo(:,12), COMM)
  endif

  if (ident==solve_a_ice) then ! a_ice
     Pmode=PET_BLOCKP+PET_SOLVE+PET_BICGSTAB +PET_QUIET !+PET_REPORT
     call PETSC_S(Pmode, &
          ident-1, &
          icestiff%dim, icestiff%nza, myrows, &
          maxiter, & 
          restart, &  
          fillin,  &  
          droptol, &  
          soltol,  &   
          part2D, icestiff%rowptr, icestiff%colind, icestiff%values, &
          rhs_a, da_ice,  rinfo(:,13), COMM)
  endif

  if (ident==solve_m_snow) then ! m_snow
     Pmode=PET_BLOCKP+PET_SOLVE+PET_BICGSTAB +PET_QUIET !+PET_REPORT
     call PETSC_S(Pmode, &
          ident-2, &
          icestiff%dim, icestiff%nza, myrows, &
          maxiter, & 
          restart, &  
          fillin,  &  
          droptol, &  
          soltol,  &  
          part2D, icestiff%rowptr, icestiff%colind, icestiff%values, &
          rhs_ms, dm_snow,  rinfo(:,14), COMM)
  endif

#elif defined(PARMS)
write(*,*) 'START: calling ice solver'      
  if (ident==solve_m_ice) then ! m_ice
    if(iter_first) then	
      call solver_init(0, icestiff, SOLBICGS, PCBJ, PCILUK, &
                       lutype, fillin, droptol, maxiter, restart, soltol, 1, COMM)
      call psolve(0, rhs_m, 0., dm_ice,0)
    else
      call psolve(0, rhs_m, icestiff%values, dm_ice,1)
    end if	 
  end if	 
  if (ident==solve_a_ice) then ! a_ice
    call psolve(0, rhs_a, 0., da_ice, 0)
  endif

  if (ident==solve_m_snow) then ! m_snow
    call psolve(0, rhs_ms, 0., dm_snow, 0)  		
 end if	 
write(*,*) 'END: calling ice solver'      
#endif

end subroutine solveIce

! ==================================================================
