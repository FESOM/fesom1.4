!Performs a time step of the ice model

subroutine ice_step
  use o_mesh
  use o_param
  use i_therm_parms
  use i_dyn_parms
  use i_array
  use i_solver
  use g_parfe
  use g_config
  implicit none
  !
  integer        :: m, row
  real(kind=8)   :: t0, t1, t2, t3

  t0=MPI_Wtime()

  !1) Dynamic part

  ! to avoid zero ice velocities (to avoid the solver's sudden death)
  do row=1,myDim_nod2d+eDim_nod2d    
     if (m_ice(row) < 0.001) then
        u_ice(row)=1e-10   !! To reach correspondence with 
        v_ice(row)=1e-10   !! global memory code
     endif
  enddo

  ! to solve u,v
  call rheology      

#ifdef use_cavity
  call clean_cavity_vel
#endif    
  t1=MPI_Wtime()

  ! to solve m, a, ms due to advection
#ifdef use_ice_gls
  call icestiff_fill_gls
  call solveIce(solve_m_ice)       
  call solveIce(solve_a_ice)       
  call solveIce(solve_m_snow)      
  m_ice=m_ice+dm_ice
  a_ice=a_ice+da_ice
  m_snow=m_snow+dm_snow 
  call com_2D(m_ice)
  call com_2D(a_ice)
  call com_2D(m_snow)  
#else

#ifdef use_ice_fct
  call ice_rhs_tg

  call fct_ice_solve

#else
  call ice_rhs_tg                 
  if(lump_ice_matrix) then
     call ice_solve
  else                            
     call solveIce(solve_m_ice)       
     call solveIce(solve_a_ice)       
     call solveIce(solve_m_snow)      
  end if
  m_ice=m_ice+dm_ice
  a_ice=a_ice+da_ice
  m_snow=m_snow+dm_snow
  call com_2D(m_ice)
  call com_2D(a_ice)
  call com_2D(m_snow)  

#endif
#endif

  call cut_off

#ifdef use_cavity
  call clean_cavity_m
#endif
  t2=MPI_Wtime()   

  !2) Thermodynamic part (update m, a, ms)

  call thermodynamics
  t3=MPI_Wtime()

#ifdef VERBOSE
  if (mod(istep,logfile_outfreq)==0 .and. mype==0) then
     write(*,*) 
     write(*,*) 'ice took      ', t3-t0
     write(*,*) 'ice uv        ', t1-t0
     write(*,*) 'ice mass      ', t2-t1
     write(*,*) 'ice thermo    ', t3-t2
  endif
#endif
end subroutine ice_step
!
!------------------------------------------------------------------------------
!
subroutine clean_cavity_vel
  use o_mesh
  use i_array
  use g_parfe
  implicit none
  integer        :: m, row

  do row=1,myDim_nod2d+eDim_nod2d           
     if(cavity_flag_nod2d(row)==1) then  
        u_ice(row)=0.
        v_ice(row)=0.
     endif
  enddo
end subroutine clean_cavity_vel
!
!------------------------------------------------------------------------------
!
subroutine clean_cavity_m
  use o_mesh
  use i_array
  use g_parfe
  implicit none
  integer        :: m, row

  do row=1,myDim_nod2d+eDim_nod2d            
     if(cavity_flag_nod2d(row)==1) then 
        m_ice(row)=0.
        a_ice(row)=0.
        m_snow(row)=0.
     endif
  enddo
end subroutine clean_cavity_m
!
!------------------------------------------------------------------------------
!
