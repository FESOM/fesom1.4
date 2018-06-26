!=============================================================================
!  Performs a time step of the ocean model
!=============================================================================

subroutine ocean_step
  use o_param
  use o_array
  use o_mixing_kpp_mod
  use o_mixing_pp_mod
  use o_mixing_my2p5_mod
  use o_mixing_tidal_mod
  use o_passive_tracer_mod
  use o_age_tracer_mod
  use o_mesh
  use o_elements
  use o_solver
  use g_config
  use g_parfe
  use cmor_variables_diag, only : compute_diag

  implicit none
  integer      :: i, m, row, row2, row3, n3
  real(kind=8) :: t0,t1,t2,t3,t4,t5,t6,t7, t8,t9,t10, t11

  t0=MPI_Wtime() 

  n3=ToDim_nod3d

  if (use_global_tides) then
     call potential
  end if

  do row=1,n3                           
     row2=row+n3                   
     uf0(row)=uf(row)                 ! uf0 & uf: u^n
     uf0(row2)=uf(row2) 
#ifdef use_non_hydrostatic
     row3=row2+n3                     
     uf0(row3)=uf(row3)
#endif
  end do

  !-----------------------------------------------------------

  call compute_density
  call compute_bfsq
  call compute_dbsfc


  if(grid_type/=2) then
     call compute_pressure            ! compute hpressure
  end if
  if(grid_type/=1) then
     call compute_pressure_force      ! compute pressure gradient
  end if
  t2=MPI_Wtime()    

  !-----------------------------------------------------------

  if(trim(mix_scheme)=='KPP') then
     call oce_mixing_kpp(Av, Kv)
     call convect_adjust
  elseif(trim(mix_scheme)=='PP') then
     call oce_mixing_pp
  elseif(trim(mix_scheme)=='MY2p5') then
     call oce_mixing_MY2p5
  end if
  if(tidal_mixing) call oce_mixing_tidal(Av, Kv)
  t3=MPI_Wtime()

  !------------------------------------------------------------

  call velocity_rhs                 ! rhs for u*
  if(use_vertvisc_impl) call uv_sfc_bott_bc
  t4=MPI_Wtime()    

  if(lump_uv_matrix) then           ! solver for du*
     call uv_solve
  else
     call solve(solve_u)   
     iteruv_first=.false.
     call solve(solve_v)
     call com_3d(duf(1:n3))              
     call com_3d(duf(1+n3:2*n3))           
  endif

#ifdef use_non_hydrostatic 
  call solve(solve_w)
  call com_3d(duf(1+2*n3:3*n3))          
#endif 
  do row=1,n3                            
     row2=row+n3                   
     uf(row)=uf(row)+duf(row)      ! uf: u*
     uf(row2)=uf(row2)+duf(row2)
#ifdef use_non_hydrostatic
     row3=row2+n3                         
     uf(row3)=uf(row3)+duf(row3)
#endif
  end do

  if(use_vertvisc_impl) then      ! apply implicit vertical viscosity
     call impl_vertvisc
  end if
  t5=MPI_Wtime() 

  !--------------------------------------------------------------

  ssh0=ssh                        ! ssh & ssh0: ssh^n  

#ifdef use_opbnd_tide
  call update_tidal_opbnd         ! update tidal open boundary
#endif

  call compute_ssh_rhs            ! ssh rhs

  call solve(solve_ssh)           ! solve dssh
  ssh=ssh0+dssh                   ! ssh: ssh^n+1
  call com_2D(ssh)                       

#ifdef use_non_hydrostatic
  nhp0=nhp
  call compute_nhp_rhs            ! nhp rhs
  call solve(solve_nhp)           ! solve nhp
  call com_3D(nhp)                       
#endif   
  t6=MPI_Wtime()      

  do row=1,n3                              
     row2=row+n3                        
     uf0(row)=uf(row)             ! uf0: u*     
     uf0(row2)=uf(row2) 
#ifdef use_non_hydrostatic
     row3=row2+n3                       
     uf0(row3)=uf(row3)
#endif	
  end do

  !---------------------------------------------------------------

  call velocity_rhs_update        ! Update rhs: contribution from ssh/nhp
  t7=MPI_Wtime()

  if(lump_uv_matrix) then         ! solve for full du 
     call uv_solve
  else
     call solve(solve_u)                          
     call solve(solve_v)
     call com_3d(duf(1:n3))             
     call com_3d(duf(1+n3:2*n3))        
  endif
#ifdef use_non_hydrostatic
  call solve(solve_w)
  call com_3d(duf(1+2*n3:3*n3))         
#endif
  do row=1,n3                            
     row2=row+n3                         
     uf(row)=uf(row)+duf(row)     ! uf: u^n+1  
     uf(row2)=uf(row2)+duf(row2)
#ifdef use_non_hydrostatic
     row3=row2+n3                       
     uf(row3)=uf(row3)+duf(row3)
#endif	
  end do
  t8=MPI_Wtime()

  !---------------------------------------------------------

#ifndef use_non_hydrostatic
  call compute_vvel_rhs           ! vertical rhs
  call solve_wpot                 ! solve for w potential 
#endif
  t9=MPI_Wtime()  

  !----------------------------------------------------------

#ifdef use_fullfreesurf
  call update_mesh
  call update_matrices
#endif 
  t10=MPI_Wtime() 

  !-----------------------------------------------------------

#ifdef use_cavity
  call cavity_heat_water_fluxes
#endif

  !-----------------------------------------------------------

  if(Redi_GM) call prepare_neutral_physis 

  !-----------------------------------------------------------
  
  if(brine_rejection_param) call cal_brine_rejection

#ifdef use_tracer_gls
  call tsstiff_fill_gls          ! tracer matrix/rhs
  call ts_sfc_bc
  if(use_passive_tracer) call ptr_sfc_bc
  if(use_age_tracer) call age_tracer_tendency
  t11=MPI_Wtime()
  do i=1,num_tracer
     call solve(solve_tra+i-1)   ! solve for tracer
     call com_3D(dtracer(:,i))          
  end do
  do i=1,num_tracer                     
     tracer(:,i)=tracer(:,i)+dtracer(:,i)
  end do

#else

#ifdef use_tracer_fct
  call tracer_rhs_tg
  call ts_sfc_bc
  if(use_passive_tracer) call ptr_sfc_bc
  if(use_age_tracer) call age_tracer_tendency
  t11=MPI_Wtime()
  call fct_tracer_solve 
#else
  call tracer_rhs_tg
  call ts_sfc_bc
  if(use_passive_tracer) call ptr_sfc_bc
  if(use_age_tracer) call age_tracer_tendency
  t11=MPI_Wtime()
  if(lump_ts_matrix) then
     call tracer_solve
  else
     do i=1,num_tracer
        call solve(solve_tra+i-1)
	call com_3D(dtracer(:,i))
     end do
  endif
  do i=1,num_tracer
     do row=1,n3                         
        tracer(row,i)=tracer(row,i)+dtracer(row,i)
     end do
  end do
#endif

  if(brine_rejection_param) call apply_brine_rejection

  if(use_vertdiff_impl) then
     call impl_vertdiff        ! apply implicit vertical diff.
  end if
#endif

  if(use_passive_tracer) call ptr_cutoff_restore
  if(use_age_tracer) call age_tracer_cutoff_restore

  !--------------------------------------------------------------

  t1=MPI_Wtime()
  iter_first = .false.
#ifdef VERBOSE
  if (mod(istep,logfile_outfreq)==0 .and. mype==0) then      
     write(*,*)
     write(*,*) 'ocean  took   ', t1-t0
     write(*,*) 'Dens/pressure ', t2-t0
     write(*,*) 'mixing scheme ', t3-t2
     write(*,*) 'v rhs         ', t4-t3
     write(*,*) 'solve v_star  ', t5-t4
     write(*,*) 'pressure_solve', t6-t5
     write(*,*) 'rhs_update    ', t7-t6
     write(*,*) 'solve full v  ', t8-t7
#ifndef use_non_hydrostatic
     write(*,*) 'vert_vel_solve', t9-t8
#endif
#ifdef use_fullfreesurf
     write(*,*) 'update mesh   ', t10-t9
#endif
     write(*,*) 'tra_assemble  ', t11-t10
     write(*,*) 'tra_solve     ', t1-t11
     write(*,*)
  endif
#endif
!write(*,*) 'mype, min(u),   max(u)  =', mype, minval(uf(1:n3)),      maxval(uf(1:n3))
!write(*,*) 'mype, min(v),   max(v)  =', mype, minval(uf(n3+1:2*n3)), maxval(uf(n3+1:2*n3))
!write(*,*) 'mype, min(ssh), max(ssh)=', mype, minval(ssh),           maxval(ssh)

  call compute_diag
end subroutine ocean_step
!=================================================================
