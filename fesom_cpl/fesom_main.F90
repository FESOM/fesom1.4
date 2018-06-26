subroutine print_define_info
#ifdef use_sw_pene
      print *, 'short wave penetration is ON'
#else
      print *, 'short wave penetration is OFF'
#endif  
#ifdef __oasis
      print *, 'coupled (oasis coupling is ON)'
#else
      print *, 'standalone (oasis coupling is OFF)'
#endif  
end

program main
  !=============================================================================!
  !
  !                     Finite Element Ocean Model
  !
  !=============================================================================!
  !                      The main driving routine
  !=============================================================================!         

  use g_meanarrays, only : all_meandata
  use g_config
  use o_param
  use o_array          
  use o_solver
  use o_mixing_kpp_mod
  use g_parfe
  use g_clock
  use g_forcing_index
  use g_forcing_param
  use g_forcing_arrays
  use g_diag
  use g_forcing_interp
#ifdef use_ice
  use i_array
#endif
#if defined (__oasis)
  use cpl_config_param
  use cpl_driver
  use cpl_exchange_mesh
!  USE PSMILe, only  : comm_psmile    
  use mod_prism
#endif
  use fesom_env_module, only : fesom_env__check_namelist_config
  implicit none
  integer :: i
  real(kind=8) :: t0, t1, wallclock_initial, wallclock_current, last_wallclock_current = 0.0
  integer :: wallclock_hrs, wallclock_min, wallclock_sec
  integer      :: ioerror
  character(len=:), allocatable :: arg
  integer arglength
  
  do i = 1, command_argument_count()    
    call get_command_argument(i, length=arglength)
    allocate(character(arglength) :: arg)
    call get_command_argument(i, value=arg)
    select case (arg)
    case ('--smoketest')
      print *, 'smoketest'
      stop
    case ('--info')
      call print_define_info
      stop
    case ('--check')
      call fesom_env__check_namelist_config
      stop
    case ('--defaults')
      call print_default_schedules
      stop
    case default
      print *, 'unknown option: ', arg
      error stop
    end select
    deallocate(arg)
  end do
  
  ! MPI initialization
  call par_init                 ! initializes MPI
  
  wallclock_initial = MPI_Wtime()

  if(mype==0) then
    call print_define_info
    call fesom_env__check_namelist_config
    write(*, *) 'Running on ', npes, ' PEs'
    write(*,*) 'using these schedules to write output:'
    call print_schedules()
    write(*,*) '*************************************************************'
  end if
  ! read namelist, initialize clock, prepare basic configuration etc.
  call setup_model              ! setup basic config, do it before clock_init
  call clock_init               ! read the clock file
  call get_run_steps            ! compute total steps to run 
  call config_remind_warning_info
  if(mype==0) write(*,*) '*************************************************************'
  ! mesh and communication buffers

  call ocean_mesh_setup         ! setup the 2D/3D meshes
  call set_par_support
  call mesh_cluster_setup       ! build cluster area and volume and save to file

  if(mype==0) write(*,*) '*************************************************************'
  ! ocean: matrices, arrays, initialization, buffer zone, tide etc.

  call ocean_matrices_setup     ! Builds matrices and call partitioning
  call ocean_array_setup        ! allocate ocean arrays 
  call ocean_init               ! initialize the oce or read restart files

  if(use_ref_density) then
     call compute_ref_density   ! Fills in ocean reference density 
  endif

#ifdef use_opbnd_restoring
  call init_restoring_vel
#endif

  if(buffer_zone) then
     call init_restoring_bufferzone
  end if

#ifdef use_opbnd_tide
  call init_tidal_opbnd         ! initialize tidal ocean open boundary
#endif

#ifdef use_ice
  if(mype==0) write(*,*) '*************************************************************'
  ! ice: matrices, arrays, initialization

  call ice_matrices_setup       ! Build ice matrices
  call ice_array_setup          ! allocate ice arrays, setup ice adv matrix
  call ice_init                 ! initialize the ice or read restart files
#endif


  if(mype==0) write(*,*) '*************************************************************'
  ! forcing: arrays, initialization, interpolation preparation  

#ifdef use_ice
  call forcing_array_setup  
  call init_forcing_interp      ! calculates the forcing interpolation weights
  call init_atm_forcing         ! initialize forcing fields
#else
#if not defined(toy_ocean) && not defined(__oasis)
  call forcing_array_setup_OnlyOcean
  call init_forcing_interp 
  call init_atm_forcing_OnlyOcean 
#endif 
#endif

  if(use_landice_water) call landice_water_init

#if defined (__oasis)
#if defined cpl_str
      call cpl_oasis3mct_define_unstr
if(mype==0)  write(*,*) '---->     cpl_oasis3mct_define_unstr nsend, nrecv:',nsend, nrecv
#else
    call cpl_oasis_define
    if(mype==0)  write(*,*) '---->     cpl_oasis_define nsend, nrecv:',nsend, nrecv
#endif /* (defined cpl_str)  */
#endif /* (__oasis)         */



  if(mype==0) write(*,*) '*************************************************************'
  ! init mean arrays, adjust clock and create output files if required  

#if defined(allow_calcmeans) || defined(allow_diag)
  call init_meanarrays          ! allocate arrays for mean fields
#endif
  call clock_newyear  		! check if it is a new year
  call init_output              ! create new output files


#ifdef use_fullfreesurf
  if(mype==0) write(*,*) '*************************************************************'
  ! updating mesh and matrices in case of full free surface setup 

  if(any(abs(ssh)>small)) then
     call update_mesh
     call update_matrices
     call update_mesh
  endif
#endif

  ! set some flags for solvers
  iter_first=.true.             ! iter_first & iteruv_first should be 'true' at start
  iteruv_first=.true.


  if(mype==0) then
     write(*,*) '*************************************************************'
     write(*,*) 'iteration starts ...'
     write(*,*) '*************************************************************'
     write(*,*)
  end if

  ! preparation done
  !----------------------------------------------------------------------------
  ! start iteration
  if(mype==0) call stopwatch_start_id(mype, 'mainloop'//CHAR(0))

  do istep=1, nsteps
!if (mype==0) write(*,*) 'reprotest  SSH:', istep, minval(ssh, 1), maxval(ssh, 1)
!if (mype==0) write(*,*) 'reprotest TEMP:', istep, minval(tracer(:,1), 1), maxval(tracer(:,1), 1)    
#ifdef __oasis
     seconds_til_now=INT(dt)*(istep-1)
#endif
     seconds_from_yearstart=seconds_from_yearstart+int(dt) !for writing output
     t0=MPI_Wtime()   
     call clock
     call init_output
     call forcing_index 
#ifdef use_ice    
     call ocean2ice
     call update_atm_forcing
     call ice_step
     call ice2ocean
     if(use_landice_water) call add_landice_water
#else
#ifndef toy_ocean
     call update_atm_forcing_OnlyOcean
#endif
#endif
#ifdef use_fullfreesurf 
     if(balance_salt_water) call check_imb_freshwater
#endif
#ifdef use_sw_pene
     call cal_shortwave_rad
#endif
     call ocean_step
#if defined(allow_calcmeans) || defined(allow_diag)
     call add2meanarrays      
#endif
     ! save (NetCDF)
     call output()
!!$     ! save (ascii, for an easy debugging during new setup test phase)
!!$     if(mod(istep,step_per_day*5)==0) then
!!$        call oce_out     
!!$#ifdef use_ice
!!$        call ice_out
!!$#endif
!!$     end if

      if(check_run_state) call check_blowup	! check if the program blows up
      t1=MPI_Wtime()
      if(mod(istep,logfile_outfreq)==0) then	
        ! log file output (for debugging during new setup test phase)
        !write(*,*) 'uf max', maxval(abs(uf)), istep
        !write(*,*) 's min', minval(tracer(:,2)), istep

        if(mype==0) then
          wallclock_current = MPI_Wtime() - wallclock_initial
          wallclock_hrs = int(wallclock_current/3600.0)
          wallclock_sec = int(wallclock_current)-wallclock_hrs*3600
          wallclock_min = int(real(wallclock_sec,8)/60.0)
          wallclock_sec = wallclock_sec-wallclock_min*60          
          write(*,'(a,i6,a,i3,a,i0,a,i0,f0.5,a,i0,f0.5,a,i0.2,a,i0.2,a,i0.2)') 'step: ', istep, ', day: ', daynew, ', year: ', yearnew, ', duration: ', int(t1-t0),t1-t0-int(t1-t0), 's, total wallclock: +', int(wallclock_current-last_wallclock_current),wallclock_current-last_wallclock_current-int(wallclock_current-last_wallclock_current), 's = ', wallclock_hrs, ':', wallclock_min, ':', wallclock_sec
          last_wallclock_current = wallclock_current          
        end if
      end if
  end do

  ! iteration done
  !--------------------------------------------------------------------
  ! some finishing-up routines follow
  do i=1, size(all_meandata)
    call all_meandata(i)%ptr%finalize
  end do
  if(mype==0) call stopwatch_stop()
  if(mype==0) call stopwatch_print_all()

  !#ifndef use_ice
  ! call close_ocean_forcing	! required in special cases
  !#endif

  if(mix_scheme=='MY2p5') call save_MY_vara  ! save MY2.5 variables for next restart

  ! save (ascii, for an easy debugging during new setup test phase)
  !call oce_out
#ifdef use_ice
  !call ice_out
#endif

  call clock_finish		! save clock

  if(mype==0) then
     open(unit=50, file='goodfile')
     write(50,*)'go on'
     close(50)
  end if
   
#ifdef VERBOSE
  write(*,*) 'Before call to par_ex'
#endif
  call par_ex 			! finalizes MPI
  write(*,*) 'Experiment '//runid//' successfully completed'

end program main
