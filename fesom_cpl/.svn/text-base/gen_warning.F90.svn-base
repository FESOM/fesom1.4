subroutine config_remind_warning_info
  ! check configuration and options and provide reminding or warning
  use o_param
  use o_array          
  use o_solver
  use g_config
  use g_parfe
  use g_clock
  use g_forcing_param
  use g_forcing_index
  use g_forcing_arrays
  implicit none

  if(mype/=0) return

#ifdef use_opbnd_restoring
  write(*,*) '---------------------------------------------------'
  write(*,*) 'Reminding:'
  write(*,*) 'It is configured to restore velocity at open boundarie.'
  write(*,*) 'Routine init_restoring_vel should be specifically prepared.'
  write(*,*) '---------------------------------------------------'
#endif

  if(buffer_zone) then
     write(*,*) '---------------------------------------------------'
     write(*,*) 'Reminding:'
     write(*,*) 'It is configured to restore tracers at open boundaries.'
     write(*,*) 'Routine init_restoring_bufferzone should be specifically'
     write(*,*) 'prepared.'
     write(*,*) '---------------------------------------------------'
  end if

#if defined(use_opbnd_tide) || defined(use_opbnd_restoring)
  if(.not.buffer_zone) then
     write(*,*) '---------------------------------------------------'
     write(*,*) 'Warning:'
     write(*,*) 'It is recommended to restore tracers in a buffer zone near o.b.'
     write(*,*) '---------------------------------------------------'
  end if
#endif

#ifdef use_opbnd_tide
#ifndef use_semiimplicit_scheme
  write(*,*) '---------------------------------------------------'
  write(*,*) 'Warning:'
  write(*,*) 'It is recommended to use the semi-implicit scheme when'
  write(*,*) 'simulating tides. To change, set it in the Makefile.'
  write(*,*) '---------------------------------------------------'
#endif
#endif

#ifdef use_ice
  if(restore_s_surf>0.0) then
     write(*,*) '---------------------------------------------------'
     write(*,*) 'Reminding:'
     write(*,*) 'It is specified to restore SSS. Check which climatology is'
     write(*,*) 'to be used (yearly or monthly, which source) and modify the'
     write(*,*) 'code in file gen_forcing_couple.F90 for your purpose.'
     write(*,*) 'Also check how restoring under ice is done in the code;'
     write(*,*) 'you might need to change the code for your application.'
     write(*,*) '---------------------------------------------------'
  end if
#endif

  if(balance_salt_water) then
     write(*,*) '---------------------------------------------------'
     write(*,*) 'Reminding:'
     write(*,*) 'You specified to balance salt and freshwater fluxes.'
     write(*,*) 'The default is to do correction for every time step.'
     write(*,*) 'You might need your own correction strategy. Check'
     write(*,*) 'the code to be sure it does what you want.'
     write(*,*) '---------------------------------------------------'
   end if  

#ifdef use_sw_pene
  write(*,*) '---------------------------------------------------'
  write(*,*) 'Reminding:'
  write(*,*) 'It is specified to consider shortwave penetration. The'
  write(*,*) 'chlorophyll climatology on meshes should have been prepared'
  write(*,*) 'offline (unformated). Automatic interpolation to model grids is'
  write(*,*) 'not supported in the code!'
  write(*,*) '---------------------------------------------------'
#endif

#ifndef use_ice
  write(*,*) '---------------------------------------------------'
  write(*,*) 'Warning:'
  write(*,*) 'You are running the ocean-alone model. The surface forcing routine'
  write(*,*) 'should be adjusted/checked to properly apply your particular surface'
  write(*,*) 'forcing to the ocean.'
  write(*,*) '---------------------------------------------------'
#endif

  if(wind_data_source=='NCEP' .and. ncar_bulk_formulae) then
      write(*,*) '---------------------------------------------------'
      write(*,*) 'Warning:'
      write(*,*) 'you are using NCEP 2m air temperature and humidity. The current'
      write(*,*) 'formulae for calculating drag and heat exchange coefficients'
      write(*,*) 'only support 10m data. If you plan to use these formulae, a small'
      write(*,*) 'update is required. Before doing it, turn off ncar_bulk_formulae!'
      write(*,*) '---------------------------------------------------'
   end if

  if(mix_scheme=='MY2p5') then
     write(*,*) '---------------------------------------------------'
     write(*,*) 'Reminding:'
     write(*,*) 'MY2.5 mixing scheme is to be used. Its variables will'
     write(*,*) 'be saved at the end of the run for the next restart job.'
     write(*,*) 'This file will be replaced in each restart run. So it'
     write(*,*) 'is not possible to restart from intermediate snapshots,'
     write(*,*) 'or from previous runs if the file is not backed up'
     write(*,*) 'manually. To get rid of this limit, code need update.'
     write(*,*) '---------------------------------------------------'
  end if

  if(use_passive_tracer) then
     if(ptr_start_year>yearnew) then
        write(*,*) '---------------------------------------------------'
        write(*,*) 'Warning:'
        write(*,*) 'You specify to use passive tracers not at the beginning of'
        write(*,*) 'this job. This is not supported. The model will start to'
        write(*,*) 'include the passive tracers from the beginning of this job.'
        write(*,*) 'If you do not want this, cancel this job, turn off'
        write(*,*) 'use_passive_tracer, run the model until the moment when'
        write(*,*) 'you want to have the passive tracers, and then re-start the'
        write(*,*) 'simulation with passive tracers used.'
        write(*,*) '---------------------------------------------------'
     end if
  end if

  if(use_age_tracer) then
     if(age_tracer_start_year>yearnew) then
        write(*,*) '---------------------------------------------------'
        write(*,*) 'Warning:'
        write(*,*) 'You specify to use age tracers not at the beginning of'
        write(*,*) 'this job. This is not supported. The model will start to'
        write(*,*) 'include the age tracers from the beginning of this job.'
        write(*,*) 'If you do not want this, cancel this job, turn off'
        write(*,*) 'use_age_tracer, run the model until the moment when'
        write(*,*) 'you want to have the age tracers, and then re-start the'
        write(*,*) 'simulation with age tracers used.'
        write(*,*) '---------------------------------------------------'
     end if
  end if
end subroutine config_remind_warning_info
!
!--------------------------------------------------------------------------
!
subroutine check_blowup
  !check if the model blows up and cancel the job if it is the case
  !Salinity is used as an indicator.
  !ALLREDUCE is used, so it slows down the code. 
  !A better way is requried!!! Qiang
  use o_array          
  use g_config
  use g_parfe
  use o_MESH
  use o_param
  implicit none

  integer     :: flag, g_flag
  integer     :: n

  flag=0
  if(any(tracer(:,2)<0.0)) then
     flag=1
     write(*,*) 'Step: ', istep
     do n=1, todim_nod3D
	if (tracer(n,2)<0.0) then
		write(*,*) '***************************************************'
		write(*,*) 's<0 at nod3D=', mylist_nod3D(n)
		write(*,*) 'x=', coord_nod3D(1,n)/rad, ' ; ', 'y=', coord_nod3D(2,n)/rad
		write(*,*) 'z=', coord_nod3D(3,n)
		write(*,*) '***************************************************'
	end if
     end do
  end if
   
  g_flag=0
  call MPI_AllREDUCE(flag, g_flag, 1, MPI_INTEGER, MPI_SUM, &
       MPI_COMM_FESOM, MPIerr)
  
  if(g_flag>0) then
     if(mype==0) then
        write(*,*) 'Negative salinity found. The model blows up!'
        write(*,*) 'The program will be forced to stop here.'
     end if
     blowed_up=.true.
     call par_ex
     stop
  end if
end subroutine check_blowup
