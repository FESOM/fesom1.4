module g_config
  implicit none
  save

  ! *** Modelname ***
  character(5)             	:: runid='test1'                ! a model/setup name

  namelist /modelname/ runid

  ! *** time step ***
  integer                  	:: step_per_day=12           	!number of steps per day
  integer                  	:: run_length=1	                !run length
  character                     :: run_length_unit='y'          !unit: y, d, s

  namelist /timestep/ step_per_day, run_length, run_length_unit

  ! *** Paths for all in and out ***
  character(2000)                :: MeshPath='./mesh/'
  character(2000)                :: OpbndPath='./opbnd/'
  character(2000)                :: ClimateDataPath='./hydrography/'
  character(2000)                :: ForcingDataPath='./forcing/'
  character(2000)                :: TideForcingPath='./tide_forcing/'
  character(2000)                :: ResultPath='./result/'
  character(2000)                :: ECHAM6ForcingDataPath='./cpl_work/'

  namelist /paths/  MeshPath, OpbndPath, ClimateDataPath, ForcingDataPath, &
       TideForcingPath, ResultPath, ECHAM6ForcingDataPath

  ! *** ocean climatology data name ***
  character(100)                :: OceClimaDataName='annual_woa01_ts.out'
  logical                       :: use_prepared_init_ice=.false.     !how to initial. ice at the beginning 

  namelist /initialization/ OceClimaDataName, use_prepared_init_ice

  ! *** in out ***
  character*4               	:: restartflag='last'  	             !restart from which saved record,'#','last'
  integer                       :: output_length=1                   !valid for d,h,s
  character                	:: output_length_unit='m'      	     !output period: y, m, d, h, s
  character                	    :: output_length_unit_restart='m'    !output period: y, m, d, h, s  
  integer                       :: logfile_outfreq=1                 !in logfile info. output frequency, # steps
  logical :: levelwise_output = .false.

  namelist /inout/ restartflag, output_length, output_length_unit, output_length_unit_restart, logfile_outfreq, levelwise_output

  ! *** mesh ***
  integer                       :: grid_type=1              	! z-level, 2 sigma, 3 sigma + z-level

  namelist /mesh_def/ grid_type

  ! *** model geometry
  logical                  	:: cartesian=.false.
  logical                  	:: fplane=.false.
  logical                  	:: betaplane=.false.
  real(kind=8)             	:: f_fplane=-1.4e-4        	![1/s]
  real(kind=8)             	:: beta_betaplane=2.0e-11  	![1/s/m]
  real(kind=8)             	:: domain_length=360.    	![degree]
  !
  logical                  	:: rotated_grid=.true.    	!option only valid for coupled model case now
  logical                  	:: force_rotation=.false.	!set to .true. for some unrotated meshes
  real(kind=8)             	:: alphaEuler=50. 		![degree] Euler angles, convention:
  real(kind=8)             	:: betaEuler=15.  		![degree] first around z, then around new x,
  real(kind=8)			:: gammaEuler=-90.		![degree] then around new z.

  namelist /geometry/  cartesian, fplane, betaplane, f_fplane, beta_betaplane, &
       domain_length, rotated_grid, force_rotation, alphaEuler, betaEuler, gammaEuler

  ! *** fleap_year ***
  logical                       :: include_fleapyear=.false.
  
  namelist /calendar/ include_fleapyear
  
   ! *** machine ***
  integer                       :: system=2                     ! XD1 2(byte), HLRN 1(word)
  
  namelist /machine/ system


  ! *** others ***
  real(kind=8)             	:: dt, dt_inv
  integer                  	:: istep, nsteps
  integer                       :: save_count, save_count_restart, seconds_from_yearstart
  logical                       :: r_restart
  logical                       :: blowed_up=.false.

end module g_config
