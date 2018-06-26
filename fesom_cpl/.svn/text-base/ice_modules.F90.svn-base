module i_dyn_parms
  implicit none
  save

  ! *** ice_stress ***
  real(kind=8)             :: Pstar = 30000.           ! [N/m^2] 15000. 27500. 30000.
  real(kind=8)             :: ellipse =2.              !
  real(kind=8)             :: c_pressure =20.0         !
  real(kind=8)             :: delta_min=1.0e-11        ! [s^(-1)]

  namelist /ice_stress/ Pstar, ellipse, c_pressure, delta_min

  ! *** friction setting ***
  real*8                   :: Cd_oce_ice = 5.5e-3      ! drag coefficient ocean-ice 3.e-3 5.5e-3
  real*8                   :: Kh_ice=0.0               ! numerical ice/snow diffusivity

  namelist /ice_fric/ Cd_oce_ice, Kh_ice

  ! *** rheology ***
  logical                  :: EVP_rheology=.true.
  integer                  :: evp_rheol_steps=120      ! substeps for EVP rheology
  integer                  :: evp_Tdamp_ratio=3        ! ratio dt/T_damp
  integer                  :: vp_rheol_steps=500       ! substeps for VP rheology

  namelist /ice_rheology/  EVP_rheology, evp_rheol_steps, evp_Tdamp_ratio, vp_rheol_steps

  ! *** numerical scheme ***
  real(kind=8)             :: ice_gamma_fct=0.1
  logical                  :: lump_ice_matrix=.true.   !for mass ice matrix case
  integer                  :: num_iter_solve_ice=3 

  namelist /ice_scheme/ ice_gamma_fct

  ! *** others ***
  real(kind=8)             :: Tevp_inv

end module i_dyn_parms

!=======================================================================

module i_therm_parms
  implicit none
  real*8, parameter  :: rhoair=  1.3      ! Air density ! AOMIP
  real*8, parameter  :: rhowat= 1025.     ! Water density
  real*8, parameter  :: rhoice=  910.     ! Ice density
  real*8, parameter  :: rhosno=  290.     ! Snow density

  real*8, parameter  :: cpair=1005.       ! Specific heat of air [J/(kg * K)] / 1004 
  real*8, parameter  :: cc=rhowat*4190.0  ! Volumetr. heat cap. of water [J/m**3/K](cc = rhowat*cp_water)
  real*8, parameter  :: cl=rhoice*3.34e5  ! Volumetr. latent heat of ice fusion [J/m**3](cl=rhoice*Lf) 
  real*8, parameter  :: clhw=2.501e6      ! Specific latent heat [J/kg]: water	-> water vapor
  real*8, parameter  :: clhi=2.835e6      !                              sea ice-> water vapor

  real*8, parameter  :: tmelt=273.15      ! 0 deg C expressed in K 
  real*8, parameter  :: boltzmann=5.67E-8 ! S. Boltzmann const.*longw. emissivity

  real*8, parameter  :: con   = 2.1656    ! Thermal conductivities: ice; W/m/K
  real*8, parameter  :: consn = 0.31	  !                         snow

  real*8, parameter  :: hmin= 0.05        ! Cut-off ice thickness     !!
  real*8, parameter  :: Armin=0.10        ! Minimum ice concentration !!

  integer,parameter  :: iclasses=7        ! Number of ice thickness gradations for ice growth calcs.

  real*8    :: Sice = 4.0        ! Ice salinity 3.2--5.0 ppt.
  real*8    :: h0=0.5	         ! Lead closing parameter [m] ! 0.5

  real*8    :: emiss_ice=0.97    ! Emissivity of Snow/Ice
  real*8    :: emiss_wat=0.97    ! of open water  

  real*8    :: albsn=	0.81 	 ! Albedo: frozen snow
  real*8    :: albsnm=  0.77   	 !         melting snow
  real*8    :: albi=	0.70 	 !         frozen ice
  real*8    :: albim=	0.68 	 !         melting ice
  real*8    :: albw=	0.1	 !         open water

  namelist /ice_therm/ Sice, h0, emiss_ice, emiss_wat, albsn, albsnm, albi, albim, albw

end module i_therm_parms

!=======================================================================

module i_array
  use o_data_types
  implicit none
  save
  real(kind=8), allocatable, dimension(:)         :: u_ice, v_ice
  real(kind=8), allocatable, dimension(:)         :: m_ice, a_ice, m_snow  
  real(kind=8), allocatable, dimension(:)         :: rhs_u, rhs_v
  real(kind=8), allocatable, dimension(:)         :: rhs_m, rhs_a, rhs_ms
  real(kind=8), allocatable, dimension(:)         :: dm_ice, da_ice, dm_snow
 
  real(kind=8), allocatable, dimension(:)         :: t_skin

  type(sparse_matrix)                             :: icestiff
  real(kind=8), allocatable, dimension(:)         :: ice_lump

  real(kind=8), allocatable, dimension(:)         :: m_icel, a_icel, m_snowl
  real(kind=8), allocatable, dimension(:)         :: icepplus, icepminus 
  real(kind=8), allocatable, dimension(:,:)       :: icefluxes

  real(kind=8), allocatable, dimension(:)         :: sigma11, sigma12, sigma22  

  real(kind=8), allocatable, dimension(:)         :: u_w, v_w
  real(kind=8), allocatable, dimension(:)         :: fresh_wa_flux
  real(kind=8), allocatable, dimension(:)         :: net_heat_flux
#if defined (__oasis) || defined (__uncplecham6)   
  real(kind=8),target, allocatable, dimension(:)  :: oce_heat_flux, ice_heat_flux  
  real(kind=8),target, allocatable, dimension(:)  :: tmp_oce_heat_flux, tmp_ice_heat_flux 
							!temporary flux fields
							!(for flux correction)
#endif /* (__oasis) || (__uncplecham6)*/
  real(kind=8), allocatable, dimension(:)         :: S_oc_array, T_oc_array
  real(kind=8), allocatable, dimension(:)         :: elevation

  real(kind=8), target, allocatable, dimension(:) :: stress_iceoce_x         
  real(kind=8), target, allocatable, dimension(:) :: stress_iceoce_y
  real(kind=8), target, allocatable, dimension(:) :: stress_atmice_x         
  real(kind=8), target, allocatable, dimension(:) :: stress_atmice_y
  real(kind=8), target, allocatable, dimension(:) :: stress_atmoce_x         
  real(kind=8), target, allocatable, dimension(:) :: stress_atmoce_y

end module i_array

!=====================================================================

module i_solver
implicit none
save
  integer                  	:: solve_m_ice=101
  integer       	 	:: solve_a_ice=102
  integer			:: solve_m_snow=103
end module i_solver
