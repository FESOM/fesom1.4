subroutine forcing_array_setup
  !inializing forcing fields 
  use o_param
  use o_mesh
  use i_array
  use g_forcing_param
  use g_forcing_arrays
  use g_parfe
  use g_config
#if defined (__oasis) || defined (__uncplecham6)
  use cpl_driver, only : nrecv
#endif   
  implicit none

  integer    :: n2

  n2=myDim_nod2D+eDim_nod2D  

  ! Allocate memory for atmospheric forcing 
  allocate(shortwave(n2), longwave(n2))
  allocate(prec_rain(n2), prec_snow(n2))
  allocate(u_wind(n2), v_wind(n2))
  allocate(Tair(n2), shum(n2))
  allocate(runoff(n2), evaporation(n2))

#if defined (__oasis) || defined (__uncplecham6)
  allocate(sublimation(n2), evap_no_ifrac(n2))
  allocate(tmp_sublimation(n2),tmp_evap_no_ifrac(n2), tmp_shortwave(n2))
  allocate(atm_net_fluxes_north(nrecv), atm_net_fluxes_south(nrecv))
  allocate(oce_net_fluxes_north(nrecv), oce_net_fluxes_south(nrecv))
  allocate(flux_correction_north(nrecv), flux_correction_south(nrecv))
  allocate(flux_correction_total(nrecv))
  sublimation=0.
  evap_no_ifrac=0.
  tmp_sublimation = 0.
  tmp_evap_no_ifrac = 0.
  tmp_shortwave = 0.
  atm_net_fluxes_north=0.
  atm_net_fluxes_south=0.
  oce_net_fluxes_north=0.
  oce_net_fluxes_south=0.
  flux_correction_north=0.
  flux_correction_south=0.
  flux_correction_total=0.  
#endif 

  shortwave=0.
  longwave=0.
  prec_rain=0.
  prec_snow=0.
  u_wind=0.
  v_wind=0.
  Tair=0.
  shum=0.
  runoff=0.
  evaporation=0.

  if(use_landice_water) then
     allocate(runoff_landice(n2))
     runoff_landice=0.0
  end if
 
  ! shortwave penetration
#ifdef use_sw_pene
  allocate(chl(n2))
  allocate(sw_3d(myDim_nod3d+eDim_nod3D))
  chl=0.0
#endif

  !for ice diagnose
#ifdef use_ice
  allocate(thdgr(n2), thdgrsn(n2), flice(n2))
  allocate(olat_heat(n2), osen_heat(n2), olwout(n2))
  thdgr=0.
  thdgrsn=0.
  flice=0.
  olat_heat=0.
  osen_heat=0.
  olwout=0.
#endif 

  ! drag coefficient and transfer coefficients for latent and sensible heat
  allocate(Cd_atm_oce_arr(n2))      
  allocate(Ce_atm_oce_arr(n2))
  allocate(Ch_atm_oce_arr(n2))
  Cd_atm_oce_arr=Cd_atm_oce
  Ce_atm_oce_arr=Ce_atm_oce 
  Ch_atm_oce_arr=Ch_atm_oce
#ifdef use_ice
  allocate(Cd_atm_ice_arr(n2)) 
  Cd_atm_ice_arr=Cd_atm_ice   
#endif

  if(mype==0) write(*,*) 'forcing arrays have been set up'   

end subroutine forcing_array_setup
!
!----------------------------------------------------------------------
!
subroutine forcing_array_setup_OnlyOcean
  !inializing forcing fields for an ocean-alone case
  !currently only wind is applied.
  use o_param
  use o_mesh
  use i_array
  use g_forcing_param
  use g_forcing_arrays
  use g_parfe
  implicit none

  integer    :: n2

  n2=myDim_nod2D+eDim_nod2D  

  ! Allocate memory for atmospheric forcing 
  allocate(u_wind(n2), v_wind(n2))
  u_wind=0.
  v_wind=0.

  ! drag coefficient and transfer coefficients for latent and sensible heat
  allocate(Cd_atm_oce_arr(n2))      
  Cd_atm_oce_arr=Cd_atm_oce

  if(mype==0) write(*,*) 'forcing arrays (for an ocean-alone case) have been set up'   

end subroutine forcing_array_setup_OnlyOcean
