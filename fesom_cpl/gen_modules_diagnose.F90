! modules of diagnose control paramers and diagnose arrays

module g_diag
  ! diagnose flag
  implicit none
  save

  ! run state
  logical                                  :: check_run_state=.true.       ! use salinity to check blowup

  ! ocean
  logical                                  :: diag_oce=.true.
  logical                                  :: diag_oce_KE=.true.
  logical                                  :: diag_oce_energy_conv=.true.
  logical                                  :: diag_oce_mix_layer=.true.
  logical                                  :: diag_oce_transp=.true.
  logical                                  :: diag_oce_GM_vel=.true.
  logical                                  :: diag_oce_SGS_transp=.true.
  logical                                  :: diag_oce_Kv=.true.
  
  ! ice
  logical                                  :: diag_ice=.true.
  ! forcing
  logical                                  :: diag_forcing=.true.
  ! mesh
  logical                                  :: diag_mesh=.true.


  namelist /diag_flag/ check_run_state, &
       diag_oce, diag_oce_KE, diag_oce_energy_conv, &
       diag_oce_mix_layer, diag_oce_transp, &
       diag_oce_GM_vel, diag_oce_SGS_transp, diag_oce_Kv, &
       diag_ice, diag_forcing, diag_mesh

end module g_diag
!
!------------------------------------------------------------------------------
!
module g_meanarrays
  use g_meandata, only : Meandata2D, Meandata3D, all_meandata, LocalMeandata, Snapshotdata2D
  ! mean and (mean)diagnose arrays
  implicit none
  save

  ! counter
  integer                                  :: meancounter

  ! mean of prediction variables

  ! ocean
  type(Meandata3D), target, allocatable, dimension(:) :: passive_tracers_mean, age_tracers_mean
  type(Meandata3D), target :: thetao_mean, so_mean
  type(Meandata2D), target :: zos_mean
  type(Meandata3D), target :: uo_mean, vo_mean, wo_mean

  ! ice
  type(Meandata2D), target :: sisnthick_mean
  type(Meandata2D), target :: ice_area_fraction_mean, sea_ice_thickness_mean, sea_ice_volume_mean
  type(Meandata2D), target :: sea_ice_x_velocity_mean, sea_ice_y_velocity_mean, sea_ice_speed_mean
  type(Meandata2D), target :: sea_ice_time_fraction_mean

  ! (mean) diagnose variables

  ! ocean
  type(Meandata3D), target :: uto_mean, vto_mean
  type(Meandata3D), target :: uso_mean, vso_mean
  real(kind=8), allocatable, dimension(:)  :: sgs_u, sgs_v
  real(kind=8), allocatable, dimension(:)  :: sgs_ut, sgs_vt
  real(kind=8), allocatable, dimension(:)  :: sgs_us, sgs_vs
  type(Meandata2D), target :: mlotst_mean
  type(Meandata3D), target :: u2o_mean, v2o_mean
  type(Meandata3D), target :: rho_mean, urho_mean
  type(Meandata3D), target :: vrho_mean, uv_mean

  ! ice
  type(Meandata2D), target :: thdgr_mean, thdgrsn_mean
  type(Meandata2D), target :: uhice_mean, vhice_mean
  type(Meandata2D), target :: uhsnow_mean, vhsnow_mean
  type(Meandata2D), target :: flice_mean, fsitherm_mean
  type(Meandata2D), target :: sisnconc_mean

  ! forcing
#ifndef __oasis
  type(Meandata2D), target :: tair_mean, shum_mean, lwrd_mean, olat_mean, osen_mean, olwout_mean, relax_salt_mean
#endif
  type(Meandata2D), target :: uwind_mean, vwind_mean
  type(Meandata2D), target :: prlq_mean, prsn_mean
  type(Meandata2D), target :: runoff_mean, evap_mean
  type(Meandata2D), target :: rsdo_mean
  type(Meandata2D), target :: hfds_mean, wnet_mean
  type(Meandata2D), target :: virtual_salt_mean
  type(Meandata2D), target :: tauuo_mean, tauvo_mean
  
  ! PRIMAVERA data request ocean and ice
  type(Meandata3D), target :: w2o_mean, wso_mean, wto_mean, opottemptend_mean
  type(Meandata2D), target :: omldamax_mean, zossq_mean, pbo_mean, tos_mean, sos_mean, sistrxdtop_mean, sistrydtop_mean, wfo_mean
  type(Snapshotdata2D), target :: tso_mean
#if defined (__oasis)
  type(Meandata2D), target :: evs_mean, sidmassevapsubl_mean, sifllatstop_mean
#endif
  type(LocalMeandata), target :: volo_const, soga_mean, thetaoga_mean, siarean_mean, siareas_mean, siextentn_mean, siextents_mean, sivoln_mean, sivols_mean

  ! PRIMAVERA data request ice
  type(Meandata2D), target :: sisnmass_mean, sidmassth_mean, sidmasssi_mean, sidmasstranx_mean, sidmasstrany_mean, sistrxubot_mean, sistryubot_mean
end module g_meanarrays
!
!----------------------------------------------------------------------------

