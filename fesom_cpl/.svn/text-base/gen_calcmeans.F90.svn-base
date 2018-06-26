subroutine init_meanarrays
  ! allocates and initializes fields used for mean arrays
  use o_mesh
  use o_param
  use g_config
  use g_diag
  use g_meanarrays
  use g_parfe

  implicit none
  integer       :: n2, n3, i
  character(1) :: tracer_suffix

  n2=myDim_nod2D+eDim_nod3D        
  n3=myDim_nod3D+eDim_nod3D       

  ! for mean fields

#ifdef allow_calcmeans

  ! ocean part                    
  call volo_const%init(1, 'volo', 'total volume of liquid seawater', 'm3') ! this could be optimized, as the value does not change. we use the standard output writer procedure for simplicity
  call soga_mean%init(1, 'soga', 'global average sea water salinity', 'psu')
  call thetaoga_mean%init(1, 'thetaoga', 'global average sea water potential temperature ', 'degC')
  call siarean_mean%init(1, 'siarean', 'total area of sea ice in the Northern hemisphere', '1e6 km2')
  call siareas_mean%init(1, 'siareas', 'total area of sea ice in the Southern hemisphere', '1e6 km2')
  call siextentn_mean%init(1, 'siextentn', 'total area of all Northern-Hemisphere grid cells that are covered by at least 15 % areal fraction of sea ice', '1e6 km2')
  call siextents_mean%init(1, 'siextents', 'total area of all Southern-Hemisphere grid cells that are covered by at least 15 % areal fraction of sea ice', '1e6 km2')
  call sivoln_mean%init(1, 'sivoln', 'total volume of sea ice in the Northern hemisphere', '1e3 km3')
  call sivols_mean%init(1, 'sivols', 'total volume of sea ice in the Southern hemisphere', '1e3 km3')

  call zos_mean%init(n2, 'zos', 'dynamic sea level', 'm')
  call zossq_mean%init(n2, 'zossq', 'sea surface height squared. Surface ocean geoid defines z=0.', 'm2')
  call omldamax_mean%init(n2, 'omldamax', 'daily maximum ocean mixed layer thickness defined by mixing scheme', 'm')
  call pbo_mean%init(n2, 'pbo', 'pressure at ocean bottom', 'Pa') ! sea_water_pressure_at_sea_floor
  call tos_mean%init(n2, 'tos', 'sea surface temperature of liquid ocean', 'degC')
  call tso_mean%init(n2, 'tso', 'sea surface temperature of liquid ocean, sampled synoptically', 'K')
  call sos_mean%init(n2, 'sos', 'sea surface salinity ', 'psu')
#if defined (__oasis)
  call evs_mean%init(n2, 'evs', 'computed as the total mass of water vapor evaporating from the ice-free portion of the ocean  divided by the area of the ocean portion of the grid cell', 'kg m-2 s-1')
  call sidmassevapsubl_mean%init(n2, 'sidmassevapsubl', 'The rate of change of sea-ice mass change through evaporation and sublimation divided by grid-cell area', 'kg m-2 s-1')
  call sifllatstop_mean%init(n2, 'sifllatstop', 'Dummy: Not computed in FESOM, rather in ECHAM. the net latent heat flux over sea ice', 'W m-2')
#endif
  call sistrxdtop_mean%init(n2, 'sistrxdtop', 'x-component of atmospheric stress on sea ice', 'N m-2')
  call sistrydtop_mean%init(n2, 'sistrydtop', 'y-component of atmospheric stress on sea ice', 'N m-2')
  call wfo_mean%init(n2, 'wfo', 'computed as the water  flux into the ocean divided by the area of the ocean portion of the grid cell.  This is the sum of the next two variables in this table', 'kg m-2 s-1')

  call uo_mean%init(n3, 'uo', 'Prognostic x-ward velocity component resolved by the model.', 'm s-1') ! u is ufmean(1:n3)
  call vo_mean%init(n3, 'vo', 'Prognostic x-ward velocity component resolved by the model.', 'm s-1') ! v is ufmean(1+n3:2*n3)
  call wo_mean%init(n3, 'wo', 'vertical component of ocean velocity', 'm s-1') 

  call thetao_mean%init(n3, 'thetao', 'sea water potential temperature', 'degC')
  call so_mean%init(n3, 'so', 'sea water salinity', 'psu')
  call opottemptend_mean%init(n3, 'opottemptend', 'tendency of sea water potential temperature expressed as heat content', 'W m-2')

  if(use_passive_tracer) then
    allocate(passive_tracers_mean(num_passive_tracer))
    do i=1, size(passive_tracers_mean)
      write(tracer_suffix, '(i1)') i
      call passive_tracers_mean(i)%init(n3, 'ptr'//tracer_suffix, 'passive tracer '//tracer_suffix, '')
    end do
  end if
  if(use_age_tracer) then
    allocate(age_tracers_mean(num_age_tracer))
    do i=1, size(age_tracers_mean)
      write(tracer_suffix, '(i1)') i
      call age_tracers_mean(i)%init(n3, 'age'//tracer_suffix, 'age tracer '//tracer_suffix, 'Year')
    end do
  end if


  ! ice part
#ifdef use_ice
  call sisnthick_mean%init(n2, 'sisnthick', 'actual thickness of snow (snow volume divided by snow-covered area)', 'm')
  call ice_area_fraction_mean%init(n2, 'sic', 'fraction of grid cell covered by sea ice.', '1.0')
  call sea_ice_thickness_mean%init(n2, 'sithick', 'actual (floe) thickness of sea ice (NOT volume divided by grid area as was done in CMIP5)', 'm')
  call sea_ice_volume_mean%init(n2, 'sivol', 'total volume of sea ice divided by grid-cell area (this used to be called ice thickness in CMIP5)', 'm')
  call sea_ice_x_velocity_mean%init(n2, 'siu', 'x-velocity of ice on native model grid', 'm s-1')
  call sea_ice_y_velocity_mean%init(n2, 'siv', 'y-velocity of ice on native model grid', 'm s-1')
  call sea_ice_speed_mean%init(n2, 'sispeed', 'speed of ice (i.e. mean absolute velocity) to account for back-and-forth movement of the ice', 'm s-1')
  call sea_ice_time_fraction_mean%init(n2, 'sitimefrac', 'fraction of time steps of the averaging period during which sea ice is present (sic > 0) in a grid cell', '1.0')
  call sisnmass_mean%init(n2, 'sisnmass', 'total mass of snow on sea ice divided by grid-cell area', 'kg m-2')
  call sistrxubot_mean%init(n2, 'sistrxubot', 'x-component of ocean stress on sea ice', 'N m-2')
  call sistryubot_mean%init(n2, 'sistryubot', 'y-component of ocean stress on sea ice', 'N m-2')
#endif

#endif

  ! for diagnostics

#ifdef allow_diag

  ! ocean
  if(diag_oce) then
     if(diag_oce_KE) then
        call u2o_mean%init(n3, 'u2o', 'square of x-component of ocean velocity', 'm2 s-2')
        call v2o_mean%init(n3, 'v2o', 'square of y-component of ocean velocity', 'm2 s-2')
        call w2o_mean%init(n3, 'w2o', 'square of vertical component of ocean velocity', 'm2 s-2')
        call wso_mean%init(n3, 'wso', 'salinity times vertical component of ocean velocity', 'm s-1')
        call wto_mean%init(n3, 'wto', 'temperature times vertical component of ocean velocity', 'degC m s-1')
     end if
     if(diag_oce_energy_conv) then
        call rho_mean%init(n3, 'rho', 'insitu density', 'kg/m3')
        call urho_mean%init(n3, 'urho', 'u * insitu density', 'm/s kg /m3')
        call vrho_mean%init(n3, 'vrho', 'v * insitu density', 'm/s kg /m3')
        call uv_mean%init(n3, 'uv', 'u*v', 'm2/s2')
     end if
     if(diag_oce_transp) then
        call uto_mean%init(n3, 'uto', 'temperature times x-component of ocean velocity', 'degC m s-1')
        call vto_mean%init(n3, 'vto', 'temperature times y-component of ocean velocity', 'degC m s-1')
        call uso_mean%init(n3, 'uso', 'salinity times x-component of ocean velocity', 'm s-1')
        call vso_mean%init(n3, 'vso', 'salinity times y-component of ocean velocity', 'm s-1')
     end if
     if(Redi_GM .and. diag_oce_GM_vel) then
        allocate(sgs_u(myDim_elem3d), sgs_v(myDim_elem3d))
     end if
     if(diag_oce_SGS_transp) then
        allocate(sgs_ut(myDim_elem3d), sgs_vt(myDim_elem3d))
        allocate(sgs_us(myDim_elem3d), sgs_vs(myDim_elem3d))
     end if
     if(diag_oce_mix_layer) then
        call mlotst_mean%init(n2, 'mlotst', 'mixed layer depth, computed with Levitus method (with 0.125 kg/m3 criterion)', 'm')
     end if
  end if

#ifdef use_ice
  ! ice
  if(diag_ice) then
     call thdgr_mean%init(n2, 'thdgr', 'thermodynamic growth rate of eff. ice thickness', 'm/s')
     call thdgrsn_mean%init(n2, 'thdgrsn', 'melting rate of snow thickness', 'm/s')
     call uhice_mean%init(n2, 'uhice', 'zonal advective flux of eff. ice thickness', 'm.m/s')
     call vhice_mean%init(n2, 'vhice', 'meridional advective flux of eff. ice thickness', 'm.m/s')
     call uhsnow_mean%init(n2, 'uhsnow', 'zonal advective flux of eff. snow thickness', 'm.m/s')
     call vhsnow_mean%init(n2, 'vhsnow', 'meridional advective flux of eff. snow thickness', 'm.m/s')
     call flice_mean%init(n2, 'flice', 'rate of flooding snow to ice', 'm/s')
     call fsitherm_mean%init(n2, 'fsitherm', 'sea ice thermodynamic water flux into the ocean divided by the area of the ocean portion of the grid cell', 'kg m-2 s-1')
     call sisnconc_mean%init(n2, 'sisnconc', 'fraction of sea ice, by area, which is covered by snow, giving equal weight to every square metre of sea ice', '1')
     call sidmassth_mean%init(n2, 'sidmassth', 'Total change in sea-ice mass from thermodynamic processes divided by grid-cell area', 'kg m-2 s-1')
     call sidmasssi_mean%init(n2, 'sidmasssi', 'The rate of change of sea ice mass due to transformation of snow to sea ice divided by grid-cell area', 'kg m-2 s-1')
     call sidmasstranx_mean%init(n2, 'sidmasstranx', 'Includes transport of both sea ice and snow by advection', 'kg s-1 m-1')
     call sidmasstrany_mean%init(n2, 'sidmasstrany', 'Includes transport of both sea ice and snow by advection', 'kg s-1 m-1')
  end if

  ! forcing
  if(diag_forcing) then
#ifndef __oasis
     call tair_mean%init(n2, 'tair', 'air temperature', 'degC')
     call shum_mean%init(n2, 'shum', 'air specific humidity', 'kg/kg')
     call lwrd_mean%init(n2, 'lwrd', 'atmosphere longwave radiation', 'W/m^2')
     call olat_mean%init(n2, 'olat', 'latent heat flux to ocean, downward positive', 'W/m^2')
     call olwout_mean%init(n2, 'olwout', 'longwave radiation from ocean, downward positve', 'W/m^2')
     call osen_mean%init(n2, 'osen', 'sensible heat flux to ocean, downward positive', 'W/m^2')
     call relax_salt_mean%init(n2, 'relax_salt', 'ocean surface salinity relaxation, >0 increase salinity', 'psu m/s')
#endif
     call uwind_mean%init(n2, 'uwind', 'Dummy: FESOM does not see it; zonal wind speed', 'm/s')
     call vwind_mean%init(n2, 'vwind', 'Dummy: FESOM does not see it ;meridional wind speed', 'm/s')
     call prlq_mean%init(n2, 'prlq', 'computed as the total mass of liquid water falling as liquid rain  into the ice-free portion of the ocean divided by the area of the ocean portion of the grid cell.', 'kg m-2 s-1')
     call prsn_mean%init(n2, 'prsn', 'at surface; includes precipitation of all forms of water in the solid phase', 'kg m-2 s-1')
     call runoff_mean%init(n2, 'runoff', 'runoff', 'm/s')
     call evap_mean%init(n2, 'evap', 'evaporation', 'm/s')
     call rsdo_mean%init(n2, 'rsdo', 'n/a', 'W/m^2')
     call hfds_mean%init(n2, 'hfds', 'This is the net flux of heat entering the liquid water column through its upper surface (excluding any "flux adjustment") .', 'W m-2')
     call wnet_mean%init(n2, 'wnet', 'net freshwater flux to ocean, downward positive', 'm/s')
     call virtual_salt_mean%init(n2, 'virtual_salt', 'virtual salt flux to ocean, >0 increase salinity', 'psu m/s')
     call tauuo_mean%init(n2, 'tauuo', 'This is the stress on the liquid ocean from overlying atmosphere, sea ice, ice shelf, etc.', 'N m-2')
     call tauvo_mean%init(n2, 'tauvo', 'This is the stress on the liquid ocean from overlying atmosphere, sea ice, ice shelf, etc.', 'N m-2')
  end if
#endif

#endif

  call clean_meanarrays

  if(mype==0) write(*,*) 'Mean arrays have been set up'

  return
end subroutine init_meanarrays
!=============================================================================!


!=============================================================================!
subroutine add2meanarrays
  ! adds values to the mean-arrays
  use o_mesh
  use o_param
  use o_array
  use i_array
  use g_config
  use g_diag
  use g_meanarrays
  use g_forcing_arrays
  use g_parfe
  use i_therm_parms, only : rhoice, rhosno
  use o_mixing_kpp_mod, only : hbl
  use cmor_variables_diag
  use g_rotate_grid
  use i_therm_parms, only : rhosno, rhowat
  implicit none
  !
  integer       :: m, row, row2, row3, j, i, tracer_index_offset
  real(kind=8) u_geo, v_geo, u_wind_geo, v_wind_geo, u_ice_geo, v_ice_geo
  real(kind=8) stress_x_geo, stress_y_geo, stress_iceoce_x_geo, stress_iceoce_y_geo, stress_atmice_x_geo, stress_atmice_y_geo


  meancounter=meancounter+1
  do i=1, size(all_meandata)
    all_meandata(i)%ptr%addcounter = all_meandata(i)%ptr%addcounter + 1
  end do

#ifdef allow_calcmeans
#ifndef use_non_hydrostatic
  call vvel_nodes 
#endif
#endif

  volo_const%local_values(1) = volo_const%local_values(1) + volo
  soga_mean%local_values(1) = soga_mean%local_values(1) + soga
  thetaoga_mean%local_values(1) = thetaoga_mean%local_values(1) + thetaoga
  siarean_mean%local_values = siarean_mean%local_values + siarean
  siareas_mean%local_values = siareas_mean%local_values + siareas
  siextentn_mean%local_values = siextentn_mean%local_values + siextentn
  siextents_mean%local_values = siextents_mean%local_values + siextents
  sivoln_mean%local_values = sivoln_mean%local_values + sivoln
  sivols_mean%local_values = sivols_mean%local_values + sivols
  
  do row=1,myDim_nod3d                      

     row2=row+myDim_nod3d+eDim_nod3D
     u_geo=uf(row)
     v_geo=uf(row2)
     call vector_r2g(u_geo, v_geo, coord_nod3D(1, row), coord_nod3D(2, row), 0)

#ifdef allow_calcmeans

     ! ocean
     uo_mean%local_values(row) = uo_mean%local_values(row) + u_geo
     vo_mean%local_values(row) = vo_mean%local_values(row) + v_geo
#ifdef use_non_hydrostatic
     row3=row2+myDim_nod3D+eDim_nod3D       
     wo_mean%local_values(row) = wo_mean%local_values(row) + uf(row3)
#else
     wo_mean%local_values(row) = wo_mean%local_values(row) + wrhs(row)
#endif
    thetao_mean%local_values(row)       = thetao_mean%local_values(row) + tracer(row,1)
    so_mean%local_values(row)       = so_mean%local_values(row) + tracer(row,2)
    opottemptend_mean%local_values(row) = opottemptend_mean%local_values(row) + opottemptend(row)
    if(use_passive_tracer) then
      tracer_index_offset = 2 ! 1st is temp, 2nd is salt
      do i=1, size(passive_tracers_mean)
        passive_tracers_mean(i)%local_values(row) = passive_tracers_mean(i)%local_values(row) + tracer(row,i+tracer_index_offset)
      end do
    end if
    if(use_age_tracer) then
      tracer_index_offset = 2 + size(passive_tracers_mean) ! 1st is temp, 2nd is salt, then all passive tracers
      do i=1, size(age_tracers_mean)
        age_tracers_mean(i)%local_values(row) = age_tracers_mean(i)%local_values(row) + tracer(row,i+tracer_index_offset)
      end do
    end if
#endif

#ifdef allow_diag
     
     if(diag_oce) then
        if(diag_oce_KE) then
           u2o_mean%local_values(row) = u2o_mean%local_values(row) + u_geo*u_geo
           v2o_mean%local_values(row) = v2o_mean%local_values(row) + v_geo*v_geo
           w2o_mean%local_values(row) = w2o_mean%local_values(row) + wrhs(row)*wrhs(row)
           wso_mean%local_values(row) = wso_mean%local_values(row) + wrhs(row)*tracer(row,2)
           wto_mean%local_values(row) = wto_mean%local_values(row) + wrhs(row)*tracer(row,1)
        end if

        if(diag_oce_energy_conv) then
           rho_mean%local_values(row) = rho_mean%local_values(row) + density_insitu(row)
           urho_mean%local_values(row) = urho_mean%local_values(row) + u_geo*density_insitu(row)
           vrho_mean%local_values(row) = vrho_mean%local_values(row) + v_geo*density_insitu(row)
           uv_mean%local_values(row) = uv_mean%local_values(row) + u_geo*v_geo
        end if
        
        if(diag_oce_transp) then
        uto_mean%local_values(row) = uto_mean%local_values(row) + u_geo*tracer(row,1)
        vto_mean%local_values(row) = vto_mean%local_values(row) + v_geo*tracer(row,1)
        uso_mean%local_values(row) = uso_mean%local_values(row) + u_geo*tracer(row,2)
        vso_mean%local_values(row) = vso_mean%local_values(row) + v_geo*tracer(row,2)
        
     endif
     end if
#endif
  end do

  !-----------------------------------------------------------

  do row=1,myDim_nod2d  

     ! ocean
#ifdef allow_calcmeans 
     zos_mean%local_values(row)       = zos_mean%local_values(row) + ssh(row)
     zossq_mean%local_values(row)       = zossq_mean%local_values(row) + ssh(row)*ssh(row)
     omldamax_mean%local_values(row) = omldamax_mean%local_values(row) + hbl(row)
     pbo_mean%local_values(row) = pbo_mean%local_values(row) + pbo(row)
     tos_mean%local_values(row) = tos_mean%local_values(row) + tos(row)
     tso_mean%local_values(row) = tos(row)+273.15 ! to Kelvin, we only need time-point data here, i.e. snapshot, this could be optimized to only store the data if the write frequency is due
     sos_mean%local_values(row) = sos_mean%local_values(row) + sos(row)
#if defined (__oasis)
     evs_mean%local_values(row) = evs_mean%local_values(row) + evap_no_ifrac(row)*(1.0-a_ice(row))*1000.0
     sidmassevapsubl_mean%local_values(row) = sidmassevapsubl_mean%local_values(row) + sublimation(row)*a_ice(row)
     sifllatstop_mean%local_values(row) = sifllatstop_mean%local_values(row) + ice_heat_flux(row)
#endif
     if(a_ice(row)>0.001) then
     stress_atmice_x_geo=stress_atmice_x(row)
     stress_atmice_y_geo=stress_atmice_y(row)
     call vector_r2g(stress_atmice_x_geo, stress_atmice_y_geo, coord_nod2D(1, row), coord_nod2D(2, row), 0)
     sistrxdtop_mean%local_values(row) = sistrxdtop_mean%local_values(row) + stress_atmice_x_geo
     sistrydtop_mean%local_values(row) = sistrydtop_mean%local_values(row) + stress_atmice_y_geo
     end if
     wfo_mean%local_values(row) = wfo_mean%local_values(row) + water_flux(row)*1000.0
#endif

     ! ice
#ifdef use_ice
#ifdef allow_calcmeans
     u_ice_geo=u_ice(row)
     v_ice_geo=v_ice(row)     
     call vector_r2g(u_ice_geo, v_ice_geo, coord_nod2D(1, row), coord_nod2D(2, row), 0)
     ice_area_fraction_mean%local_values(row) = ice_area_fraction_mean%local_values(row) + a_ice(row)
     sea_ice_thickness_mean%local_values(row) = sea_ice_thickness_mean%local_values(row) + m_ice(row)/(a_ice(row)+1.0e-3)
     sea_ice_volume_mean%local_values(row) = sea_ice_volume_mean%local_values(row) + m_ice(row)
     sisnthick_mean%local_values(row) = sisnthick_mean%local_values(row) + m_snow(row)/(a_ice(row)+1.0e-3)
     sea_ice_x_velocity_mean%local_values(row) = sea_ice_x_velocity_mean%local_values(row) + u_ice_geo
     sea_ice_y_velocity_mean%local_values(row) = sea_ice_y_velocity_mean%local_values(row) + v_ice_geo
     sea_ice_speed_mean%local_values(row) = sea_ice_speed_mean%local_values(row) + sqrt(u_ice_geo*u_ice_geo+v_ice_geo*v_ice_geo)
     if(a_ice(row) > 0.001) then
       sea_ice_time_fraction_mean%local_values(row) = sea_ice_time_fraction_mean%local_values(row) + 1
     endif
     sisnmass_mean%local_values(row) =  sisnmass_mean%local_values(row) + m_snow(row)*rhosno
     if(a_ice(row)>0.001) then
     stress_iceoce_x_geo=stress_iceoce_x(row)
     stress_iceoce_y_geo=stress_iceoce_y(row)
     call vector_r2g(stress_iceoce_x_geo, stress_iceoce_y_geo, coord_nod2D(1, row), coord_nod2D(2, row), 0)
     sistrxubot_mean%local_values(row) = sistrxubot_mean%local_values(row) + stress_iceoce_x_geo
     sistryubot_mean%local_values(row) = sistryubot_mean%local_values(row) + stress_iceoce_y_geo
     end if
#endif
#ifdef allow_diag
     if(diag_ice) then
        thdgr_mean%local_values(row) = thdgr_mean%local_values(row) + thdgr(row) !warning, thdgr is the same as fsitherm, sidmassth/rhoice
        thdgrsn_mean%local_values(row)  = thdgrsn_mean%local_values(row) + thdgrsn(row)
        uhice_mean%local_values(row)    = uhice_mean%local_values(row) + u_ice_geo*m_ice(row)
        vhice_mean%local_values(row)    = vhice_mean%local_values(row) + v_ice_geo*m_ice(row)
        uhsnow_mean%local_values(row)   = uhsnow_mean%local_values(row) + u_ice_geo*m_snow(row)
        vhsnow_mean%local_values(row)   = vhsnow_mean%local_values(row) + v_ice_geo*m_snow(row)
        flice_mean%local_values(row)    = flice_mean%local_values(row) + flice(row)
        fsitherm_mean%local_values(row) = fsitherm_mean%local_values(row)+thdgr(row)*rhoice !warning, thdgr is the same as fsitherm, sidmassth/rhoice
        if (m_snow(row)>1.e-3) then
           sisnconc_mean%local_values(row) = sisnconc_mean%local_values(row)+a_ice(row)
        endif
        sidmassth_mean%local_values(row) = sidmassth_mean%local_values(row) + thdgr(row)*rhoice!warning, thdgr is the same as fsitherm, sidmassth/rhoice
        sidmasssi_mean%local_values(row) = sidmasssi_mean%local_values(row) + flice(row)*rhoice
        sidmasstranx_mean%local_values(row) = sidmasstranx_mean%local_values(row) + u_ice_geo*(rhoice*m_ice(row)+rhosno*m_snow(row))
        sidmasstrany_mean%local_values(row) = sidmasstrany_mean%local_values(row) + v_ice_geo*(rhoice*m_ice(row)+rhosno*m_snow(row))
     endif
#endif
#endif

     ! forcing
#ifdef allow_diag
#ifdef use_ice
     if(diag_forcing) then
#ifndef __oasis
        tair_mean%local_values(row)     = tair_mean%local_values(row) + Tair(row)
        shum_mean%local_values(row)     = shum_mean%local_values(row) + shum(row)
        lwrd_mean%local_values(row)     = lwrd_mean%local_values(row) + longwave(row)
        olat_mean%local_values(row)     = olat_mean%local_values(row) + olat_heat(row)
        osen_mean%local_values(row)     = osen_mean%local_values(row) + osen_heat(row)
        olwout_mean%local_values(row)   = olwout_mean%local_values(row) + olwout(row)
        relax_salt_mean%local_values(row)   = relax_salt_mean%local_values(row) + relax_salt(row)
#endif
        u_wind_geo=u_wind(row)
        v_wind_geo=v_wind(row)
        call vector_r2g(u_wind_geo, v_wind_geo, coord_nod2D(1, row), coord_nod2D(2, row), 0)
        uwind_mean%local_values(row)    = uwind_mean%local_values(row) + u_wind_geo
        vwind_mean%local_values(row)    = vwind_mean%local_values(row) + v_wind_geo  
        prlq_mean%local_values(row)     = prlq_mean%local_values(row) + prec_rain(row)*rhowat
        prsn_mean%local_values(row)     = prsn_mean%local_values(row) + prec_snow(row)*rhosno
        runoff_mean%local_values(row)   = runoff_mean%local_values(row) + runoff(row)
        evap_mean%local_values(row)     = evap_mean%local_values(row) + evaporation(row)
        rsdo_mean%local_values(row)     = rsdo_mean%local_values(row) + shortwave(row)
        hfds_mean%local_values(row)     = hfds_mean%local_values(row) + net_heat_flux(row)
        wnet_mean%local_values(row)     = wnet_mean%local_values(row) + fresh_wa_flux(row)
        virtual_salt_mean%local_values(row) = virtual_salt_mean%local_values(row) + virtual_salt(row)
        stress_x_geo=stress_x(row)
        stress_y_geo=stress_y(row)
        call vector_r2g(stress_x_geo, stress_y_geo, coord_nod2D(1, row), coord_nod2D(2, row), 0)
        tauuo_mean%local_values(row) = tauuo_mean%local_values(row) + stress_x_geo
        tauvo_mean%local_values(row) = tauvo_mean%local_values(row) + stress_y_geo
     endif
#endif
#endif

  enddo

  ! other diagnostics

  ! mixed layer thickness
#ifdef allow_diag
  if(diag_oce .and. diag_oce_mix_layer) then
     call compute_mixlay
  end if
#endif

  ! diagnose for SGS parameterization is done in the ts_rhs routine

  return
end subroutine add2meanarrays
!=============================================================================!


!=============================================================================!
subroutine compute_means
  !computes the mean values for output
  use o_param
  use g_config
  use g_diag
  use g_meanarrays
  implicit none

  real(kind=8)   :: cnt

  cnt=float(max(meancounter,1))

  ! diagnose
#ifdef allow_diag

  ! ocean
  if(diag_oce) then
     if(Redi_GM .and. diag_oce_GM_vel) then
        sgs_u                 = sgs_u              /cnt
        sgs_v                 = sgs_v              /cnt
     endif
     if(diag_oce_SGS_transp) then
        sgs_ut                = sgs_ut             /cnt
        sgs_vt                = sgs_vt             /cnt
        sgs_us                = sgs_us             /cnt
        sgs_vs                = sgs_vs             /cnt
     endif
  endif

#endif

end subroutine compute_means
!=============================================================================!


!=============================================================================!
subroutine clean_meanarrays
  ! puts zeros into the mean-arrays
  use o_param
  use g_config
  use g_diag
  use g_meanarrays
  implicit none

  meancounter=0

  ! diagnose
#ifdef allow_diag

  if(diag_oce) then
     if(Redi_GM .and. diag_oce_GM_vel) then
        sgs_u=0.
        sgs_v=0.
     endif
     if(diag_oce_SGS_transp) then
        sgs_ut=0.
        sgs_vt=0.
        sgs_us=0.
        sgs_vs=0.
     endif
  endif

#endif

  return
end subroutine clean_meanarrays
!=============================================================================!
