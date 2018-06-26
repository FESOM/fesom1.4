subroutine ice_init
  !sets inital values or reads restart file for ice model
  use i_ARRAY
  use o_MESH    
  use o_PARAM   
  use o_array        
  use g_parfe 
  use g_config
  implicit none
  !
  integer        :: i

  if (.not.r_restart) then
     m_ice =0.
     a_ice =0.
     u_ice =0.
     v_ice =0.
     m_snow=0.

     if(use_prepared_init_ice) then
        if(mype==0) write(*,*) 'initialize the sea ice with previous simulation results'
        call read_prepared_initial_ice
     else           
        if(mype==0) write(*,*) 'initialize the sea ice according to initial SST'
 	do i=1,ToDim_nod2d          
           if(cavity_flag_nod2d(i)==1) cycle              
           if (tracer(nod3D_below_nod2D(1,i),1) < -1.0) then  
	      if(coord_nod2d(2,i)>0.) then
	         m_ice(i)=2.0
	      else
                 m_ice(i) = 1.0
	      endif   
              a_ice(i) = 0.9
              u_ice(i) = 0.
              v_ice(i) = 0.
              m_snow(i)= 0. 
           endif
        enddo
     end if

  else

     if(mype==0) write(*,*) 'read ice restart file'
     call ice_input     

  endif
end subroutine ice_init
!
!----------------------------------------------------------------------------
!
subroutine ice_array_setup
  !inializing sea ice model 
  use o_param
  use o_mesh
  use o_elements
  use i_array
  use g_forcing_arrays
  use g_rotate_grid
  use g_parfe
  implicit none
  !
  integer       :: i, n2
  real(kind=8)  :: midlat, lon, lat, rlon, rlat

  n2=ToDim_nod2D           

  ! Allocate memory for variables of ice model      
  allocate(u_ice(n2), v_ice(n2))
  allocate(m_ice(n2), a_ice(n2), m_snow(n2))
  allocate(dm_ice(n2), da_ice(n2), dm_snow(n2))
  allocate(sigma11(myDim_elem2D), sigma12(myDim_elem2D), sigma22(myDim_elem2D))
  allocate(rhs_m(n2), rhs_a(n2), rhs_u(n2), rhs_v(n2))
  allocate(rhs_ms(n2))
  allocate(t_skin(n2))
  rhs_m=0.
  rhs_ms=0.
  rhs_a=0.
  rhs_u=0.
  rhs_v=0.
  u_ice=0.
  v_ice=0.
  m_ice=0.
  a_ice=0.
  m_snow=0.
  dm_ice=0.
  da_ice=0.
  dm_snow=0.
  sigma11=0.
  sigma12=0.
  sigma22=0.
  t_skin=0.

  ! Allocate memory used for coupling
  allocate(S_oc_array(n2), T_oc_array(n2))
  allocate(fresh_wa_flux(n2), net_heat_flux(n2))
#if defined (__oasis) || defined (__uncplecham6)
  allocate(oce_heat_flux(n2), ice_heat_flux(n2))
  allocate(tmp_oce_heat_flux(n2), tmp_ice_heat_flux(n2))
#endif /* (__oasis) || (__uncplecham6) */
  allocate(stress_atmice_x(n2), stress_atmice_y(n2))    
  allocate(stress_atmoce_x(n2), stress_atmoce_y(n2))    
  allocate(stress_iceoce_x(n2), stress_iceoce_y(n2)) 
  allocate(u_w(n2), v_w(n2))
  allocate(elevation(n2))

  S_oc_array=0.
  T_oc_array=0.
  fresh_wa_flux=0.
  net_heat_flux=0.
  stress_atmice_x=0.
  stress_atmice_y=0.
  stress_atmoce_x=0.
  stress_atmoce_y=0.
  stress_iceoce_x=0.
  stress_iceoce_y=0.

#ifdef use_ice_fct
  call fct_ice_init
#endif

#if defined (__oasis) || defined (__uncplecham6)
  oce_heat_flux=0.
  ice_heat_flux=0.
  tmp_oce_heat_flux=0.
  tmp_ice_heat_flux=0.
#endif /* (__oasis) || (__uncplecham6) */  

  if(mype==0) write(*,*) 'ice arrays have been set up'   

end subroutine ice_array_setup



