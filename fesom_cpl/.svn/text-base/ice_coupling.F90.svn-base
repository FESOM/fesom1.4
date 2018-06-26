!=======================================================================
subroutine ice2ocean
  !transmits the relevant fields from the ice to the ocean model
  !Ralph Timmermann, 15.9.2004

  use o_PARAM
  use o_array
  use o_mesh
  use i_array
  use g_parfe
  use i_dyn_parms
  use i_therm_parms
  use g_config
  implicit none

  integer          :: row
  real*8           :: aux

  do row=1,myDim_nod2d+eDim_nod2d            
#ifdef use_cavity
     if(cavity_flag_nod2d(row)==1) cycle
#endif   

     !heat flux:
     heat_flux(row) = -net_heat_flux(row)     ! W/m2  !heat_flux >0 f. upward 
     ! net_heat_flux > 0 f. downward

     !fresh water flux:
     water_flux(row) = -fresh_wa_flux(row)    ! m freshwater/s  !water_flux > 0 f. upward
     ! fresh_wa_flux > 0 f. downward

     !momentum flux:
     aux=sqrt((u_ice(row)-u_w(row))**2+(v_ice(row)-v_w(row))**2)*rho0*Cd_oce_ice
     stress_iceoce_x(row) = aux * (u_ice(row)-u_w(row))
     stress_iceoce_y(row) = aux * (v_ice(row)-v_w(row))
     stress_x(row)=stress_iceoce_x(row)*a_ice(row) + stress_atmoce_x(row)*(1.-a_ice(row))
     stress_y(row)=stress_iceoce_y(row)*a_ice(row) + stress_atmoce_y(row)*(1.-a_ice(row))
     if(a_ice(row)<0.001) then
        stress_iceoce_x(row)=0.0
        stress_iceoce_y(row)=0.0
     end if
  end do

end subroutine ice2ocean
!=======================================================================


!=======================================================================
subroutine ocean2ice
  !transmits the relevant fields from the ocean to the ice model
  !Ralph Timmermann, 15.9.2004

  use o_param
  use o_array
  use i_array
  use o_MESH
  use g_parfe
  use g_config
  implicit none

  integer :: m, row

  ! the arrays in the ice model are renamed

  do row=1, myDim_nod2d+eDim_nod2d       
     m=nod3D_below_nod2D(1,row)       
     T_oc_array(row)=tracer(m,1)  
     S_oc_array(row)=tracer(m,2)  
     u_w(row) = uf(m)                        
     v_w(row) = uf(m+myDim_nod3d+eDim_nod3D)  
     elevation(row)= ssh(row)
  enddo
end subroutine ocean2ice
!=======================================================================

