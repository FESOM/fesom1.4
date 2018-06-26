subroutine cal_shortwave_rad
  ! Calculate shortwave radiation in the ocean according to chlorophyll climatology
  ! Assuming under ice no penetration. A decent way for ice region is to be discussed.
  ! This routine should be called after ice2oce coupling done if ice model is used.
  ! Ref.: Morel and Antoine 1994, Sweeney et al. 2005
  use o_mesh
  use o_param
  use o_array
  use g_forcing_arrays
  use g_parfe
  use g_config
  use i_therm_parms
  use i_array

  implicit none

  integer      :: m, n2, n3, k
  real(kind=8) :: swsurf, aux, z
  real(kind=8) :: c, c2, c3, c4, c5
  real(kind=8) :: v1, v2, sc1, sc2

  sw_3d=0.0

  do n2=1,myDim_nod2D+eDim_nod2D     
#ifdef use_cavity
     if(cavity_flag_nod2d(n2)==1) cycle   
#endif
#ifdef use_ice  
     if(a_ice(n2)>0.0) cycle !assume in ice region no penetration
#endif     

     ! shortwave rad.

#if defined (__oasis) || defined (__uncplecham6)     
     swsurf=shortwave(n2)
#else
     swsurf=(1.0-albw)*shortwave(n2)
#endif     
     ! the visible part (300nm-750nm)
     swsurf=swsurf*0.54
     ! subtract visible sw rad. from heat_flux, which is '+' for upward
     heat_flux(n2)=heat_flux(n2)+swsurf 

     ! attenuation func. for vis. sw rad. according to Morel/Antoine param.
     ! the four parameters in the func.

     ! limit chl from below
     if(chl(n2)<0.02) chl(n2)=0.02

     c=log10(chl(n2))  
     c2=c*c
     c3=c2*c
     c4=c3*c
     c5=c4*c
     v1=0.008*c+0.132*c2+0.038*c3-0.017*c4-0.007*c5
     v2=0.679-v1
     v1=v1+0.321
     sc1=1.54-0.197*c+0.166*c2-0.252*c3-0.055*c4+0.042*c5
     sc2=7.925-6.644*c+3.662*c2-1.815*c3-0.218*c4+0.502*c5

     ! convert from heat flux [W/m2] to temperature flux [K m/s]
     swsurf=swsurf/vcpw

     ! vis. sw. rad. in the colume
     sw_3d(nod3d_below_nod2d(1,n2))=swsurf
     do k=2,num_layers_below_nod2d(n2)+1
        n3=nod3d_below_nod2d(k,n2)
        z=coord_nod3d(3,n3)
        aux=(v1*exp(z/sc1)+v2*exp(z/sc2))
        sw_3d(n3)=swsurf*aux
        if(aux<0.00001) exit
     end do
  end do

end subroutine cal_shortwave_rad
