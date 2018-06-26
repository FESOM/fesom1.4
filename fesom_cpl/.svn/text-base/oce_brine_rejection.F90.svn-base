subroutine cal_brine_rejection
  ! Calculate total salt released via brine rejection during ice formation.
  ! Unit: psu m3
  ! This need to be done before SSS is updated through adv-diff equation.
  !
  ! Coded by Qiang Wang, 10.02.2012
  ! Reviewed by ??
  !--------------------------------------------------------
  
  use o_mesh
  use o_elements
  use o_array, only: salt_brine_rejection
  use i_therm_parms, only: rhowat, rhoice, Sice
  use i_array, only: S_oc_array
  use g_config
  use g_forcing_arrays, only: thdgr
  use g_parfe
  implicit none

  integer         :: row
  real(kind=8)    :: aux

  aux=rhoice/rhowat*dt
  
  do row=1,myDim_nod2d   ! myDim is sufficient
     if(thdgr(row)>0.0) then
        salt_brine_rejection(row)= &
             (S_oc_array(row)-Sice)*thdgr(row)*aux*cluster_area_2d(row) 
        !unit: psu m3
     else
        salt_brine_rejection(row)=0.0
     end if
  end do
end subroutine cal_brine_rejection
!
!----------------------------------------------------------------------------
!
subroutine apply_brine_rejection
  ! Apply the subgrid scale brine rejection parameterization.
  ! Distribute salt released through brine rejection inside the mixed layer
  ! accorcing to the specified vertical distribution function.
  !
  ! Currently only support the non-linear free surface option!! 
  ! Call this routine after the advection equation and before the vertical
  ! diffusion equation.
  !
  ! Ref: Duffy1997, Duffy1999, Nguyen2009
  !
  ! Coded by Qiang Wang, 10.02.2012
  ! Reviewed by ??
  !--------------------------------------------------------

  use o_mesh
  use o_elements
  use o_array
  use g_parfe
  implicit none

  integer         :: row, nsurf, k, nod, nup, nlo, kml
  real(kind=8)    :: zsurf, rhosurf, drhodz, spar(100)

  integer         :: n_distr
  real(kind=8)    :: drhodz_cri, rho_cri
  data drhodz_cri /0.01/  !kg/m3/m  !NH   !Nguyen2011
  data n_distr /5/
  data rho_cri /0.4/      !kg/m3    !SH   !Duffy1999

  do row=1,myDim_nod2d   ! myDim is sufficient
     nsurf=nod3D_below_nod2D(1, row)     
     if (salt_brine_rejection(row)<=0.0) cycle
     ! do not parameterize brine rejection in regions with low salinity
     ! 1. it leads to further decrease of SSS
     ! 2. in case of non zero salinity of ice (the well accepted value is 5psu) the SSS might become negative
     if (tracer(nsurf,2) < 10.0) cycle
  if (geolat2d(row)>0.0) then  !NH
        kml=1
        zsurf=coord_nod3d(3, nsurf)
        spar(1)=0.0
        do k=1,num_layers_below_nod2d(row)
           nup=nod3d_below_nod2d(k,row)
           nlo=nod3d_below_nod2d(k+1,row)
           drhodz=abs((density_insitu(nlo)-density_insitu(nup))/ &
                (coord_nod3d(3,nup)-coord_nod3d(3,nlo)))
           if (drhodz>=drhodz_cri .or. coord_nod3d(3,nlo)<-80.0) exit
           kml=kml+1
           spar(k+1)=cluster_vol_3d(nlo)*(zsurf-coord_nod3d(3,nlo))**n_distr
        end do

        if (kml>1) then
           tracer(nsurf,2)=tracer(nsurf,2)-salt_brine_rejection(row)/cluster_vol_3d(nsurf)
           spar(2:kml)=spar(2:kml)/sum(spar(2:kml))
           do k=2,kml
              nod=nod3d_below_nod2d(k,row)
              tracer(nod,2)=tracer(nod,2)+salt_brine_rejection(row)*spar(k)/cluster_vol_3d(nod)
           end do
        endif
!  else  !SH
!        tracer(nsurf,2)=tracer(nsurf,2)-salt_brine_rejection(row)/cluster_vol_3d(nsurf)
!        kml=1
!        rhosurf=density_insitu(nsurf)
!        spar(1)=cluster_vol_3d(nsurf)
!        do k=1,num_layers_below_nod2d(row)
!           nlo=nod3d_below_nod2d(k+1,row)
!           if(density_insitu(nlo)>=rhosurf+rho_cri) exit
!           kml=kml+1
!           spar(k+1)=cluster_vol_3d(nlo)
!        end do
!        spar(1:kml)=spar(1:kml)/sum(spar(1:kml))
!     
!        do k=1,kml
!           nod=nod3d_below_nod2d(k,row)
!           tracer(nod,2)=tracer(nod,2)+salt_brine_rejection(row)*spar(k)/cluster_vol_3d(nod)
!        end do
!
  endif

  end do

end subroutine apply_brine_rejection
