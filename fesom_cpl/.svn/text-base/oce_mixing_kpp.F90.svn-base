module o_mixing_kpp_mod
  ! Original numerical algorithm by Bill Large at NCAR, 1994
  ! Equation numbers in the code refer to the Large etal paper. 
  !
  ! Modified from MOM4 to FESOM by Qiang Wang, Feb. 2011
  ! checked by ??
  !---------------------------------------------------------------  

  use o_MESH
  use o_ELEMENTS
  use g_config
  use o_PARAM
  use o_array
  use g_forcing_arrays
  use g_PARFE

  implicit none
  save

  private

  public oce_mixing_kpp_init
  public oce_mixing_kpp

  public hbl, ghats, blmc

  private bldepth
  private wscale
  private ri_iwmix
  private ddmix
  private blmix_kpp
  private enhance
  private cal_nodal_alpha_beta


  real, dimension(:), allocatable      :: bfsfc    ! surface buoyancy forcing    (m^2/s^3)
  real, dimension(:), allocatable      :: caseA    ! = 1 in case A; =0 in case B
  real, dimension(:), allocatable      :: stable   ! = 1 in stable forcing; =0 in unstable
  real, dimension(:,:), allocatable    :: dkm1     ! boundary layer diff_cbt at kbl-1 level
  real, dimension(:,:), allocatable    :: blmc     ! boundary layer mixing coefficients
  real, dimension(:), allocatable      :: Talpha   ! -d(rho)/ d(pot.temperature)  (kg/m^3/C)
  real, dimension(:), allocatable      :: Sbeta    ! d(rho)/ d(salinity)       (kg/m^3/PSU)
  real, dimension(:), allocatable      :: ustar    ! surface friction velocity       (m/s)
  real, dimension(:), allocatable      :: Bo       ! surface turb buoy. forcing  (m^2/s^3)
  real, dimension(:), allocatable      :: dVsq     ! (velocity shear re sfc)^2   (m/s)^2

  integer, dimension(:), allocatable   :: kbl      ! index of first grid level below hbl

  real, dimension(:), allocatable      :: ghats    ! nonlocal transport (s/m^2)
  real, dimension(:), allocatable      :: hbl      ! boundary layer depth

  real, parameter :: epsln   = 1.0e-40_08 ! a small value

  real, parameter :: epsilon = 0.1
  real, parameter :: vonk    = 0.4
  real, parameter :: conc1   = 5.0

  real, parameter :: zmin    = -4.e-7  ! m3/s3 limit for lookup table of wm and ws
  real, parameter :: zmax    = 0.0     ! m3/s3 limit for lookup table of wm and ws
  real, parameter :: umin    = 0.0     ! m/s limit for lookup table of wm and ws
  real, parameter :: umax    = 0.04    ! m/s limit for lookup table of wm and ws
  logical         :: limit_forcing=.false.

  real :: Ricr   = 0.3  ! critical bulk Richardson Number
  real :: concv  = 1.6  ! constant for pure convection (eqn. 23) (Large 1.5-1.6; MOM default 1.8)

  real :: cg            ! non-dimensional coefficient for counter-gradient term
  real :: Vtc           ! non-dimensional coefficient for velocity 
  ! scale of turbulant velocity shear        
  ! (=function of concv,concs,epsilon,vonk,Ricr)

  real :: deltaz        ! delta zehat in table
  real :: deltau        ! delta ustar in table

  integer, parameter :: nni = 890         ! number of values for zehat in the look up table
  integer, parameter :: nnj = 480         ! number of values for ustar in the look up table
  real, dimension(0:nni+1,0:nnj+1) :: wmt ! lookup table for wm, the turbulent velocity scale for momentum
  real, dimension(0:nni+1,0:nnj+1) :: wst ! lookup table for ws, the turbulent velocity scale scalars

  ! namelist parameters in the future: 
  ! concv, Ricr

contains


  !#######################################################################
  ! <SUBROUTINE NAME="oce_mixing_kpp_init">
  !
  ! Initialization for the KPP vertical mixing scheme
  !
  !     output:
  !       Vtc = non-dimensional constant used in calc. bulk Ri              
  !       cg  = constant used in calc.nonlocal transport term                
  !       wmt = turbulent velocity scale for momentum                         
  !       wst = turbulent velocity scale for scaler                          

  subroutine oce_mixing_kpp_init

    implicit none

    real, parameter :: cstar  = 10.0   ! proportionality coefficient for nonlocal transport
    real, parameter :: conam  = 1.257
    real, parameter :: concm  = 8.380 
    real, parameter :: conc2  = 16.0
    real, parameter :: zetam  = -0.2
    real, parameter :: conas  = -28.86
    real, parameter :: concs  = 98.96
    real, parameter :: conc3  = 16.0
    real, parameter :: zetas  = -1.0

    real :: zehat ! = zeta * ustar**3
    real :: zeta  ! = stability parameter d/L
    real :: usta

    integer :: i, j

    allocate (ghats(ToDim_nod3d))
    allocate (hbl(ToDim_nod2d))

    allocate (bfsfc(ToDim_nod2d))         ! surface buoyancy forcing    (m^2/s^3)
    allocate (caseA(ToDim_nod2d))         ! = 1 in case A; =0 in case B
    allocate (stable(ToDim_nod2d))        ! = 1 in stable forcing; =0 in unstable
    allocate (dkm1(ToDim_nod2d,3))        ! boundary layer diff at kbl-1 level
    allocate (blmc(ToDim_nod3d,3))        ! boundary layer mixing coefficients
    allocate (talpha(Todim_nod3d))        ! -d(rho)/ d(pot.temperature)  (g/m^3/C)
    allocate (sbeta(ToDim_nod3d))         ! d(rho)/ d(salinity)       (g/m^3/PSU)
    allocate (ustar(ToDim_nod2d))         ! surface friction velocity       (m/s)
    allocate (Bo(ToDim_nod2d))            ! surface turb buoy. forcing  (m^2/s^3)
    allocate (dVsq(ToDim_nod3d))          ! (velocity shear re sfc)^2   (m/s)^2
    allocate (kbl(ToDim_nod2d))           ! index of first grid level below hbl

    ghats       = 0.0
    hbl         = 0.0

    !-----------------------------------------------------------------------
    !     initialize some constants for kmix subroutines, and initialize
    !     for kmix subroutine "wscale" the 2D-lookup table for wm and ws
    !     as functions of ustar and zetahat (=vonk*sigma*hbl*bfsfc).
    !-----------------------------------------------------------------------

    !-----------------------------------------------------------------------
    ! define some non-dimensional constants  (recall epsilon=0.1)
    !-----------------------------------------------------------------------

    !     Vtc used in eqn. 23
    Vtc     = concv * sqrt(0.2/concs/epsilon) / vonk**2 / Ricr

    !     cg = cs in eqn. 20
    cg      = cstar * vonk * (concs * vonk * epsilon)**(1./3.)

    !-----------------------------------------------------------------------
    ! construct the wm and ws lookup tables (eqn. 13 & B1)
    !-----------------------------------------------------------------------

    deltaz = (zmax-zmin)/real(nni+1) 
    deltau = (umax-umin)/real(nnj+1)

    do i=0,nni+1
       zehat = deltaz*(i) + zmin
       do j=0,nnj+1
          usta = deltau*(j) + umin
          zeta = zehat/(usta**3+epsln)

          if(zehat >= 0.) then
             wmt(i,j) = vonk*usta/(1.+conc1*zeta)
             wst(i,j) = wmt(i,j)
          else
             if(zeta > zetam) then
                wmt(i,j) = vonk* usta * (1.-conc2*zeta)**(1./4.)
             else
                wmt(i,j) = vonk* (conam*usta**3-concm*zehat)**(1./3.)
             endif
             if(zeta > zetas) then
                wst(i,j) = vonk* usta * (1.-conc3*zeta)**(1./2.)
             else
                wst(i,j) = vonk* (conas*usta**3-concs*zehat)**(1./3.)
             endif
          endif
       enddo
    enddo

  end subroutine oce_mixing_kpp_init


  !#######################################################################
  ! <SUBROUTINE NAME="oce_mixing_kpp">
  !
  ! This subroutine computes the vertical diffusivity and viscosity according
  ! to the KPP scheme of Large etal.  In brief, the scheme does the 
  ! following:
  !
  ! --Compute interior mixing everywhere:                               
  !   interior mixing gets computed due to constant internal wave 
  !   background activity ("visc_cbu_iw" and "diff_cbt_iw").
  !   Mixing is enhanced in places of static instability (local Ri < 0).
  !   Additionally, mixing can be enhanced by contribution from shear 
  !   instability which is a function of the local Ri.
  !
  ! --Double diffusion:
  !   Interior mixing can be enhanced by double diffusion due to salt
  !   fingering and diffusive convection.
  !
  ! --Boundary layer:
  !
  !   (A) Boundary layer depth:
  !       at every gridpoint the depth of the oceanic boundary layer 
  !       ("hbl") gets computed by evaluating bulk richardson numbers.
  !
  !   (B) Boundary layer mixing:
  !       within the boundary layer, above hbl, vertical mixing is 
  !       determined by turbulent surface fluxes, and interior mixing at
  !       the lower boundary, i.e. at hbl.
  !
  ! outputs
  !
  !  hbl   = boundary layer depth (meters)
  !  ghats = nonlocal transport coefficient (s/m^2)
  !  viscA = viscosity coefficient  (m^2/s) 
  !  diffK = diffusion coefficient (m^2/s) 
  !
  subroutine oce_mixing_kpp(viscA, diffK)

    implicit none

    real, dimension(ToDim_nod3d),   intent(inout) :: viscA !for momentum
    real, dimension(ToDim_nod3d,2), intent(inout) :: diffK !for T and S
    integer                             :: i, k, row, kn
    integer                             :: num, nodlo, nodup, ns, lay, lay_mi
    real                                :: smftu, smftv, aux
    real                                :: dens_up, minmix
    real                                :: usurf, vsurf, tsurf_loc, ssurf_loc
! the variables below will be used to limit the forcing in case it out of bounds
! the limits are: umin =< ustar(i) <=umax ; Bo(i)>=zmin/(vonk*0.1*abs(z))
    real                                :: z, bfsfc_lim
    integer                             :: nlims_bf, gnlims_bf

    do i=1,ToDim_nod2d
       kn=num_layers_below_nod2D(i)+1

       !-----------------------------------------------------------------------
       !     compute squared vertical difference of velocity
       !-----------------------------------------------------------------------

       usurf=0.5*(uf(nod3d_below_nod2d(1,i))+uf(nod3d_below_nod2d(2,i)))
       vsurf=0.5*(uf(nod3d_below_nod2d(1,i)+ToDim_nod3d) &
            +uf(nod3d_below_nod2d(2,i)+ToDim_nod3d))
       dVsq(nod3d_below_nod2d(1,i))=0.0
       do k=2, kn
          row=nod3d_below_nod2d(k,i)
          dVsq(row)=(usurf-uf(row))**2+(vsurf-uf(ToDim_nod3d+row))**2
       end do
       
       ! Squared buoyancy frequency and the buoyancy differenc between surface 
       ! and grid points below, which are required in this module,
       ! are calculated in separate routines!!
    end do

    !-----------------------------------------------------------------------
    !     compute thermal and haline expansion coefficients (without factor of rho).
    !     thermal expansion coefficient without 1/rho factor       (kg/m3/C)
    !           talpha= -d(rho{k,k})/d(T(k))
    !     salt expansion coefficient without 1/rho factor        (kg/m3/PSU)
    !           sbeta = d(rho{k,k})/d(S(k))
    !-----------------------------------------------------------------------

    call cal_nodal_alpha_beta

    !-----------------------------------------------------------------------
    ! friction velocity, turbulent sfc buoyancy forcing
    ! ustar = sqrt( sqrt( stress_x^2 + stress_y^2 ) / rho ) (m/s)
    ! bo =  -g * ( Talpha*heat_flux/vcpw + Sbeta * salinity*water_flux ) (m^2/s^3)
    !-----------------------------------------------------------------------

    nlims_bf=0
    do i=1,ToDim_nod2d
       row=nod3d_below_nod2d(1,i)
       ustar(i) = sqrt( sqrt(stress_x(i)**2 + stress_y(i)**2)*rho0r)

       if (limit_forcing) ustar(i)=min(ustar(i), umax)

       Bo(i)    = -g * (Talpha(row) * heat_flux(i)/vcpw  &   !heat_flux & water_flux: positive up
            + Sbeta(row) * water_flux(i) * tracer(row,2))

       if (limit_forcing) then
          kn=num_layers_below_nod2D(i)+1
          row=nod3d_below_nod2d(kn, i)
          z=coord_nod3D(3, row)
          bfsfc_lim=zmin/(vonk*0.1*abs(z))
          if (Bo(i) < bfsfc_lim) then
             Bo(i) = bfsfc_lim
             nlims_bf=nlims_bf+1
             if (i<=MyDim_nod2D) nlims_bf=nlims_bf+1
          end if
       end if
    enddo
    if (limit_forcing) then
       call MPI_AllREDUCE(nlims_bf, gnlims_bf, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_FESOM, MPIerr)
       if (gnlims_bf>0) then
          if (mype==0) write(*,*) 'KPP Warning: buoyancy forcing has been limited from below at ', gnlims_bf, ' points!'
       end if
    end if
    !-----------------------------------------------------------------------
    !     compute interior mixing coefficients everywhere, due to constant 
    !     internal wave activity, static instability, and local shear 
    !     instability.
    !-----------------------------------------------------------------------
    call ri_iwmix(viscA, diffK)

    !-----------------------------------------------------------------------
    !     add double diffusion
    !-----------------------------------------------------------------------

    if (double_diffusion) then
       call ddmix(diffK)
    endif

    !-----------------------------------------------------------------------
    !     boundary layer mixing coefficients: diagnose new b.l. depth
    !-----------------------------------------------------------------------

    call bldepth

    !-----------------------------------------------------------------------
    !     boundary layer diffusivities
    !-----------------------------------------------------------------------

    call blmix_kpp(viscA, diffK)

    !-----------------------------------------------------------------------
    !     enhance diffusivity at interface kbl - 1
    !-----------------------------------------------------------------------

    call enhance(viscA, diffK) 

    !-----------------------------------------------------------------------
    !     combine interior and b.l. coefficients and nonlocal term;
    !     add horizontal smoothing
    !-----------------------------------------------------------------------

    ! first do horizontal smoothing
    ! In structured mesh case: smooth diffusivities with a 5 point filter
    ! to remove 2 delta x noise
    ! In FEOM case: do smoothing twice over the cluster over each grid layer

    if (smooth_blmc) then 
       ! use dVsq and Talpha to save memory
       dVsq=blmc(:,2)
       Talpha=blmc(:,3)
       do i=1, myDim_nod2d  ! only over myDim_nod2d
          num=nghbr_nod2d(i)%nmb
!DS 18.11.2016
!         lay=minval(num_layers_below_nod2d(nghbr_nod2d(i)%addresses))+1
!         lay_mi=maxval(kbl(nghbr_nod2d(i)%addresses))
          lay=max_num_layers
	  lay_mi=1
          do k=1, num
             row=nghbr_nod2d(i)%addresses(k)
             lay=min(lay, num_layers_below_nod2d(row))
             lay_mi=max(lay_mi, kbl(row))
          end do
          lay=lay+1
          lay=min(lay_mi, lay)
!DS 18.11.2016
          do k=1, lay
             aux=sum(blmc(nod3d_below_nod2d(k,nghbr_nod2d(i)%addresses),2))
             aux=aux+blmc(nod3d_below_nod2d(k,i),2)*real(num-2)
             dVsq(nod3d_below_nod2d(k,i))=aux/2./real(num-1)
             aux=sum(blmc(nod3d_below_nod2d(k,nghbr_nod2d(i)%addresses),3))
             aux=aux+blmc(nod3d_below_nod2d(k,i),3)*real(num-2)
             Talpha(nod3d_below_nod2d(k,i))=aux/2./real(num-1)                   
          end do
       end do

       call com_3D(dVsq)
       call com_3D(Talpha)

       blmc(:,2)=dVsq
       blmc(:,3)=Talpha
       
       ! the second time
       do i=1, myDim_nod2d  ! only over myDim_nod2d
          num=nghbr_nod2d(i)%nmb
          lay=minval(num_layers_below_nod2d(nghbr_nod2d(i)%addresses))+1
          lay_mi=maxval(kbl(nghbr_nod2d(i)%addresses))
          lay=min(lay_mi,lay)
          do k=1, lay
             aux=sum(blmc(nod3d_below_nod2d(k,nghbr_nod2d(i)%addresses),2))
             aux=aux+blmc(nod3d_below_nod2d(k,i),2)*real(num-2)
             dVsq(nod3d_below_nod2d(k,i))=aux/2./real(num-1)
             aux=sum(blmc(nod3d_below_nod2d(k,nghbr_nod2d(i)%addresses),3))
             aux=aux+blmc(nod3d_below_nod2d(k,i),3)*real(num-2)
             Talpha(nod3d_below_nod2d(k,i))=aux/2./real(num-1)                   
          end do
       end do
       blmc(:,2)=dVsq
       blmc(:,3)=Talpha
    end if

    ! then combine blmc and viscA/diffK

    do i=1, myDim_nod2d
       do k=1,num_layers_below_nod2d(i)+1
          if (k < kbl(i)) then
             !viscA(nod3d_below_nod2d(k,i))   = &
             !     max(Av0, blmc(nod3d_below_nod2d(k,i),1))
             !diffK(nod3d_below_nod2d(k,i),1) = &
             !     max(Kv0, blmc(nod3d_below_nod2d(k,i),2))
             !diffK(nod3d_below_nod2d(k,i),2) = &
             !     max(Kv0, blmc(nod3d_below_nod2d(k,i),3))
             viscA(nod3d_below_nod2d(k,i))   = &
                  max(viscA(nod3d_below_nod2d(k,i)), blmc(nod3d_below_nod2d(k,i),1))
             diffK(nod3d_below_nod2d(k,i),1) = &
                  max(diffK(nod3d_below_nod2d(k,i),1), blmc(nod3d_below_nod2d(k,i),2))
             diffK(nod3d_below_nod2d(k,i),2) = &
                  max(diffK(nod3d_below_nod2d(k,i),2), blmc(nod3d_below_nod2d(k,i),3))
		  		  
          else
             ghats(nod3d_below_nod2d(k,i))=0.0
          endif
       enddo
       !do k=1, min(3,num_layers_below_nod2d(i)+1)
       !	  viscA(nod3d_below_nod2d(k,i)) = &
       !	max(minmix,viscA(nod3d_below_nod2d(k,i)))
       !enddo
    enddo
    

    ! Set the mixing coeff. in the first layer above some limiting value
    ! this is very helpful to avoid huge surface velocity when vertical
    ! viscosity is very small derived from the KPP scheme.
    ! I strongly recommend this trick, at least in the current FESOM version.	
    minmix=3.0e-3
    where(viscA(1:myDim_nod2D)<minmix) 
       viscA(1:myDim_nod2d)=minmix
    end where
    !where(diffK(1:myDim_nod2D,1)<minmix) 
    !   diffK(1:myDim_nod2d,1)=minmix 
    !end where
    !-----------------------------------------------------------------------
    
    ! non-local contribution will be added to tracer_rhs directly

  end subroutine oce_mixing_kpp


  !#######################################################################
  ! <SUBROUTINE NAME="bldepth">
  !
  ! The oceanic planetray boundary layer depth, hbl, is determined as
  ! the shallowest depth where the bulk richardson number is
  ! equal to the critical value, Ricr.
  !
  ! Bulk Richardson numbers are evaluated by computing velocity and
  ! buoyancy differences between values at level k and surface
  ! reference values.
  !
  ! In this configuration, the reference values are equal to the
  ! mean values between the first two levels.
  ! When using a very fine vertical grid, these values should be 
  ! computed as the vertical average of velocity and buoyancy from 
  ! the surface down to epsilon*zcoord(k).
  !
  ! When the bulk richardson number at k exceeds Ricr, hbl is
  ! linearly interpolated between grid levels zcoord(k) and zcoord(k-1).
  !
  ! The water column and the surface forcing are diagnosed for 
  ! stable/ustable forcing conditions, and where hbl is relative 
  ! to grid points (caseA), so that conditional branches can be 
  ! avoided in later subroutines.
  !
  !
  !  input         
  !      real ustar(t2d)      = surface friction velocity     (m/s)      
  !      real dVsq(t3d)       = (velocity shear re sfc)^2   (m/s)^2
  !      real Bo(t2d)         = surface turbulent buoyancy forcing(m^2/s^3)  
  !      real sw_3d(t3d)      = radiative buoyancy forcing (C m/s)          
  !      real f_c(t2d)        = Coriolis parameter            (1/s)            
  !
  !  output
  !      real hbl(t2d)        ! boundary layer depth              (m)      
  !      real bfsfc(t2d)      ! Bo+radiation absorbed (m^2/s^3)      
  !      real stable(t2d)     ! =1 in stable forcing; =0 unstable          
  !      real caseA(t2d)      ! =1 in case A, =0 in case B                
  !      integer kbl(t2d)     ! index of first grid level below hbl        
  !
  subroutine bldepth
    implicit none

    real            :: Ritop, bvsq, Vtsq
    real            :: hekman, hmonob, hlimit
    real            :: Rib_km1, Rib_k, coeff_sw, zk, zkm1
    real            :: sigma,zehat, wm, ws, dzup, dzlo
    integer         :: i, k, nk, ns, nod, nodup, nodlo

    real, parameter :: cekman = 0.7  ! constant for Ekman depth
    real, parameter :: cmonob = 1.0  ! constant for Monin-Obukhov depth

    ! initialize hbl and kbl to bottomed out values
    kbl=num_layers_below_nod2d(:)+1
    do i=1,ToDim_nod2d
       nod=bt_nds(i)
       hbl(i)=abs(coord_nod3d(3, nod))
    end do

    do i=1,ToDim_nod2d
#ifdef use_sw_pene      
       ns=nod3d_below_nod2d(1,i)
       coeff_sw=g*Talpha(ns)
#endif
       Rib_km1=0.0
       nk=num_layers_below_nod2d(i)+1
       bfsfc(i)=Bo(i) 

       do k=2,nk
          nod=nod3d_below_nod2d(k,i)
          nodup=nod3d_below_nod2d(k-1,i)
          zk=-coord_nod3d(3,nod)     
          zkm1=-coord_nod3d(3,nodup)

          ! bfsfc = Bo + sw contricution
#ifdef use_sw_pene
          bfsfc(i)=Bo(i) + coeff_sw*(sw_3d(ns)-sw_3d(nod))  ! sw_3d: [C m/s], positive downward
#endif
          stable(i) = 0.5 + sign( 0.5, bfsfc(i))
          sigma  = stable(i) + (1.0-stable(i)) * epsilon

          !-----------------------------------------------------------------------
          ! compute velocity scales at sigma, for z=abs(coord_nod3d(3,nod)):
          !-----------------------------------------------------------------------

          zehat= vonk * sigma * zk * bfsfc(i)
          call wscale(zehat, ustar(i), wm, ws)

          !-----------------------------------------------------------------------
          ! compute the turbulent shear contribution to Rib
          ! eqn. (23)
          !-----------------------------------------------------------------------        

          if(k<nk) then
             bvsq=(bfsq_3D(nodup)+bfsq_3D(nod))*0.5
          else
             bvsq = bfsq_3D(nodup)
          end if
          Vtsq = zk * ws * sqrt(abs(bvsq)) * Vtc

          !-----------------------------------------------------------------------
          !  compute bulk Richardson number at new level
          !  eqn. (21)
          !-----------------------------------------------------------------------

          Ritop = zk*dbsfc_3D(nod)  !dbsfc_3D m/s2
          Rib_k = Ritop / (dVsq(nod) + Vtsq + epsln )
          dzup=zk-zkm1

          if(Rib_k > Ricr) then
             ! linearly interpolate to find hbl where Rib = Ricr
             hbl(i) = zkm1 + dzup*(Ricr-Rib_km1)/(Rib_k-Rib_km1+epsln)
             kbl(i) = k
             exit
          else
             Rib_km1=Rib_k
          endif
       end do

       !-----------------------------------------------------------------------
       ! find stability and buoyancy forcing for boundary layer
       !-----------------------------------------------------------------------
#ifdef use_sw_pene
       ! Linear interpolation of sw_3d to depth of hbl
       bfsfc(i) =Bo(i) +  & 
            coeff_sw*(sw_3d(ns)-(sw_3d(nodup)+(sw_3d(nod)-sw_3d(nodup))*(hbl(i)-zkm1)/dzup))
       stable(i)=0.5 + sign( 0.5, bfsfc(i) )
       bfsfc(i) =bfsfc(i) + stable(i) * epsln  ! ensures bfsfc never=0
#endif

       !-----------------------------------------------------------------------
       !        check hbl limits for hekman or hmonob
       !        eqn. (24)
       !-----------------------------------------------------------------------
       if (bfsfc(i) > 0.0) then
          hekman = cekman * ustar(i) / max(abs(coriolis_param_nod2d(i)), epsln)
          hmonob = cmonob * ustar(i)*ustar(i)*ustar(i)     &
               /vonk / (bfsfc(i)+epsln) 
          hlimit = stable(i) * min(hekman,hmonob)
          hbl(i) = min(hbl(i), hlimit)
          hbl(i) = max(hbl(i), -coord_nod3d(3,nod3d_below_nod2d(2,i))) 
       endif

       !-----------------------------------------------------------------------
       !     find new kbl 
       !-----------------------------------------------------------------------
       kbl(i)=nk
       do k=2,nk
          if(-coord_nod3d(3,nod3d_below_nod2d(k,i)) > hbl(i)) then
             kbl(i) = k
             exit
          endif
       enddo

       !-----------------------------------------------------------------------
       !     find stability and buoyancy forcing for final hbl values
       !-----------------------------------------------------------------------      
#ifdef use_sw_pene
       ! Linear interpolation of sw_3d to depth of hbl
       nod=nod3d_below_nod2d(kbl(i),i)
       nodup=nod3d_below_nod2d(kbl(i)-1,i)
       bfsfc(i) =Bo(i) +  & 
            coeff_sw*(sw_3d(ns)-(sw_3d(nodup)+(sw_3d(nod)-sw_3d(nodup)) &
            *(hbl(i)+coord_nod3d(3,nodup)) &
            /(coord_nod3d(3,nodup)-coord_nod3d(3,nod))))
       stable(i)=0.5 + sign( 0.5, bfsfc(i) )
       bfsfc(i) =bfsfc(i) + stable(i) * epsln 
#endif

       !-----------------------------------------------------------------------
       !     determine caseA and caseB
       !     (if hbl is below (deeper than) the mid point of level kbl
       !     then caseA=0  else caseA=1)
       !-----------------------------------------------------------------------

       nod=nod3d_below_nod2d(kbl(i),i)
       nodup=nod3d_below_nod2d(kbl(i)-1,i)
       dzup=coord_nod3d(3,nodup)-coord_nod3d(3,nod)
       caseA(i)  = 0.5 + sign( 0.5, -coord_nod3d(3,nod)-0.5*dzup-hbl(i))

    enddo  !i

  end subroutine bldepth



  !#######################################################################
  ! <SUBROUTINE NAME="wscale">
  !
  ! Compute turbulent velocity scales.
  ! Use a 2D-lookup table for wm and ws as functions of ustar and
  ! zetahat (=vonk*sigma*hbl*bfsfc=zeta*ustar**3).
  !
  ! Note: the lookup table is only used for unstable conditions
  ! (zehat <= 0), in the stable domain wm (=ws) gets computed
  ! directly.
  !
  !  input
  !      real zehat    zeta*us**3
  !      real us       utar(nod), surface friction velocity    (m/s)    
  !  output                                                             
  !      real wm, ws   turbulent velocity scales 
  !
  subroutine wscale(zehat, us, wm, ws)

    implicit none

    real, intent(in)     :: zehat, us
    real, intent(out)    :: wm, ws
    real                 :: zdiff, udiff, zfrac, ufrac, fzfrac
    real                 :: wam, wbm, was, wbs, u3
    integer              :: iz, izp1, ju, jup1

    !-----------------------------------------------------------------------
    !     use lookup table for zehat < zmax only; otherwise use
    !     stable formulae
    !-----------------------------------------------------------------------

    ! zehat = vonk * sigma * abs(depth) * bfsfc

    if (zehat <= zmax) then
! Here was a KPP bug: max was forgotten and zdiff could be negative. This was giving negative values to zfrac.
! If limit_forcing is set to TRUE, this limiting may not be required. Same for udiff.
!      zdiff  = max(zehat-zmin, 0.)
       zdiff  = zehat-zmin
       iz = int( zdiff/deltaz)
       iz = min( iz , nni )
       iz = max( iz , 0  )
       izp1=iz+1

!      udiff  = min(us-umin, umax-umin)
!      ju = int(udiff/deltau)
!      ju = min( ju , nnj)

       udiff  = us-umin
       ju = int( min(udiff/deltau,float(nnj)))
       ju = max( ju ,   0)
       jup1=ju+1

       zfrac = zdiff/deltaz - real(iz)
       ufrac = udiff/deltau - real(ju)

       fzfrac= 1.-zfrac

       wam   = (fzfrac)  * wmt(iz,jup1) + zfrac*wmt(izp1,jup1)
       wbm   = (fzfrac)  * wmt(iz,ju  ) + zfrac*wmt(izp1,ju  )
       wm    = (1.-ufrac)* wbm          + ufrac*wam

       was   = (fzfrac)  * wst(iz,jup1) + zfrac*wst(izp1,jup1)
       wbs   = (fzfrac)  * wst(iz,ju  ) + zfrac*wst(izp1,ju  )
       ws    = (1.-ufrac)* wbs          + ufrac*was
    else
       u3    = us*us*us
       wm = vonk * us * u3 / (u3 + conc1*zehat+epsln)
       ws = wm
    endif

  end subroutine wscale



  !#######################################################################
  ! <SUBROUTINE NAME="ri_iwmix">
  !
  ! Compute interior viscosity and diffusivity due 
  ! to shear instability (dependent on a local richardson number),
  ! to background internal wave activity, and 
  ! to static instability (local richardson number < 0).                    
  ! outputs:
  !    visc = viscosity coefficient (m**2/s)       
  !    diff = diffusion coefficient (m**2/s)     
  !
  subroutine ri_iwmix(viscA, diffK)

    implicit none

    real, dimension(ToDim_nod3d),   intent(out) :: viscA
    real, dimension(ToDim_nod3d,2), intent(out) :: diffK

    integer           :: i, row, k, kn, nodup, nodlo, mr
    real, parameter   :: Riinfty = 0.8  ! local Richardson Number limit for shear instability
    real              :: dz, ri_prev, tmp, fx, dens_up
    real              :: Rigg, ratio, frit, fcont
    real              :: aux, dep, lat
    logical :: smooth_richardson_number = .false.
    integer :: num_smoothings = 1 ! for vertical smoothing of Richardson number

    fx = -g*rho0r

    ! compute Richardson number and store it as viscA to save memory
    do i=1,ToDim_nod2d
       kn=num_layers_below_nod2D(i)
       do k=1, kn
          nodup=nod3d_below_nod2D(k,i)
          nodlo=nod3D_below_nod2D(k+1,i)        
          call fcn_density(tracer(nodup,1),tracer(nodup,2),coord_nod3d(3,nodlo),dens_up)
          dz=coord_nod3d(3,nodup)-coord_nod3d(3,nodlo)
          viscA(nodup)=fx*dz*(dens_up-density_insitu(nodlo)) / &
               ((uf(nodup)-uf(nodlo))**2+(uf(ToDim_nod3d+nodup)-uf(ToDim_nod3d+nodlo))**2+epsln)
       end do

       ! smooth Richardson number in the vertical using a 1-2-1 filter
       if(smooth_richardson_number .and. kn>2) then
          do mr=1,num_smoothings
             ri_prev=0.25*viscA(nod3d_below_nod2d(1,i))
             do k=2,kn-1
                nodup=nod3d_below_nod2d(k,i)
                nodlo=nod3d_below_nod2d(k+1,i)
                tmp=viscA(nodup)
                viscA(nodup)=ri_prev+0.5*viscA(nodup)+0.25*viscA(nodlo)
                ri_prev=0.25*tmp
             end do
          end do
       end if
    end do

    ! compute viscA and diffK
    do row=1, ToDim_nod3D

       !-----------------------------------------------------------------------
       !           evaluate function of Ri# for shear instability eqn. (28b&c)
       !-----------------------------------------------------------------------

       Rigg  = max( viscA(row) , 0.0)
       ratio = min( Rigg/Riinfty , 1.0 )
       frit  = (1. - ratio*ratio)
       frit  = frit*frit*frit

       !-----------------------------------------------------------------------
       !           evaluate function of Ri# for convection eqn. (28a)
       !-----------------------------------------------------------------------

       !fcont  = 0.5 * ( abs(viscA(row)) - viscA(row) )
       !fcont  = fcont / (AMAX1(fcont,epsln))

       !-----------------------------------------------------------------------
       !           add mixing due to internal wave activity (equ 29),
       !           static instability and shear instability 
       !-----------------------------------------------------------------------

       !viscA(row)     = Av0 + fcont*visc_conv_limit + visc_sh_limit*frit   
       !diffK(row,1)   = Kv0 + fcont*diff_conv_limit + diff_sh_limit*frit 
       !diffK(row,2)   = diffK(row,1)         

       ! viscosity
       viscA(row)     = Av0 + visc_sh_limit*frit 
       
       ! diffusivity
       if(Kv0_const) then
          diffK(row,1) = Kv0 + diff_sh_limit*frit 
       else
          dep=coord_nod3d(3,row)
          lat=abs(geolat(row))/rad
          if(lat<5.0) then
             ratio=1.0
          else
             ratio=min(1.0+9.0*(lat-5.0)/10.0, 10.0)
          end if
          if (geolat(row)/rad < 70.0) then
             aux=(0.6+1.0598/3.1415926*atan(4.5e-3*(abs(dep)-2500.0)))*1.e-5
          else
             aux=(0.6+1.0598/3.1415926*atan(4.5e-3*(abs(dep)-2500.0)))*1.e-6
             ratio=3.0
             if (coord_nod3d(3,row)>=-80.0) then
                ratio=1.0
             elseif(coord_nod3d(3,row)>=-100.0) then
                ratio=2.0   
             endif   
          end if
          
          diffK(row,1) = aux*ratio + diff_sh_limit*frit
       end if
       diffK(row,2)   = diffK(row,1)  

    enddo

  end subroutine ri_iwmix


  !#######################################################################
  ! <SUBROUTINE NAME="ddmix">
  !
  ! Rrho dependent interior flux parameterization.
  ! Add double-diffusion diffusivities to Ri-mix values at blending
  ! interface and below. 
  ! Salt fingering code modified july 2003
  ! by stephen.griffies@noaa.gov based on NCAR CCSM2.x
  !
  ! output: update diffu
  !
  subroutine ddmix(diffK)

    implicit none
    real, dimension(ToDim_nod3d, 2), intent(inout) :: diffK
    real, parameter :: Rrho0               = 1.9    ! limit for double diffusive density ratio
    real, parameter :: dsfmax              = 1.e-4  ! (m^2/s) max diffusivity in case of salt fingering
    real, parameter :: viscosity_molecular = 1.5e-6 ! (m^2/s)

    integer :: i, k, nodlo, nodup
    real    :: alphaDT, betaDS
    real    :: diffdd, Rrho, prandtl

    do i=1,ToDim_nod2d
       do k=1, num_layers_below_nod2d(i)

          nodup=nod3d_below_nod2d(k,i)
          nodlo=nod3d_below_nod2d(k+1,i)

          alphaDT=0.5 * (Talpha(nodup)+Talpha(nodlo)) &
               * (tracer(nodup,1)-tracer(nodlo,1))
          betaDS =0.5 * (Sbeta(nodup)+Sbeta(nodlo)) &
               * (tracer(nodup,2)-tracer(nodlo,2))

          if (alphaDT > betaDS .and. betaDS > 0.) then

             ! salt fingering case,  eqn. (31)

             Rrho   = min(alphaDT/betaDS, Rrho0)

             !        diffdd = dsfmax*(1.0-((Rrho-1)/(Rrho0-1))**2)**pexp2  ! (very old code) 
             !        diffdd = 1.0-((Rrho-1)/(Rrho0-1))**2                  ! (less old code)
             diffdd = 1.0-((Rrho-1.0)/(Rrho0-1.0))                    ! (new code)
             diffdd = dsfmax*diffdd*diffdd*diffdd

             diffK(nodup,1) = diffK(nodup,1) + 0.7*diffdd !for temperature
             diffK(nodup,2) = diffK(nodup,2) + diffdd !for salinity

          else if ( alphaDT < 0. .and. alphaDT > betaDS ) then

             ! diffusive convection eqn. (32)

             Rrho    = alphaDT/betaDS 
             diffdd  = viscosity_molecular*0.909*exp(4.6*exp(-0.54*(1./Rrho-1.)))  

             ! eqn. (34)
             prandtl = 0.15*Rrho
             if (Rrho > 0.5) prandtl = (1.85-0.85/Rrho)*Rrho

             diffK(nodup,1) = diffK(nodup,1) + diffdd  !for temperature
             diffK(nodup,2) = diffK(nodup,2) + prandtl*diffdd  !for salinity

          endif

       enddo
    enddo

  end subroutine ddmix



  !#######################################################################
  ! <SUBROUTINE NAME="blmix_kpp">
  !
  ! Mixing coefficients within boundary layer depend on surface
  ! forcing and the magnitude and gradient of interior mixing below
  ! the boundary layer ("matching").
  !
  !     inputs:
  !
  !      real ustar(2d)    ! surface friction velocity         (m/s) 
  !      real bfsfc(2d)    ! surface buoyancy forcing     (m^2/s^3)  
  !      real hbl(2d)      ! boundary layer depth              (m)   
  !      real stable(2d)   ! = 1 in stable forcing                   
  !      real caseA(2d)    ! = 1 in case A                           
  !      integer kbl(2d)   ! index of first grid level below hbl
  !
  !     outputs:
  !
  !      real dkm1(2d,3) = boundary layer diff_cbt at kbl-1 level 
  !      real blmc(3d,3) = boundary layer mixing coeff.(m**2/s)   
  !      real ghats(3d)  = nonlocal scalar transport              
  !
  subroutine blmix_kpp(viscA,diffK)

    implicit none

    real, dimension(ToDim_nod3d),   intent(inout) :: viscA
    real, dimension(ToDim_nod3d,2), intent(inout) :: diffK
    integer  :: i, nl, kn, knm1, knp1, ki, row
    integer  :: nod_col(1:max_num_layers)
    real     :: delhat, R, dvdzup, dvdzdn, wm, ws, sigma
    real     :: viscp, difsp, diftp, visch, difsh, difth, f1
    real     :: sig, a1, a2, a3, Gm, Gs, Gt
    real     :: zehat
    real     :: gat1m, gat1t, gat1s, dat1m, dat1s, dat1t
    real     :: z_col(1:max_num_layers),  z_col_mid(1:max_num_layers)
    real     :: dthick(1:max_num_layers), diff_col(1:max_num_layers,3)

    blmc = 0.0
    
    do i=1,ToDim_nod2d

       nl=num_layers_below_nod2d(i)+1

       if(nl<3) cycle  ! a temporary solution

       nod_col(1:nl)=nod3d_below_nod2d(1:nl,i)
       z_col(1:nl)=-coord_nod3d(3,nod_col(1:nl))
       z_col_mid(1:nl-1)=0.5*(z_col(1:nl-1)+z_col(2:nl))
       z_col_mid(nl)=z_col(nl)  !!
       dthick(2:nl-1)=0.5*(z_col(3:nl)-z_col(1:nl-2))
       dthick(1)=dthick(2)
       dthick(nl)=dthick(nl-1)
       diff_col(1:nl-1,1)=viscA(nod_col(1:nl-1))
       diff_col(1:nl-1,2:3)=diffK(nod_col(1:nl-1),:)
       diff_col(nl,:)=diff_col(nl-1,:)

       !-----------------------------------------------------------------------
       !       compute velocity scales at hbl. (recall epsilon=0.1)
       !-----------------------------------------------------------------------
       sigma = stable(i) * 1.0 + (1.-stable(i)) * epsilon
       zehat= vonk * sigma * hbl(i) * bfsfc(i)
       call wscale(zehat, ustar(i), wm, ws)

       kn = int(caseA(i)+epsln) *(kbl(i) -1) +   &
            (1-int(caseA(i)+epsln)) * kbl(i)

       knm1 = max(kn-1,1)
       knp1 = min(kn+1,nl)

       !-----------------------------------------------------------------------
       !         find the interior viscosities and derivatives at hbl(i) 
       !         eqn. (18)
       !-----------------------------------------------------------------------

       delhat = z_col_mid(kn)-hbl(i)
       R      = 1.0 - delhat / dthick(kn)

       dvdzup = (diff_col(knm1,1) - diff_col(kn,1))/dthick(kn)
       dvdzdn = (diff_col(kn,1) - diff_col(knp1,1))/dthick(knp1)
       viscp  = 0.5 * ( (1.-R) * (dvdzup + abs(dvdzup))+         &
            R  * (dvdzdn + abs(dvdzdn)) )

       dvdzup = (diff_col(knm1,3) - diff_col(kn,3))/dthick(kn)
       dvdzdn = (diff_col(kn,3) - diff_col(knp1,3))/dthick(knp1)     
       difsp  = 0.5 * ( (1.-R) * (dvdzup + abs(dvdzup))+         &
            R  * (dvdzdn + abs(dvdzdn)) )

       dvdzup = (diff_col(knm1,2) - diff_col(kn,2))/dthick(kn)
       dvdzdn = (diff_col(kn,2) - diff_col(knp1,2))/dthick(knp1)
       diftp  = 0.5 * ( (1.-R) * (dvdzup + abs(dvdzup))+         &
            R  * (dvdzdn + abs(dvdzdn)) )

       visch  = diff_col(kn,1) + viscp * delhat
       difsh  = diff_col(kn,3) + difsp * delhat
       difth  = diff_col(kn,2) + diftp * delhat

       f1 = stable(i) * conc1 * bfsfc(i) / (ustar(i)**4+epsln)

       gat1m = visch / (hbl(i)+epsln) / (wm+epsln)
       dat1m = -viscp / (wm+epsln) + f1 * visch
       dat1m = min(dat1m,0.) 

       gat1s = difsh  / (hbl(i)+epsln) / (ws+epsln)

       dat1s = -difsp / (ws+epsln) + f1 * difsh 
       dat1s = min(dat1s,0.) 

       gat1t = difth /  (hbl(i)+epsln) / (ws+epsln)
       dat1t = -diftp / (ws+epsln) + f1 * difth 
       dat1t = min(dat1t,0.) 

       do ki=1,nl-1

          if (ki >= kbl(i)) exit
          row=nod_col(ki)

          !-----------------------------------------------------------------------
          !         compute turbulent velocity scales on the interfaces
          !-----------------------------------------------------------------------

          sig   = z_col_mid(ki) / (hbl(i)+epsln)
          sigma = stable(i)*sig                          &
               + (1.-stable(i))*min(sig,epsilon)
          zehat= vonk * sigma * hbl(i) * bfsfc(i)
          call wscale(zehat, ustar(i), wm, ws)

          !-----------------------------------------------------------------------
          !         compute the dimensionless shape functions at the interfaces
          !         eqn. (11)
          !-----------------------------------------------------------------------

          a1 = sig - 2.
          a2 = 3.-2.*sig
          a3 = sig - 1.

          Gm = a1 + a2 * gat1m + a3 * dat1m
          Gs = a1 + a2 * gat1s + a3 * dat1s
          Gt = a1 + a2 * gat1t + a3 * dat1t

          !-----------------------------------------------------------------------
          !          compute boundary layer diffusivities at the interfaces
          !          eqn. (10)
          !-----------------------------------------------------------------------

          blmc(row,1) = hbl(i) * wm * sig * (1.+sig*Gm) 
          blmc(row,2) = hbl(i) * ws * sig * (1.+sig*Gt) 
          blmc(row,3) = hbl(i) * ws * sig * (1.+sig*Gs) 

          !-----------------------------------------------------------------------
          !          nonlocal transport term = ghats * <ws>o (eqn. 20)
          !-----------------------------------------------------------------------

          ghats(row) = (1.-stable(i)) * cg    &
               / (ws * hbl(i) + epsln)

       enddo

       !-----------------------------------------------------------------------
       !     find diffusivities at kbl-1 grid level 
       !-----------------------------------------------------------------------

       sig   =  z_col(kbl(i)-1)  / (hbl(i)+epsln)
       sigma =stable(i) * sig                  &
            + (1.-stable(i)) * min(sig,epsilon)
       zehat= vonk * sigma * hbl(i) * bfsfc(i)
       call wscale(zehat, ustar(i), wm, ws)

       a1= sig - 2.
       a2 = 3.-2.*sig
       a3 = sig - 1.
       Gm = a1 + a2 * gat1m + a3 * dat1m
       Gs = a1 + a2 * gat1s + a3 * dat1s
       Gt = a1 + a2 * gat1t + a3 * dat1t
       dkm1(i,1) = hbl(i) * wm * sig * (1. + sig * Gm)
       dkm1(i,2) = hbl(i) * ws * sig * (1. + sig * Gt)
       dkm1(i,3) = hbl(i) * ws * sig * (1. + sig * Gs)

    enddo

  end subroutine blmix_kpp


  !#######################################################################
  ! <SUBROUTINE NAME="enhance">
  !
  ! Enhance the diffusivity at the kbl-.5 interface
  !
  ! input
  !      integer kbl(n2)   =  grid above hbl                     
  !      real hbl(n2)      =  boundary layer depth (m)           
  !      real dkm1(n2,3)   =  bl diffusivity at kbl-1 grid level 
  !      real caseA(n2)    =  1 in caseA, = 0 in case B
  !
  ! input/output
  !      real ghats(ij_bounds,nk)  =  nonlocal transport     (s/m**2)
  !      modified ghats at kbl(i)-1 interface        
  ! output
  !      real blmc(n3,3) = enhanced boundary layer mixing coefficient
  !
  subroutine enhance(viscA, diffK)
    implicit none
    real, dimension(ToDim_nod3d),   intent(inout) :: viscA
    real, dimension(ToDim_nod3d,2), intent(inout) :: diffK

    integer :: i, k, nodk, nodkbl
    real    :: delta, dkmp5, dstar

    do i=1,ToDim_nod2d
       k = kbl(i) - 1
       nodk=nod3d_below_nod2d(k,i)
       nodkbl=nod3d_below_nod2d(kbl(i),i)

       delta = (hbl(i)+coord_nod3d(3,nodk))/(coord_nod3d(3,nodk)-coord_nod3d(3,nodkbl))

       ! momentum
       dkmp5 = caseA(i) * viscA(nodk)      &
            + (1.-caseA(i)) * blmc(nodk,1)
       dstar = (1.-delta)**2 * dkm1(i,1) + delta**2 * dkmp5
       blmc(nodk,1) = (1.-delta) * viscA(nodk)  &
            + delta * dstar

       ! temperature:
       dkmp5 = caseA(i) * diffK(nodk,1)  &
            + (1.-caseA(i)) * blmc(nodk,2)
       dstar = (1.-delta)**2 * dkm1(i,2) + delta**2 * dkmp5    
       blmc(nodk,2) = (1.-delta) * diffK(nodk,1)  &
            + delta * dstar

       ! salinity:   
       dkmp5 = caseA(i) * diffK(nodk,2)  &
            + (1.-caseA(i)) * blmc(nodk,3)
       dstar = (1.-delta)**2 * dkm1(i,3) + delta**2 * dkmp5
       blmc(nodk,3) = (1.-delta) * diffK(nodk,2)  &
            + delta * dstar

       ghats(nodk) = (1.-caseA(i)) * ghats(nodk)
    enddo

  end subroutine enhance
  !
  !----------------------------------------------------------------------------
  !
  subroutine cal_nodal_alpha_beta
    !   A function to calculate the thermal expansion coefficient
    !   and saline contraction coefficient. 
    !
    ! REFERENCE:
    !    McDougall, T.J. 1987.  Neutral Surfaces
    !    Journal of Physical Oceanography, vol 17, 1950-1964,
    !
    ! INPUT:
    !   tracer(:,1) = potential temperature [degree C (ITS-90)]
    !   tracer(:,2) = salinity              [psu      (PSS-78)]
    !   z           = pressure (or -depth)  [db]
    !
    ! OUTPUT:
    !   Talpha = Thermal expansion coeff (alpha) [degree_C.^-1]
    !   Sbeta  = Saline contraction coeff (beta) [psu.^-1]
    !
    ! Qiang Wang, 25,11,2004
    ! Modified to compute nodal values for KPP, Feb. 2011, Qiang
    !-----------------------------------------------------------------
    ! CHECK VALUE:
    !    sw_beta=0.72088e-3 psu^-1 @ S=40.0psu, ptmp=10.0C (ITS-90), p=4000db
    !    a_over_b=0.34765 psu*C^-1 @ S=40.0psu, ptmp=10.0C, p=4000db
    !-----------------------------------------------------------------

    implicit none

    integer    :: i, ni
    real       :: t1,t1_2,t1_3,t1_4,p1,p1_2,p1_3,s1,s35,s35_2 
    real       :: a_over_b    

    ! we should ensure that we are evoluting potential temperature in the model.

    ! cycle over node
    do i = 1, ToDim_nod3d
       ! prepare values that will be used below
       t1 = tracer(i,1)*1.00024_8
       s1 = tracer(i,2)
       p1 = abs(coord_nod3D(3,i)) 

       t1_2 = t1*t1
       t1_3 = t1_2*t1
       t1_4 = t1_3*t1
       p1_2 = p1*p1
       p1_3 = p1_2*p1
       s35 = s1-35.0_8
       s35_2 = s35*s35

       ! calculate beta
       Sbeta(i) = 0.785567e-3 - 0.301985e-5*t1 &
            + 0.555579e-7*t1_2 - 0.415613e-9*t1_3 &
            + s35*(-0.356603e-6 + 0.788212e-8*t1 &
            + 0.408195e-10*p1 - 0.602281e-15*p1_2) &
            + s35_2*(0.515032e-8) & 
            + p1*(-0.121555e-7 + 0.192867e-9*t1 - 0.213127e-11*t1_2) &
            + p1_2*(0.176621e-12 - 0.175379e-14*t1) &
            + p1_3*(0.121551e-17)

       ! calaculate the thermal expansion / saline contraction ratio
       a_over_b = 0.665157e-1 + 0.170907e-1*t1 &
            - 0.203814e-3*t1_2 + 0.298357e-5*t1_3 &
            - 0.255019e-7*t1_4 &
            + s35*(0.378110e-2 - 0.846960e-4*t1 &
            - 0.164759e-6*p1 - 0.251520e-11*p1_2) &
            + s35_2*(-0.678662e-5) &
            + p1*(0.380374e-4 - 0.933746e-6*t1 + 0.791325e-8*t1_2) &
            + p1_2*t1_2*(0.512857e-12) &
            - p1_3*(0.302285e-13)

       ! calculate alpha
       Talpha(i) = a_over_b*Sbeta(i)
    end do
  end subroutine cal_nodal_alpha_beta

end module o_mixing_KPP_mod

! These routines were created to assist the debugging of KPP. They check whether an array has been communicated or not. May be helpfull in the future so we leave them here.
subroutine is_communicated2(arr)
  use o_MESH
  use o_ELEMENTS
  use g_config
  use o_PARAM
  use o_array
  use g_forcing_arrays
  use g_PARFE

  implicit none
  save

  real(kind=8), dimension(ToDim_nod2D), intent(in)  :: arr
  real(kind=8), dimension(:), allocatable           :: aux
  real(kind=8)                                      :: gmin, gmax

  allocate(aux(ToDim_nod2D))
  aux=arr
  call com_2D(aux)
  call MPI_AllREDUCE(minval(aux-arr), gmin, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_FESOM, MPIerr)
  call MPI_AllREDUCE(maxval(aux-arr), gmax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_FESOM, MPIerr)
  if (mype==0) write(*,*) 'is communicated 2D min/max=', gmin, gmax
  deallocate(aux)

end subroutine is_communicated2

subroutine is_communicated3(arr)

  use o_MESH
  use o_ELEMENTS
  use g_config
  use o_PARAM
  use o_array
  use g_forcing_arrays
  use g_PARFE

  implicit none
  save

  real(kind=8), dimension(ToDim_nod3D), intent(in)    :: arr
  real(kind=8), dimension(:), allocatable             :: aux
  real(kind=8)                                        :: gmin, gmax

  allocate(aux(ToDim_nod3D))

  aux=arr
  call com_3D(aux)
  call MPI_AllREDUCE(minval(aux-arr), gmin, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_FESOM, MPIerr)
  call MPI_AllREDUCE(maxval(aux-arr), gmax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_FESOM, MPIerr)
  if (mype==0) write(*,*) 'is communicated 3D min/max=', gmin, gmax
  deallocate(aux)

end subroutine is_communicated3
