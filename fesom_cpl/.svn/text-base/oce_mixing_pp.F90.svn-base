module o_mixing_pp_mod
  !1)Richardson number dependent Av and Kv following
  !  Pacanowski and Philander, 1981.
  !2)TB04 mixing scheme
  !
  ! Coded by Ralph Timmermann and Sergey Danilov
  ! Reviewed by Qiang Wang
  !----------------------------------------------------------

  use o_MESH
  use o_ELEMENTS
  use o_PARAM
  use o_array
  use i_array
  use g_config
  use g_PARFE

  implicit none
  save
  private

  public oce_mixing_pp
  public oce_mixing_pp_init
 
  real, allocatable, dimension(:)  :: mo_mixlength

contains
  

  subroutine oce_mixing_pp_init
    !initialize the pp scheme
    !Here only need to allocate mo_mixlength used by the tb04 scheme
    implicit none

    allocate(mo_mixlength(toDim_nod2d))
    mo_mixlength=0.0
  end subroutine oce_mixing_pp_init
  !
  !------------------------------------------------------------------------------
  !
  subroutine oce_mixing_pp
    !Compute Richardson number dependent Av and Kv following
    !Pacanowski and Philander, 1981
    !Av = c * factor**2 + Av0,
    !Kv = Av * factor + Kv0, 
    !factor=1/(1+10Ri),  Ri=N**2/(dU/dz)**2 is the Richardson number
    !                    N is buoyancy frequency
    !                   dU/dz is vertical velocity shear
    !c is a tunable parameter (named mix_coeff_PP in the code)
    !Kv0 ~ 1.e-5, Av0 ~ 1.e-4

    implicit none
    !
    integer     :: col,n,j,node_low,node_up, node, k, n3
    real   	:: density_up, drho, g_rhoinv
    real	:: dz_inv, dzu, dzv, zm
    real	:: velocity_shear, dzb, factor
    logical     :: mo_on 

    g_rhoinv=g*rho0r 
    n3=ToDim_nod3d

    do col=1,ToDim_nod2d

       n=num_layers_below_nod2D(col)

       do j=1,n       
          node=nod3D_below_nod2D(j,col)
          node_up = node
          node_low = nod3D_below_nod2D(j+1,col)

          dz_inv=1./(coord_nod3d(3,node_up)-coord_nod3d(3,node_low))

          call fcn_density(tracer(node_up,1),tracer(node_up,2), &
               coord_nod3d(3,node_low),density_up)

          ! density vert. grad.
          drho = density_up - density_insitu(node_low)

          if(drho>=0.0) then
             factor=1.0 
          else
             ! BV frequency squared
             dzb=drho*dz_inv*g_rhoinv  
             ! velocity shear
             dzu=(uf(node_up)- uf(node_low)) * dz_inv
             dzv=(uf(node_up+n3) - uf(node_low+n3)) * dz_inv 
             velocity_shear=dzu*dzu + dzv*dzv
             factor=velocity_shear/( velocity_shear - 5.0 * dzb)      
          end if

          Av(node)=PP_max_mix_coeff*(factor**2)
          Kv(node,1)=Av(node)*factor

          ! add background mixing
          Av(node)=Av(node)+Av0
          Kv(node,1)=Kv(node,1)+Kv0

          ! turn on convection if unstable 
          if(drho>0.0) then
             if(allow_convect_global) then
                Kv(node,1)=diff_conv_limit
                Av(node)=visc_conv_limit
             elseif(coord_nod2d(2,col)>0.0) then     
                Kv(node,1)=diff_conv_limit
                Av(node)=visc_conv_limit
             endif
          end if

          ! apply the tb04 scheme
#ifdef use_ice
#ifdef use_cavity
          mo_on=(add_tb04_to_PP .and. cavity_flag_nod2d(col)==0)
#else
          mo_on=add_tb04_to_PP
#endif
          if(mo_on) then
             call mo_length(water_flux(col),heat_flux(col), &
                  stress_x(col),stress_y(col), &
                  u_ice(col),v_ice(col),a_ice(col), &
                  dt, mo_mixlength(col))
             zm=0.5*(coord_nod3d(3,node_up)+coord_nod3d(3,node_low))
             if(abs(zm)<=mo_mixlength(col)) then
                Av(node)=Av(node)+modiff
                Kv(node,1)=Kv(node,1)+modiff
             end if
          end if
#endif

       end do
    end do

    ! approximation for high freq. wind mixing near the surface
    ! (in case not in the forcing)
    where(Av(1:myDim_nod2D)<wndmix) 
       Av(1:myDim_nod2d)=wndmix 
    end where
    where(Kv(1:myDim_nod2D,1)<wndmix) 
       Kv(1:myDim_nod2d,1)=wndmix 
    end where

  end subroutine oce_mixing_PP
  !
  !---------------------------------------------------------------------------
  !
  subroutine mo_length(water_flux,heat_flux,stress_x,stress_y,  &
       u_ice,v_ice,a_ice,dt,mixlength)
    ! Vertical mixing scheme of Timmermann and Beckmann, 2004.
    ! Computes the mixing length derived from the Monin-Obukhov length
    ! Coded by Ralph Timmermann, 14.06.2006

    implicit none

    real*8              :: water_flux, heat_flux, stress_x, stress_y
    real*8              :: u_ice, v_ice, a_ice, uabs
    real*8              :: dt, ret, rtc, mixlength
    real*8              :: qfm, qtm, qw
    real*8              :: ustar,tau, obuk
    real(kind=8), parameter :: cosgam = 0.913632  ! cos(24.*3.14/180.)

    qfm            = water_flux * 34.       ! note that water_flux>0
    ![psu * m/s]   [m/s]   [psu]    ! f. upward fresh water flux

    qtm            = - 2.38e-7 * heat_flux  ! heat_flux>0 f. upward heat flux
    ![K * m/s]

    tau = sqrt(stress_x**2+stress_y**2)
    ustar = sqrt(tau/1030.)
    uabs = sqrt(u_ice**2+v_ice**2)

    qw = 1.25 * ustar**3 * (1-a_ice) + 0.005 * uabs**3 * cosgam * a_ice  !Eq. 8 of TB04

    call pmlktmo(qfm,qtm,qw,obuk)

    rtc=10.0*86400.0/dt     ! time constant of mixed layer retreat

    if (obuk.lt.mixlength) then
       ret=(obuk-mixlength)/rtc
       mixlength=mixlength+ret
    else
       mixlength=obuk
    endif

  end subroutine mo_length
  !
  !=========================================================================	    	    
  !
  subroutine pmlktmo(qfm,qtm,qw,obuk)
    ! gives the Monin-Obukhov length
    ! qtm  = Heat Flux into ML                                 [K m/s]
    ! qfm  = salinity flux into ML                             [psu m/s]
    ! qw   = production of turbulent kinetic energy 
    ! Coded by Ralph Timmermann, 14.06.2006
    !-----------------------------------------------------------------------
    implicit none

    integer           :: iter
    real*8            :: qtm, qfm, qw, obuk
    real*8, parameter :: qhw   = 7.0                           ! [m]
    real*8, parameter :: betas = 0.0008
    real*8, parameter :: betat = 0.00004
    real*8            :: a1, f0, f1, ttmp, qrho

    qrho=betas*qfm-betat*qtm
    ttmp=60.

    if(qrho>0.) then
       ttmp=0.
    else
       do iter=1,5
          a1 = exp(-ttmp/qhw) 
          f0 = 2.* qw * a1 + 9.81 * qrho * ttmp
          f1 =-2.* qw * a1 / qhw + 9.81 * qrho
          ttmp = ttmp - f0 / f1
          ttmp = max(ttmp,10.) 
       enddo
    end if

    obuk=max(ttmp,10.)

  end subroutine pmlktmo

end module o_mixing_pp_mod


