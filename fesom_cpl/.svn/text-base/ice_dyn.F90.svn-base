subroutine stress_tensor
  ! Internal stress tensor
  ! This version does not differentiate the metric terms. 
  ! This should be changed! Qiang
  !
  ! Coded by Sergey Danilov and N. Yakovlev
  ! Reasoning and irony by Ralph Timmermann, March 3. 06
  ! Reviewed by Qiang Wang
  !===================================================================
  
  use i_dyn_parms
  use i_therm_parms
  use i_array
  use o_mesh
  use o_elements 
  use o_param
  use g_parfe
  use g_config
  implicit none

  integer        :: i, elem, elnodes(3)
  real(kind=8)   :: dx(3), dy(3)
  real(kind=8)   :: eps11, eps12, eps22, eta, xi, pressure, delta, aa
  real(kind=8)   :: elcos(3), val3, meancos, usum, vsum, vale, xi_limit
  real(kind=8)   :: det1, det2, r1, r2, r3, si1, si2, dte

  Tevp_inv=dt_inv*real(evp_Tdamp_ratio)
  val3=1.0_8/3.0_8
  vale=1.0_8/(ellipse**2)

  if(EVP_rheology) then 
     dte=dt/(1.0_8*evp_rheol_steps)
     det1=1.0_8+0.5_8*Tevp_inv*dte
     !det2=1.0_8+0.5_8*Tevp_inv*dte*ellipse**2    !RTSD corrected 8.3.2006
     det2=1.0_8+0.5_8*Tevp_inv*dte 
     det1=1.0_8/det1
     det2=1.0_8/det2
  else
     xi_limit=0.5*rhoice*9.e8*vp_rheol_steps/dt 
     ! Limit moduli to satisfy the CFL - type criterium in Euler
     ! forward time stepping;  9.0e8 is the squared mesh step ??
  end if


  do elem=1,myDim_elem2d                    
     elnodes=elem2D_nodes(:,elem)

     aa=product(m_ice(elnodes))*product(a_ice(elnodes))
     if (aa==0.0) cycle

     dx=bafux_2d(:,elem)
     dy=bafuy_2d(:,elem)     
     !     vsum=sum(v_ice(elnodes))
     !     usum=sum(u_ice(elnodes))
     !     elcos=cos(coord_nod2D(2,elnodes))
     !     meancos=cos_elem2D(elem)

     ! Deformation rate tensor on element elem:
     eps11=sum(dx*u_ice(elnodes))
     !     eps11=eps11+val3*vsum*sum(dy*elcos)/meancos
     eps22=sum(dy*v_ice(elnodes))
     eps12=0.5_8*sum(dy*u_ice(elnodes) + dx*v_ice(elnodes))
     !     eps12=eps12-0.5_8*val3*usum*sum(dy*elcos)/meancos  

     ! moduli:
     delta=(eps11**2+eps22**2)*(1.0_8+vale)+4.0_8*vale*eps12**2 + &
          2.0_8*eps11*eps22*(1.0_8-vale)
     delta=sqrt(delta)
    
     vsum=sum(m_ice(elnodes))*val3
     usum=sum(a_ice(elnodes))*val3

     pressure=pstar*vsum*exp(-c_pressure*(1.0_8-usum))
    
     if (.not.EVP_rheology) then   
        xi=0.5_8*pressure/(delta+delta_min)

        if(xi>xi_limit*vsum) then              ! Limiting moduli 
           xi=xi_limit*vsum
        end if

        eta=xi*vale
        pressure=0.5_8*pressure*delta/(delta+delta_min)

        ! stress arrays 
        sigma11(elem)=(xi-eta)*(eps11+eps22)-0.5_8*pressure
        sigma12(elem)=2.0_8*eta*eps12
        sigma22(elem)=2.0_8*eta*eps22+sigma11(elem)
        sigma11(elem)=2.0_8*eta*eps11+sigma11(elem) 
     else
        pressure=0.5_8*pressure*Tevp_inv
        pressure=pressure*delta
        delta=1.0_8/(delta+delta_min)
        pressure=pressure*delta

        r1=pressure*((eps11+eps22)*delta -1.0_8) 
        r2=pressure*(eps11-eps22)*delta/ellipse**2
        r3=pressure*eps12*delta/ellipse**2
        si1=sigma11(elem)+sigma22(elem)
        si2=sigma11(elem)-sigma22(elem)

        si1=det1*(si1+dte*r1)
        si2=det2*(si2+dte*r2)
        sigma12(elem)=det2*(sigma12(elem)+dte*r3)
        sigma11(elem)=0.5_8*(si1+si2)
        sigma22(elem)=0.5_8*(si1-si2)
     end if
  end do

end subroutine stress_tensor
!
!===================================================================
!
subroutine stress2rhs
  ! add internal stress and ssh gradient terms to the rhs
  ! Metric terms neglected.
  ! This should be changed! Qiang
  ! 
  ! Coded by Sergey Danilov, March 2006
  ! Reviewed by Ralph Timmermann
  ! Modified by Qiang Wang, 16.12.2010, add the nonlinear free surface option
  !-----------------------------------------------------------------
  
  use o_MESH
  use o_ELEMENTS
  use o_param
  use i_dyn_parms
  use i_therm_parms
  use i_array
  use g_PARFE
  use g_clock
  use g_config
  implicit none
  
  integer       :: i ,k, row, elem, elnodes(3)
  real(kind=8)  :: dx(3), dy(3), vol
  real(kind=8)  :: val3, val12, meancos, elcos(3), aa, aux(3)
  real(kind=8)  :: mass, cluster_area, elevation_elem(3)

  val3=1.0_8/3.0_8
  val12=1.0_8/12.0_8

  do row=1, myDim_nod2d 
     rhs_u(row)=0.0
     rhs_v(row)=0.0
     rhs_a(row)=0.0
     rhs_m(row)=0.0
  end do

  do elem=1,myDim_elem2d         
     elnodes=elem2D_nodes(:,elem)

     aa=product(m_ice(elnodes))*product(a_ice(elnodes))
     if (aa==0.0) cycle

     vol=voltriangle(elem)
     dx=bafux_2d(:,elem)
     dy=bafuy_2d(:,elem)     
     elevation_elem=elevation(elnodes)
#ifdef use_fullfreesurf
     aux=(rhoice*m_ice(elnodes)+rhosno*m_snow(elnodes))*rho0r
     do i=1,3
        aux(i)=min(aux(i),max_ice_loading)
     end do
     elevation_elem=elevation_elem+aux
     !rho0r is used here to be consistent with the computation of ocean hydrostatic pressure
#endif

     !elcos=cos(coord_nod2D(2,elnodes))
     !meancos=cos_elem2D(elem)

     do k=1,3
        row=elnodes(k)
        rhs_u(row)=rhs_u(row) - vol* &
             (sigma11(elem)*dx(k)+sigma12(elem)*dy(k))   ! ) &
        !+val3*sum(dy*elcos)/meancos))
        rhs_v(row)=rhs_v(row) - vol* &
             (sigma12(elem)*dx(k)+sigma22(elem)*dy(k))   ! ) &
        ! +sigma11(elem)*val3*sum(dy*elcos)/meancos)
     end do

     ! use rhs_m and rhs_a for storing the contribution from elevation:
     aa=g*val3*vol
     do k=1,3
        row=elnodes(k)
        rhs_a(row)=rhs_a(row)-aa*sum(dx*elevation_elem)	    
        rhs_m(row)=rhs_m(row)-aa*sum(dy*elevation_elem)
     end do
  end do

  do row=1, myDim_nod2d               
     cluster_area=ice_lump(row)*dt  !note: ice_lump contains dt_inv
     mass=cluster_area*(m_ice(row)*rhoice+m_snow(row)*rhosno)

     if (mass>1.0e-4) then
        rhs_u(row)=rhs_u(row)/mass + rhs_a(row)/cluster_area 
        rhs_v(row)=rhs_v(row)/mass + rhs_m(row)/cluster_area 
     else
        rhs_u(row)=0.  
        rhs_v(row)=0.
     end if
  end do

  ! Clean neighbours:
  do row=myDim_nod2d+1,eDim_nod2d+myDim_nod2d  
     rhs_u(row)=0. 
     rhs_v(row)=0.
  end do

end subroutine stress2rhs
!
!===================================================================
!
subroutine rheology
  ! assemble rhs and solve ice velocity
  !
  ! Coded by Sergey Danilov and N. Yakovlev 
  ! Reviewed by Ralph Timmermann and Qiang Wang
  !---------------------------------------------------------

  use o_param
  use o_MESH
  use o_ELEMENTS
  use o_array, only: coriolis_param_nod2D
  use i_dyn_parms
  use i_therm_parms
  use i_array
  use g_parfe
  use g_config

  implicit none
  integer         :: steps, shortstep, i, j
  real(kind=8)    :: rdt, drag, det, fc
  real(kind=8)    :: thickness, inv_thickness, umod, rhsu, rhsv


  if(.not.EVP_rheology) then
     rdt=dt/(1.0*vp_rheol_steps)
     steps=vp_rheol_steps
  else
     rdt=dt/(1.0*evp_rheol_steps)
     steps=evp_rheol_steps
  end if

  do shortstep=1, steps 
     ! ===== Boundary condition      
     do i=1, toDim_nod2D    
        j=nod3D_below_nod2D(1,i)  
#ifdef use_cavity    
        if(index_nod3D(j)==11 .or. cavity_flag_extended(i)==1) then 
#else
        if(index_nod3D(j)==11) then     
#endif
           u_ice(i)=0.0
           v_ice(i)=0.0
        end if
     end do

     call stress_tensor
     call stress2rhs
     
     do i=1,myDim_nod2d 
        !in case only myDim, 2d and 3d nodes are the same!
#ifdef use_cavity        
        if (index_nod3D(i)==11 .or. cavity_flag_extended(i)==1) cycle 
#else
        if (index_nod3D(i)==11) cycle       ! Skip boundary nodes
#endif
        if (a_ice(i) > Armin) then          ! If ice is present, update velocities
           thickness=(rhoice*m_ice(i)+rhosno*m_snow(i))/a_ice(i)
           thickness=max(thickness, 90.0)   ! Limit the weighted thickness if it is too small
           inv_thickness=1.0/thickness

           umod=sqrt((u_ice(i)-u_w(i))**2+(v_ice(i)-v_w(i))**2)
           drag=Cd_oce_ice*umod*rho0*inv_thickness

           !rhs for water stress, air stress, and rhs_u/v (internal stress + ssh)
           rhsu=u_ice(i)+drag*rdt*u_w(i)+rdt*(inv_thickness*stress_atmice_x(i)+rhs_u(i))
           rhsv=v_ice(i)+drag*rdt*v_w(i)+rdt*(inv_thickness*stress_atmice_y(i)+rhs_v(i))

           !solve (coriolis and water stress are treated implicitly)
           fc=coriolis_param_nod2D(i)
           det=(1.+drag*rdt)**2+(rdt*fc)**2
           det=1.0_8/det
           u_ice(i)=det*((1.0+drag*rdt)*rhsu+rdt*fc*rhsv)
           v_ice(i)=det*((1.0+drag*rdt)*rhsv-rdt*fc*rhsu)

           !elseif(a_ice(i) > 0.001)          ! Set ice velocity equal to water velocity
           !   u_ice(i)=u_w(i)             
           !   v_ice(i)=v_w(i)
        end if
     end do
     call com_2D(u_ice)
     call com_2D(v_ice)
  end do
end subroutine rheology
!
!===================================================================
!
subroutine cut_off
  use i_array
  use i_therm_parms
  implicit none
  !
  real(kind=8)   :: m_min, a_min 
  !
  a_min=1.0e-3
  m_min=1.0e-3

  !where ((m_ice< m_min).or.(a_ice<a_min))
  !   a_ice=0.
  !   m_ice=0.
  !end where

  where(a_ice>1.0_8)
     a_ice=1.0_8
  end where

  where(a_ice<0.0_8)
     a_ice=0.0_8
  end where
end subroutine cut_off
!
!===================================================================


