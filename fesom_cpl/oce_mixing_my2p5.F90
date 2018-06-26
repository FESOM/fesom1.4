module o_mixing_my2p5_mod
  ! Mellor-Yamada 2.5 mixing scheme (Mellor and Yamada, 1982)
  ! The Cantha and Clayson (1994) stability functions are used.
  ! Constrains on turbulent scales following Galperin etal (1988) are enforced.
  ! Here we took a simple/approximate way to deal with advection and diffusion.
  !
  ! It was coded following the structure of POM with changes in detail. Further
  ! modification was then undertaken to use the Cantha and Clayson stability 
  ! functions. I am still not very sure about the treatment of the upper and 
  ! bottom boundaries. A check-through for it is required at a later stage.
  !
  ! The background values (Kv0 and Av0) should be added when operating vertical
  ! mixing for the momentum and tracer equations. They are traditionally set to
  ! Kv0~1e-5 and Av0~1e-4. Define them in the namelist file.
  !
  ! Qiang, 06.08.2010

  use o_mesh
  use o_elements
  use o_matrices
  use o_param
  use o_array 
  use g_config
  use g_parfe

  implicit none

  public oce_mixing_my2p5_init
  public oce_mixing_my2p5

  real, allocatable, dimension(:) 	  :: Kq
  real, allocatable, dimension(:) 	  :: q2, q2b, q2l, q2lb
  real, allocatable, dimension(:) 	  :: rhs_q2, rhs_q2l

contains

  subroutine oce_mixing_my2p5_init
    ! initialize the MY2p5 scheme
    implicit none

    integer	    :: m, row, k, kb, nodbot, nodk
    real(kind=8)    :: l, smallvalue

    data smallvalue /1.0e-8/

    allocate(Kq(ToDim_nod3d))
    allocate(q2(ToDim_nod3d), q2b(ToDim_nod3d))
    allocate(q2l(Todim_nod3d), q2lb(ToDim_nod3d))
    allocate(rhs_q2(ToDim_nod3d), rhs_q2l(ToDim_nod3d))
    do row=1,ToDim_nod2d
       kb=num_layers_below_nod2d(row)+1
       nodbot=bt_nds(row)
       l=0.1*(coord_nod3d(3,row)-coord_nod3d(3,nodbot))
       do k=1,kb
          nodk=nod3d_below_nod2d(k,row)
          q2b(nodk)=smallvalue
          q2lb(nodk)=l*q2b(nodk)
          Kv(nodk,1)=l*sqrt(q2b(nodk))
          Av(nodk)=Kv(nodk,1)
          Kq(nodk)=Kv(nodk,1)
          q2(nodk)=q2b(nodk)
          q2l(nodk)=q2lb(nodk)
       end do
    end do

  end subroutine oce_mixing_my2p5_init
  !
  !------------------------------------------------------------------------------
  !
  subroutine oce_mixing_my2p5
    ! update MY2.5

    implicit none

    integer   	:: m, row, n, i, j, elnodes(4), elem, elem2
    integer	:: k, ki, kbm1, kb, lay, n3, rowd, elem_type
    integer	:: nodbot, nodk, nodkm1, nodkp1
    integer	:: nod_col(max_num_layers)
    real(kind=8)        :: dt2, stress_sur, stress_bot, ustr, vstr
    real(kind=8)        :: aux, dudz, dvdz, v, vtr, mass
    real(kind=8)	:: dx(4), dy(4), dz(4), qdx, qdy, qdz
    real(kind=8)	:: um, vm, wm, u_el(4), v_el(4), q_el(4)
    real(kind=8)	:: Ah, d1, d2, S(3)
    real(kind=8)	:: rotate_coe, temp_coe, temp_coe2, temp_coe3
    real(kind=8)	:: zb, z0, z, tp, sp, p, kn, Gh
    real(kind=8)        :: dens_k, grad_rho, grad_rho_up
    real(kind=8)	:: a1, a2, b1, b2, c1, c2, c3
    real(kind=8)	:: e1, e2, l0, stf, smallvalue
    real(kind=8)	:: cbcnst, surfl, sef, smoth, inv4
    real(kind=8)        :: coef1, coef2, coef3, coef4, coef5
    real(kind=8)	:: a(max_num_layers), c(max_num_layers)
    real(kind=8)	:: ee(max_num_layers), gg(max_num_layers)
    real(kind=8)	:: rhs_q2_col(max_num_layers), rhs_q2l_col(max_num_layers)
    real(kind=8)        :: dz_col(max_num_layers), l(max_num_layers)
    real(kind=8)        :: boygr(max_num_layers), cc(max_num_layers)
    real(kind=8)        :: prod(max_num_layers), q_col(max_num_layers)
    real(kind=8)	:: dtef(max_num_layers)
    real(kind=8)	:: sh(max_num_layers), sm(max_num_layers), maxgrad
    real(kind=8)	:: kappa !von Karman's constant

    data kappa /0.4/
    data a1,b1,a2,b2,c1,c2,c3 /0.92, 16.6, 0.74, 10.1, 0.08, 0.7, 0.2/
    data e1/1.8/, e2/1.33/
    data cbcnst/100./, surfl/2.e5/, sef/1./

    data smoth/0.05/, inv4/0.25/
    data smallvalue /1.0e-8/

    dt2=2.0*dt

    coef1=a2*(1.0-6.0*a1/b1)
    coef2=3.0*a2*(6.0*a1+b2*(1.0-c3))
    coef3=a1*(1.0-3.0*c1-6.0*a1/b1)
    coef4=a1*9.0*(2.0*a1+a2*(1-c2))
    coef5=9.0*a1*a2

    n3=myDim_nod3D+eDim_nod3D              

    rhs_q2=0.0
    rhs_q2l=0.0

    ! assemble rhs for q2 and q2l (adv. and hori. diff.)
    ! use leapfrog for adv, Euler for diff.
    ! use lumped mass matrix.
    do elem=1, myDim_elem3d        
       elem2=elem2d_corresp_to_elem3d(elem)
       elem_type=grid_type_elem2d(elem2)
       elnodes=elem3D_nodes(:,elem)
       dx=bafux_3D(:,elem)
       dy=bafuy_3D(:,elem)
       dz=bafuz_3D(:,elem)
       v=voltetra(elem)

       ! mean velocity  
       u_el=uf(elnodes)
       v_el=uf(elnodes+n3)                       
       um=sum(u_el)*inv4
       vm=sum(v_el)*inv4
#ifndef use_non_hydrostatic
       wm=sum(w(elnodes)*dz)
#else
       wm=sum(uf(elnodes+2*n3))*inv4              
#endif

       ! horizontal viscosity 
       if(smagorinsky_visc) then
          vtr=voltriangle(elem2)
          d1=sum(u_el*dx)-sum(v_el*dy)
          d2=sum(u_el*dy)+sum(v_el*dx)
          Ah=sqrt(d1*d1+d2*d2)*vtr*4.0
          Ah=max(Ah, 40.0/1.0e8*vtr)              !limit from below: 20m^2/s on 10km
          if(Ah*dt/vtr > 0.05) Ah=0.05*vtr*dt_inv !limit from above
       else
          Ah=Ah0
          if(scale_mixing_h) Ah=Ah*(voltriangle(elem2)/scalevol)**(1.0/real(scale_mixing_type))
       endif

       ! rhs for q2
       ! adv
       q_el=q2(elnodes)
       qdx=sum(dx*q_el)
       qdy=sum(dy*q_el)
       qdz=sum(dz*q_el)
       rhs_q2(elnodes)=rhs_q2(elnodes) - (um*qdx+vm*qdy+wm*qdz)*inv4*v
       ! diff 
       q_el=q2b(elnodes)
       qdx=sum(dx*q_el)
       qdy=sum(dy*q_el)
       qdz=sum(dz*q_el)
       if(elem_type==1) then  !sigma diff.
          lay=elem3d_layer(elem)
          S(1)=grid_slope(1,lay,elem2)
          S(2)=grid_slope(2,lay,elem2)
          aux=S(1)**2+S(2)**2
          S(3)=sqrt(aux)
          rotate_coe=1.0/(1.0+aux) 

          !diagonal part1 (lateral)
          temp_coe=Ah*rotate_coe*v
          rhs_q2(elnodes)=rhs_q2(elnodes)-(dx*qdx*(1.0+S(2)*S(2))+ &
               dy*qdy*(1.0+S(1)*S(1)))*temp_coe
          !diagonal part2 (cross slope)
          temp_coe2=S(3)*S(3)*temp_coe
          rhs_q2(elnodes)=rhs_q2(elnodes)-dz*qdz*temp_coe2
          !off diagonal part1 (lateral) --> (1,3),(2,3)
          rhs_q2(elnodes)=rhs_q2(elnodes)-(S(1)*dx+S(2)*dy)*qdz*temp_coe
          !off diagonal part1 (lateral)--> (1,2),(2,1), '- -' -> '+'
          temp_coe3=S(1)*S(2)*temp_coe
          rhs_q2(elnodes)=rhs_q2(elnodes) + (dx*qdy+dy*qdx)*temp_coe3    
          !off diagonal part2 (cross slope) --> (3,1),(3,2)
          rhs_q2(elnodes)=rhs_q2(elnodes)-(S(1)*qdx+S(2)*qdy)*dz*temp_coe
       else
          rhs_q2(elnodes)=rhs_q2(elnodes) - Ah*(dx*qdx+dy*qdy)*v
       endif

       ! rhs for q2l
       ! adv
       q_el=q2l(elnodes)
       qdx=sum(dx*q_el)
       qdy=sum(dy*q_el)
       qdz=sum(dz*q_el)
       rhs_q2l(elnodes)=rhs_q2l(elnodes) - (um*qdx+vm*qdy+wm*qdz)*inv4*v
       ! diff 
       q_el=q2lb(elnodes)
       qdx=sum(dx*q_el)
       qdy=sum(dy*q_el)
       qdz=sum(dz*q_el)
       if(elem_type==1) then  !sigma diff.
          !diagonal part1 (lateral)
          rhs_q2l(elnodes)=rhs_q2l(elnodes)-(dx*qdx*(1.0+S(2)*S(2))+ &
               dy*qdy*(1.0+S(1)*S(1)))*temp_coe
          !diagonal part2 (cross slope)
          rhs_q2l(elnodes)=rhs_q2l(elnodes)-dz*qdz*temp_coe2
          !off diagonal part1 (lateral) --> (1,3),(2,3)
          rhs_q2l(elnodes)=rhs_q2l(elnodes)-(S(1)*dx+S(2)*dy)*qdz*temp_coe
          !off diagonal part1 (lateral)--> (1,2),(2,1), '- -' -> '+'
          rhs_q2l(elnodes)=rhs_q2l(elnodes) + (dx*qdy+dy*qdx)*temp_coe3    
          !off diagonal part2 (cross slope) --> (3,1),(3,2)
          rhs_q2l(elnodes)=rhs_q2l(elnodes)-(S(1)*qdx+S(2)*qdy)*dz*temp_coe
       else
          rhs_q2l(elnodes)=rhs_q2l(elnodes) - Ah*(dx*qdx+dy*qdy)*v
       endif
    end do ! elem

    ! solve for rhs_q2 and rhs_q2l
    rhs_q2 = q2b + dt2*rhs_q2/uv_lump    
    rhs_q2l = q2lb + dt2*rhs_q2l/uv_lump

    ! add other terms in the q2 and q2l equations
    do row=1,myDim_nod2D                    

       kbm1=num_layers_below_nod2D(row)
       kb=kbm1+1

       ! nodes and dz in this colume
       nod_col(1)=nod3d_below_nod2d(1,row)
       do k=2, kb
          nod_col(k)=nod3d_below_nod2d(k,row)
          dz_col(k-1)=coord_nod3d(3,nod_col(k-1))-coord_nod3d(3,nod_col(k))
       end do

       ! surface and bottom depth and nodes
       rowd=nod_col(1) 
       z0=coord_nod3d(3,rowd)  
       nodbot=nod3d_below_nod2d(kb,row)
       zb=coord_nod3d(3,nodbot)

       !rhs_q2 and rhs_q2l in this colume
       rhs_q2_col(1:kb)=rhs_q2(nod_col(1:kb))
       rhs_q2l_col(1:kb)=rhs_q2l(nod_col(1:kb))

       !-----------------------------------------------------------------------
       !Surface and bottom boundary conditions

       ! amplitude of surface and bottom stress
       if(cavity_flag_nod2d(row)==0) then       
          stress_sur=sqrt(stress_x(row)**2+stress_y(row)**2)*rho0r
       else            
          aux=sqrt(uf(rowd)**2+uf(n3+rowd)**2)       
          ustr=-C_d*uf(rowd)*aux                      
          vstr=-C_d*uf(n3+rowd)*aux                 
          stress_sur=sqrt(ustr**2+vstr**2)
       end if
       aux=sqrt(uf(nodbot)**2+uf(n3+nodbot)**2)      
       ustr=-C_d*uf(nodbot)*aux
       vstr=-C_d*uf(n3+nodbot)*aux                   
       stress_bot=sqrt(ustr**2+vstr**2)        

       ee(1)=0.0
       aux=(16.6**(2.0/3.0))*sef
       if(cavity_flag_nod2d(row)==0) then      
          ! Wave breaking energy- a variant of Craig & Banner (1994)
          ! see Mellor and Blumberg, 2004.
          gg(1)=(15.8*cbcnst)**(2.0/3.0)*stress_sur 
          ! Surface length scale following Stacey (1999).
          l0=surfl*stress_sur/g
       else  !under the cavity
          gg(1)=stress_sur*aux   ! similar to ocean bottom
          l0=0.0
       end if
       rhs_q2_col(kb)=stress_bot*aux     


       !-----------------------------------------------------------------------
       ! The following section solves the equation:
       ! dti2*(kq*q2')' - q2*(2.*dti2*dtef+1.) = -q2b

       ! coefficients a, c
       do k=2, kbm1
          nodk=nod_col(k)
          nodkm1=nod_col(k-1)
          nodkp1=nod_col(k+1)
          aux=dz_col(k)+dz_col(k-1)
          a(k)=-dt2*(Kq(nodk)+Kq(nodkp1)+2.0*Kv0)/(dz_col(k)*aux)
          c(k)=-dt2*(Kq(nodkm1)+Kq(nodk)+2.0*Kv0)/(dz_col(k-1)*aux)
       end do

       ! the original way of computing buoyancy frequency in POM
       ! Calculate speed of sound:
       !do k=2,kbm1
       !   nodk=nod_col(k)
       !   tp=TF(nodk)
       !   sp=SF(nodk)
       !   ! pressure in units of decibars
       !   ! p=abs(coord_nod3d(3,nodk))*g*rho0*1.0e-4
       !   p=abs(coord_nod3d(3,nodk))
       !   cc(k)=1449.1+0.00821*p+4.55*tp-0.045*tp**2+1.34*(sp-35.0)
       !   cc(k)=cc(k)/sqrt((1.0-0.01642*p/cc(k))*(1.0-0.4*p/cc(k)**2))
       !end do
       !
       ! Calculate buoyancy gradient:
       !do k=2,kbm1
       !   nodk=nod_col(k)
       !   nodkm1=nod_col(k-1)
       !   nodkp1=nod_col(k+1)
       !   boygr(k)=g*0.5*((density_insitu(nodkm1)-density_insitu(nodk))/dz_col(k-1) + &
       !	(density_insitu(nodk)-density_insitu(nodkp1))/dz_col(k))*rho0r
       !   if(.not.density_linear) boygr(k)=boygr(k)+(g**2)/(cc(k)**2)       
       !end do

       ! my way of calculating buoyancy frequency:
       nodk=nod_col(1)
       nodkp1=nod_col(2)
       call fcn_density(tracer(nodk,1), tracer(nodk,2), coord_nod3d(3,nodkp1), dens_k)
       grad_rho_up=(dens_k-density_insitu(nodkp1))/dz_col(1)
       boygr(1)=g*rho0r*grad_rho_up
       do k=2,kbm1
          nodk=nod_col(k)
          nodkp1=nod_col(k+1)
          call fcn_density(tracer(nodk,1), tracer(nodk,2), coord_nod3d(3,nodkp1), dens_k)
          grad_rho=(dens_k-density_insitu(nodkp1))/dz_col(k)
          boygr(k)=g*rho0r*0.5*(grad_rho+grad_rho_up)
          grad_rho_up=grad_rho
       end do
       boygr(kb)=g*rho0r*grad_rho_up
       ! boygr(1) and boygr(kb) are only required for limiting l

       ! Calculate production of turbulent kinetic energy:
       do k=2,kbm1
          nodk=nod_col(k)
          nodkm1=nod_col(k-1)
          nodkp1=nod_col(k+1)
          dudz=(((uf(nodkm1)-uf(nodk))/dz_col(k-1))**2 + &
               ((uf(nodk)-uf(nodkp1))/dz_col(k))**2)*0.5
          dvdz=(((uf(n3+nodkm1)-uf(n3+nodk))/dz_col(k-1))**2 + &   
               ((uf(n3+nodk)-uf(n3+nodkp1))/dz_col(k))**2)*0.5      
          prod(k)=Av(nodk)*sef*(dudz+dvdz)  
          prod(k)=prod(k) + Kv(nodk,1)*boygr(k)       
       end do

       ! length scale at time step  n
       l(1)=kappa*l0
       do k=2,kbm1
          nodk=nod_col(k)
          q2b(nodk)=max(q2b(nodk), smallvalue)    !necessary to use q2b and q2lb
          q2lb(nodk)=max(q2lb(nodk),smallvalue)   !to derive l (not q2 and q2l)
          l(k)=q2lb(nodk)/q2b(nodk)           
          if((z0-coord_nod3d(3,nodk))<0.5*(z0-zb)) l(k)=max(l(k), l(1))  
       end do
       l(kb)=0.0

       !dissipation rate: dtef*q2
       do k=2, kbm1
          nodk=nod_col(k)
          dtef(k)=sqrt(q2(nodk))/(b1*l(k)+smallvalue)
       end do

       ! compute ee, gg
       ! ee(1) gg(1) was set boundary condition; ee(kb) gg(kb) is not required.		
       do k=2,kbm1
          gg(k)=1.0/(a(k)+c(k)*(1.0-ee(k-1))-2.0*dt2*dtef(k)-1.0)
          ee(k)=a(k)*gg(k)
          gg(k)=(-2.0*dt2*prod(k)-rhs_q2_col(k)+c(k)*gg(k-1))*gg(k)
       end do

       ! compute q2 at n+1
       ! rhs_q2_col(kb) was set boundary condition
       do k=1,kbm1
          ki=kb-k
          rhs_q2_col(ki)=ee(ki)*rhs_q2_col(ki+1)+gg(ki)
       end do

       !-----------------------------------------------------------------------
       ! The following section solves the equation:
       ! dti2(kq*q2l')' - q2l*(dti2*dtef+1.) = -q2lb

       !boundary 
       ee(1)=0.0  
       gg(1)=0.0
       rhs_q2l_col(kb)=0.0
       !ee(2)=0.0
       !gg(2)=kappa*(z0-coord_nod3d(3,nod_col(2)))*q2(nod_col(2))
       !rhs_q2l_col(kbm1)=kappa*(coord_nod3d(3,nod_col(kbm1))-zb)*q2(nod_col(kbm1))

       ! dissipation rate
       do k=2,kbm1 
          z=coord_nod3d(3,nod_col(k))
          dtef(k)=dtef(k)*(1.0+e2*((1.0/abs(z-z0)+1.0/abs(z-zb))*l(k)/kappa)**2)   
       end do

       ! compute ee, gg 
       ! ee(1) gg(1) was set boundary condition; ee(kb) gg(kb) is not required.	
       do k=2, kbm1 !3
          gg(k)=1.0/(a(k)+c(k)*(1.0-ee(k-1))-dt2*dtef(k)-1.0)
          ee(k)=a(k)*gg(k)
          gg(k)=(-dt2*prod(k)*l(k)*e1-rhs_q2l_col(k)+c(k)*gg(k-1))*gg(k)
       end do

       ! compute q2l at n+1
       ! rhs_q2l_col(kb) was set boundary condition	
       do k=1,kb-1  !kb-2        
          ki=kb-k
          rhs_q2l_col(ki)=ee(ki)*rhs_q2l_col(ki+1)+gg(ki)
       end do

       ! avoid negative values of q2 and q2l
       do k=1,kb
          if(rhs_q2_col(k)<smallvalue .or. rhs_q2l_col(k)<smallvalue) then
             rhs_q2_col(k)=smallvalue
             rhs_q2l_col(k)=smallvalue
          end if
       end do

       !-----------------------------------------------------------------------
       ! The following section solves for Av and Kv:

       ! new length scale
       l(2:kbm1)=rhs_q2l_col(2:kbm1)/rhs_q2_col(2:kbm1)
       ! constrain l according to Galperin etal. 1988
       q_col(1:kb)=sqrt(rhs_q2_col(1:kb))
       do k=1,kb
          if(boygr(k)<0.0) then
             l(k)=min(l(k), 0.53*q_col(k)/(sqrt(-boygr(k))+smallvalue))
          end if
       end do

       ! stability functions: Cantha and Clayson, 1994
       do k=1,kb
          Gh=l(k)**2/rhs_q2_col(k)*boygr(k)
          Gh=min(max(Gh,-0.28), 0.029)
          sh(k)=coef1/(1.0-coef2*Gh)
          sm(k)=coef3+sh(k)*coef4*Gh
          sm(k)=sm(k)/(1.0-coef5*Gh)
       end do

       ! POM note: There are 2 options for Kq which, unlike Av and Kv, 
       ! was not derived by Mellor and Yamada but was purely
       ! empirical based on neutral boundary layer data.
       ! The choice is whether or not it should be subject to
       ! the stability factor, sh. Generally, there is not a great
       ! difference in output. 
       ! Q:The first was mentioned/used by Mellor and Blumberg, 2004; 
       ! the second was used by Cantha and Clayson 1994)

       ! POM takes the mean of the current step and the last step for
       ! Av, Kv and Kq. I only take the new value as I did no see the 
       ! reason to do that.
       ! POM uses old l and q2 to compute Gh and kn; I decided to
       ! use the new l and q2 as I did not see the reason to do that.
       ! If we recognize what we did is wrong, modification is required here!

       do k=1,kb
          nodk=nod_col(k)
          kn=l(k)*q_col(k)
          !Kq(nodk)=kn*0.41*sh(k)
          Kq(nodk)=kn*0.2
          Av(nodk)=kn*sm(k)
          Kv(nodk,1)=kn*sh(k)
       end do

       !-----------------------------------------------------------------------
       ! Average Kv and Av to values at the mid-edges (between every two nodes)
       ! because currently the implicit vertical mixing equations require such values. 
       ! Storage convention: store a mean value to the upper node location of the edge
       do k=1,kbm1
          nodk=nod_col(k)
          Av(nodk)=0.5*(Av(nodk)+Av(nod_col(k+1)))
          Kv(nodk,1)=0.5*(Kv(nodk,1)+Kv(nod_col(k+1),1))
       end do

       !----------------------------------------------------------------------- 
       ! updating arrays q2 q2b q2l q2lb and filtering (Robert-Asselin time filter)
       do k=1,kb
          nodk=nod_col(k)
          q2(nodk)=q2(nodk)+0.5*smoth*(rhs_q2_col(k)+q2b(nodk)-2.0*q2(nodk))
          q2b(nodk)=q2(nodk)
          q2(nodk)=rhs_q2_col(k)
          q2l(nodk)=q2l(nodk)+0.5*smoth*(rhs_q2l_col(k)+q2lb(nodk)-2.0*q2l(nodk))
          q2lb(nodk)=q2l(nodk)
          q2l(nodk)=rhs_q2l_col(k)
       end do

    end do

    ! communicate partition neighbours
    call com_3d(Av)
    call com_3d(Kv(:,1))
    call com_3d(q2)
    call com_3d(q2b)
    call com_3d(q2l)
    call com_3d(q2lb)

  end subroutine oce_mixing_MY2p5

end module o_mixing_MY2p5_mod
