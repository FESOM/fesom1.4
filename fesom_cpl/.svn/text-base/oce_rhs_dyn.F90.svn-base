! Assemble rhs for ocean dynamics, and solving routines

subroutine uv_sfc_bott_bc
  ! momentum surface and bottom boundary conditions
  ! This routine is only required when using implicit vertical mixing scheme.
  ! In case of the explicit scheme, the assembling is done directly within
  ! the velocity_rhs routine. This is different from the tracer bc treatment.
  !
  ! Coded by Qiang Wang
  ! Reviewed by ??
  !-------------------------------------------------------------
  
  use o_param
  use o_array
  use o_mesh
  use o_elements
  use g_config
  use g_parfe
  implicit none

  integer      :: elem2, q, row, elnodes2(3), elnodes2_3d(3)
  real(kind=8) :: vol, usum, vsum, aux, uvel(3), vvel(3)
  real(kind=8) :: inv3, inv12

  inv3=1./3.
  inv12=1./12.

  uv_sfc_force=0.
  uv_bott_force=0.

  do elem2=1, myDim_elem2d    
     vol=voltriangle(elem2)
     elnodes2=elem2D_nodes(:,elem2)

     !surface             
     !Wind/ice stress, or stress from ice-cavity    
     usum=sum(stress_x(elnodes2))
     vsum=sum(stress_y(elnodes2))
     aux=rho0r*inv12*vol
     do q=1,3
        row=elnodes2(q)
        uv_sfc_force(row,1)=uv_sfc_force(row,1) + aux*(usum+stress_x(row)) 
        uv_sfc_force(row,2)=uv_sfc_force(row,2) + aux*(vsum+stress_y(row)) 
     end do

     !bottom 
     !Bottom drag contribution  -C_d U|U| 
     elnodes2_3d=bt_nds(elnodes2)        
     uvel=uf(elnodes2_3d)
     vvel=uf(elnodes2_3d+ToDim_nod3D)                           
     usum=sum(uvel)
     vsum=sum(vvel)
     aux=-C_d*sqrt(usum**2+vsum**2)*inv3*inv12*vol
     do q=1,3
        row=elnodes2(q)
        uv_bott_force(row,1)=uv_bott_force(row,1) + aux*(usum+uvel(q)) 
        uv_bott_force(row,2)=uv_bott_force(row,2) + aux*(vsum+vvel(q))        
     end do
  end do

  do q=1, myDim_nod2d  !eDim_nod2d is not required 
     if(mod(index_nod3d(q),10)==1) then
     	uv_sfc_force(q,1)=0.
	uv_sfc_force(q,2)=0.
	!uv_bott_force(q,1)=0.
	!uv_bott_force(q,2)=0.
     endif
     if(mod(index_nod3d(bt_nds(q)),10)==1) then
        uv_bott_force(q,1)=0.
	uv_bott_force(q,2)=0.
     endif
  end do

end subroutine uv_sfc_bott_bc
!
!----------------------------------------------------------------------------
!
subroutine velocity_rhs
  ! Assemble rhs for momentum equ.
  !
  ! Coded by Qiang Wang
  ! Reviewed by ??
  !-------------------------------------------------------------

  use o_param
  use o_array
  use o_mesh
  use o_elements
  use g_config
  use g_parfe
  implicit none
  !
  integer                        :: m, i, elem, elem2, elnodes(4), elnodes23(4)
  integer                        :: q, row, row2, elnodes2(3), ind, lay
  integer                        :: n3, elem_type
  real(kind=8)                   :: fc, aux, adv, coe, Ah, Av_el, d1, d2
  real(kind=8)                   :: nhp_el(4), Ah_sma
  real(kind=8)                   :: px_ssh, py_ssh, uvel(4), vvel(4)
  real(kind=8)                   :: usum, vsum, um, vm, wm
  real(kind=8)                   :: pxsum ,pysum
  real(kind=8)                   :: udx, udy, udz, vdx, vdy, vdz
  real(kind=8)                   :: dx(4), dy(4), dz(4), v, vtr, val
  real(kind=8)                   :: inv5, inv4, inv3, inv12, inv20
  real(kind=8)                   :: urhs_temp, vrhs_temp, cori_u_tmp, cori_v_tmp
  real(kind=8)                   :: S(3), rotate_coe, temp_coe, temp_coe2
  real(kind=8)                   :: dparam, beta, gamma
#ifdef use_non_hydrostatic
  integer                        :: row3
  real(kind=8)                   :: wvel(4), wsum, npz, npx, npy
  real(kind=8)                   :: wdx, wdy, wdz, wrhs_temp
#endif
#ifdef use_fullfreesurf
  integer                        :: n_el
  real(kind=8)                   :: coe_diff, uusum, uvsum, vvsum, wmove(4)
  real(kind=8)                   :: wmove_sum, uwmove_sum, vwmove_sum
  logical                        :: flag_move=.false.
#ifdef use_non_hydrostatic
  real(kind=8)                   :: wusum, wvsum, wwsum, wwmove_sum
#endif
#endif

  inv5=0.2_8
  inv4=0.25_8
  inv3=1.0_8/3.0_8
  inv12=1.0_8/12.0_8
  inv20=0.05_8
  n3=ToDim_nod3d
  do row=1, myDim_nod3d                   
     row2=row+n3                          
     uv_rhs(row)=0.
     uv_rhs(row2)=0.
#ifdef use_non_hydrostatic
     row3=row2+n3                       
     uv_rhs(row3)=0.
#endif
  enddo

  !3D contributions	 
  do elem=1, myDim_elem3d       
     elnodes=elem3D_nodes(:,elem)
     elem2 = elem2D_corresp_to_elem3D(elem)
     elnodes2=elem2D_nodes(:,elem2)
     elnodes23=nod2D_corresp_to_nod3D(elnodes)
     dx=bafux_3D(:,elem)
     dy=bafuy_3D(:,elem)
     dz=bafuz_3D(:,elem)
     v=voltetra(elem)
     elem_type=grid_type_elem2d(elem2)

     !velocity
     uvel=uf(elnodes)
     vvel=uf(elnodes+n3)                  
     usum=sum(uvel)         
     vsum=sum(vvel)
     um=usum*inv4
     vm=vsum*inv4
#ifndef use_non_hydrostatic
     wm=sum(w(elnodes)*dz)
#else
     wvel=uf(elnodes+2*n3)                
     wsum=sum(wvel)
     wm=wsum*inv4
#endif

     !coriolis parameter
     fc=coriolis_param_elem2D(elem2)
     !terms for semi-implicit schemes
     if(use_cori_semi) then
        val=v
        dparam=dt_inv**2 + alpha_trapez**2*fc**2
        beta=dt_inv/dparam
        gamma=fc*alpha_trapez/dparam
     else
        val=v*dt
     endif

     !ssh/non-hydrostatic p
     nhp_el=g*ssh(elnodes23)
     px_ssh=sum(dx*nhp_el)
     py_ssh=sum(dy*nhp_el)
#ifdef use_non_hydrostatic
     nhp_el=nhp(elnodes)
     npx=sum(dx*nhp_el)
     npy=sum(dy*nhp_el)
     npz=sum(dz*nhp_el)
#endif

     !visc/adv
     udx=sum(dx*uvel)
     udy=sum(dy*uvel)
     udz=sum(dz*uvel)
     vdx=sum(dx*vvel)
     vdy=sum(dy*vvel)
     vdz=sum(dz*vvel)
#ifdef use_non_hydrostatic
     wdx=sum(dx*wvel)
     wdy=sum(dy*wvel)
     wdz=sum(dz*wvel)
#endif

     !viscosity
     if(.not.biharmonic_visc) then
        Ah=Ah0
        vtr=voltriangle(elem2)
        if(scale_mixing_h) Ah=Ah*(vtr/scalevol)**(1.0/real(scale_mixing_type))
        if(smagorinsky_visc) then
           d1=udx-vdy
           d2=udy+vdx
           Ah_sma=sqrt(d1*d1+d2*d2)*vtr                  !set C/pi to 1
           !Ah_sma=max(Ah_sma, 40.0/1.0e8*vtr)           !limit from below: 20m^2/s on 10km
           if(Ah_sma*dt/vtr > 0.03) Ah_sma=0.03*vtr/dt   !limit from above
           Ah=Ah+Ah_sma
        endif
        !     if(any(mod(index_nod3d(elnodes),10)==1)) then
        !        Ah=Ah*2.0
        !     end if
     end if

     if(use_vertvisc_impl) then
        Av_el=0.0
     else
        Av_el=sum(Av(elnodes))/4.0
     endif

     !required variables for free surface
#ifdef use_fullfreesurf
     uusum=sum(uvel*uvel)
     vvsum=sum(vvel*vvel)
     uvsum=sum(uvel*vvel)
#ifdef use_non_hydrostatic
     wusum=sum(wvel*uvel)
     wvsum=sum(wvel*vvel)
     wwsum=sum(wvel*wvel)
#endif
     n_el=map_elem(elem)
     flag_move=.false.
     if(n_el/=0) then
        flag_move=.true.
        coe_diff=inv20*(v-voltetra_new(n_el))
        wmove=0.
        do i=1, 4   
           if(index_nod3d(elnodes(i))>12) cycle
           row=nod2d_corresp_to_nod3d(elnodes(i))
           wmove(i)=-(ssh(row)-ssh0(row))*dt_inv 
        end do
        wmove_sum=sum(wmove)
        uwmove_sum=sum(uvel*wmove)
        vwmove_sum=sum(vvel*wmove)
#ifdef use_non_hydrostatic
        wwmove_sum=sum(wvel*wmove)
#endif
     end if
#endif

     !hp
     if(elem_type==0) then
#ifdef use_fullfreesurf
        if(n_el/=0) then
           pxsum=sum(bafux_3d_fix(:,n_el)*hpressure(elnodes))
           pysum=sum(bafuy_3d_fix(:,n_el)*hpressure(elnodes)) 
        else
           pxsum=sum(dx*hpressure(elnodes))
           pysum=sum(dy*hpressure(elnodes)) 
        end if
#else
        pxsum=sum(dx*hpressure(elnodes))
        pysum=sum(dy*hpressure(elnodes)) 
#endif
     else 
        lay=elem3d_layer(elem) 
        pxsum=PGF(1,lay,elem2)                
        pysum=PGF(2,lay,elem2)
     end if


     !for along sigma mixing
     if (elem_type==1) then
        lay=elem3d_layer(elem)      
        S(1)=grid_slope(1,lay,elem2)
        S(2)=grid_slope(2,lay,elem2)
        aux=S(1)**2+S(2)**2
        S(3)=sqrt(aux)
        rotate_coe=1.0/(1.0+aux) 
     end if


     coe=inv20*val
     do q=1,4
        row=elnodes(q)
        row2=row+n3                       
#ifdef use_non_hydrostatic
        row3=row2+n3                         
#endif
        aux=(um*dx(q)+vm*dy(q)+wm*dz(q))*dt*0.5_8*val
#ifdef use_fullfreesurf
        if(flag_move) then
           aux=aux+wmove_sum*inv4*dz(q)*dt*0.5_8*val  
        end if
#endif  
        !
#ifdef use_fullfreesurf
        if(flag_move) then
           uv_rhs(row)=uv_rhs(row)+(usum+uvel(q))*coe_diff
           uv_rhs(row2)=uv_rhs(row2)+(vsum+vvel(q))*coe_diff
#ifdef use_non_hydrostatic
           uv_rhs(row3)=uv_rhs(row3)+(wsum+wvel(q))*coe_diff
#endif   
        endif
#endif
        !coriolis force
        cori_u_tmp=fc*(vsum+vvel(q))*coe
        cori_v_tmp=-fc*(usum+uvel(q))*coe 
        !stabilization - coriolis
        cori_u_tmp=cori_u_tmp+aux*fc*vm
        cori_v_tmp=cori_v_tmp-aux*fc*um	   

        !hydrostatic pressure
        urhs_temp=-pxsum*val*inv4
        vrhs_temp=-pysum*val*inv4
        !stabilization - hp
        urhs_temp=urhs_temp-aux*pxsum
        vrhs_temp=vrhs_temp-aux*pysum

        !ssh/non-hydrostatic pressure
#ifndef use_non_hydrostatic
        urhs_temp=urhs_temp-inv4*val*px_ssh*gamma_stab
        vrhs_temp=vrhs_temp-inv4*val*py_ssh*gamma_stab
        !stabilization - ssh
        urhs_temp=urhs_temp-aux*px_ssh
        vrhs_temp=vrhs_temp-aux*py_ssh
#else
        urhs_temp=urhs_temp-inv4*val*(px_ssh*gamma_stab+npx*gamma_stab_nh)
        vrhs_temp=vrhs_temp-inv4*val*(py_ssh*gamma_stab+npy*gamma_stab_nh)
        wrhs_temp=-inv4*val*npz*gamma_stab_nh
        !stabilization - ssh/nhp
        urhs_temp=urhs_temp-aux*(px_ssh+npx)
        vrhs_temp=vrhs_temp-aux*(py_ssh+npy)
        wrhs_temp=wrhs_temp-aux*npz
#endif

        !viscosity
        if(.not.biharmonic_visc) then
           if(elem_type==1) then  !along-sigma viscosity
              !diagonal part1 (lateral)
              temp_coe=Ah*rotate_coe*val
              urhs_temp=urhs_temp-(dx(q)*udx*(1.0+S(2)*S(2))+ &
                   dy(q)*udy*(1.0+S(1)*S(1)))*temp_coe
              vrhs_temp=vrhs_temp-(dx(q)*vdx*(1.0+S(2)*S(2))+ &
                   dy(q)*vdy*(1.0+S(1)*S(1)))*temp_coe
#ifdef use_non_hydrostatic
              wrhs_temp=wrhs_temp-(dx(q)*wdx*(1.0+S(2)*S(2))+ &
                   dy(q)*wdy*(1.0+S(1)*S(1)))*temp_coe
#endif
              !diagonal part2 (cross slope)
              temp_coe2=Av_el*val+S(3)*S(3)*temp_coe
              urhs_temp=urhs_temp-dz(q)*udz*temp_coe2
              vrhs_temp=vrhs_temp-dz(q)*vdz*temp_coe2
#ifdef use_non_hydrostatic
              wrhs_temp=wrhs_temp-dz(q)*wdz*temp_coe2
#endif
              !off diagonal part1 (lateral) --> (1,3),(2,3)
              urhs_temp=urhs_temp-(S(1)*dx(q)+S(2)*dy(q))*udz*temp_coe
              vrhs_temp=vrhs_temp-(S(1)*dx(q)+S(2)*dy(q))*vdz*temp_coe
#ifdef use_non_hydrostatic
              wrhs_temp=wrhs_temp-(S(1)*dx(q)+S(2)*dy(q))*wdz*temp_coe
#endif
              !off diagonal part1 (lateral)--> (1,2),(2,1), '- -' -> '+'
              temp_coe2=S(1)*S(2)*temp_coe
              urhs_temp=urhs_temp + (dx(q)*udy+dy(q)*udx)*temp_coe2
              vrhs_temp=vrhs_temp + (dx(q)*vdy+dy(q)*vdx)*temp_coe2
#ifdef use_non_hydrostatic
              wrhs_temp=wrhs_temp + (dx(q)*wdy+dy(q)*wdx)*temp_coe2
#endif
              !off diagonal part2 (cross slope) --> (3,1),(3,2)
              urhs_temp=urhs_temp-(S(1)*udx+S(2)*udy)*dz(q)*temp_coe
              vrhs_temp=vrhs_temp-(S(1)*vdx+S(2)*vdy)*dz(q)*temp_coe 
#ifdef use_non_hydrostatic
              wrhs_temp=wrhs_temp-(S(1)*wdx+S(2)*wdy)*dz(q)*temp_coe 
#endif
           else
              urhs_temp=urhs_temp-(Ah*dx(q)*udx + &
                   Ah*dy(q)*udy + Av_el*dz(q)*udz)*val
              vrhs_temp=vrhs_temp-(Ah*dx(q)*vdx + &
                   Ah*dy(q)*vdy + Av_el*dz(q)*vdz)*val
#ifdef use_non_hydrostatic
              wrhs_temp=wrhs_temp-(Ah*dx(q)*wdx + &
                   Ah*dy(q)*wdy+Av_el*dz(q)*wdz)*val
#endif
           endif
        end if

        !advection
#ifndef use_non_hydrostatic
#ifndef use_fullfreesurf
        adv=(usum+uvel(q))*udx*inv5+(vsum+vvel(q))*udy*inv5+wm*udz
        urhs_temp=urhs_temp-adv*inv4*val
        adv=(usum+uvel(q))*vdx*inv5+(vsum+vvel(q))*vdy*inv5+wm*vdz 
        vrhs_temp=vrhs_temp-adv*inv4*val
#else
        adv=(dx(q)*(usum*usum+uusum)+dy(q)*(vsum*usum+uvsum))*inv20
        adv=adv+dz(q)*wm*usum*inv4
        if(flag_move) adv=adv+dz(q)*(wmove_sum*usum+uwmove_sum)*inv20
        urhs_temp=urhs_temp+adv*val
        adv=(dx(q)*(usum*vsum+uvsum)+dy(q)*(vsum*vsum+vvsum))*inv20
        adv=adv+dz(q)*wm*vsum*inv4
        if(flag_move) adv=adv+dz(q)*(wmove_sum*vsum+vwmove_sum)*inv20
        vrhs_temp=vrhs_temp+adv*val
#endif
        !stabilization - adv
        adv=um*udx+vm*udy+wm*udz
        urhs_temp=urhs_temp-aux*adv
        adv=um*vdx+vm*vdy+wm*vdz
        vrhs_temp=vrhs_temp-aux*adv
#ifdef use_fullfreesurf
        if(flag_move) then
           urhs_temp=urhs_temp-aux*wmove_sum*inv4*udz
           vrhs_temp=vrhs_temp-aux*wmove_sum*inv4*vdz
        end if
#endif 
#else
#ifndef use_fullfreesurf
        adv=(usum+uvel(q))*udx+(vsum+vvel(q))*udy+(wsum+wvel(q))*udz
        urhs_temp=urhs_temp-adv*coe
        adv=(usum+uvel(q))*vdx+(vsum+vvel(q))*vdy+(wsum+wvel(q))*vdz
        vrhs_temp=vrhs_temp-adv*coe
        adv=(usum+uvel(q))*wdx+(vsum+vvel(q))*wdy+(wsum+wvel(q))*wdz
        wrhs_temp=wrhs_temp-adv*coe
#else
        adv=dx(q)*(usum*usum+uusum)+dy(q)*(vsum*usum+uvsum) + &
             dz(q)*(wsum*usum+wusum)
        if(flag_move) adv=adv+dz(q)*(wmove_sum*usum+uwmove_sum)
        urhs_temp=urhs_temp+adv*coe
        adv=dx(q)*(usum*vsum+uvsum)+dy(q)*(vsum*vsum+vvsum) + &
             dz(q)*(wsum*vsum+wvsum)
        if(flag_move) adv=adv+dz(q)*(wmove_sum*vsum+vwmove_sum)
        vrhs_temp=vrhs_temp+adv*coe
        adv=dx(q)*(usum*wsum+wusum)+dy(q)*(vsum*wsum+wvsum) + &
             dz(q)*(wsum*wsum+wwsum)
        if(flag_move) adv=adv+dz(q)*(wmove_sum*wsum+wwmove_sum)
        wrhs_temp=wrhs_temp+adv*coe
#endif
        !stabilization - adv
        adv=um*udx+vm*udy+wm*udz
        urhs_temp=urhs_temp-aux*adv
        adv=um*vdx+vm*vdy+wm*vdz
        vrhs_temp=vrhs_temp-aux*adv
        adv=um*wdx+vm*wdy+wm*wdz
        wrhs_temp=wrhs_temp-aux*adv
#ifdef use_fullfreesurf
        if(flag_move) then
           urhs_temp=urhs_temp-aux*wmove_sum*inv4*udz
           vrhs_temp=vrhs_temp-aux*wmove_sum*inv4*vdz
           wrhs_temp=wrhs_temp-aux*wmove_sum*inv4*wdz
        end if
#endif 
#endif

        !put rhs temp values to the uv_rhs array
        if(use_cori_semi) then
           urhs_temp=urhs_temp+cori_u_tmp
           vrhs_temp=vrhs_temp+cori_v_tmp
           uv_rhs(row)=uv_rhs(row)+beta*urhs_temp+gamma*vrhs_temp
           uv_rhs(row2)=uv_rhs(row2)+beta*vrhs_temp-gamma*urhs_temp
#ifdef use_non_hydrostatic
           uv_rhs(row3)=uv_rhs(row3)+wrhs_temp*dt
#endif
        else
           uv_rhs(row)=uv_rhs(row)+urhs_temp
           uv_rhs(row2)=uv_rhs(row2)+vrhs_temp
           ucori(row)=ucori(row)+cori_u_tmp
           vcori(row)=vcori(row)+cori_v_tmp
#ifdef use_non_hydrostatic
           uv_rhs(row3)=uv_rhs(row3)+wrhs_temp
#endif
        endif
     end do
  end do


  !2D contributions: surface forcing and bottom stress
  if(.not.use_vertvisc_impl) then
     do elem2=1, myDim_elem2d    
        v=voltriangle(elem2)
        if(use_cori_semi) then
           val=v
           !coriolis parameter
           fc=coriolis_param_elem2D(elem2)
           !terms for semi-implicit schemes
           dparam=dt_inv**2 + alpha_trapez**2*fc**2
           beta=dt_inv/dparam
           gamma=fc*alpha_trapez/dparam
        else
           val=v*dt
        endif

        !surface             
        !Wind/ice stress contribution, or stress from ice-cavity    
        elnodes2=elem2D_nodes(:,elem2)   
#ifdef use_cavity
        if(all(cavity_flag_nod2d(elnodes2)==0)) then   
#endif
           usum=sum(stress_x(elnodes2))
           vsum=sum(stress_y(elnodes2))
           aux=rho0r*inv12*val
           do q=1,3
              row=elnodes2(q)
              urhs_temp=aux*(usum+stress_x(row)) 
              vrhs_temp=aux*(vsum+stress_y(row))
              row=nod3D_below_nod2D(1,elnodes2(q))       
              row2=row+n3                                 
              if(use_cori_semi) then
                 uv_rhs(row)=uv_rhs(row)+beta*urhs_temp+gamma*vrhs_temp
                 uv_rhs(row2)=uv_rhs(row2)+beta*vrhs_temp-gamma*urhs_temp
              else
                 uv_rhs(row)=uv_rhs(row)+urhs_temp
                 uv_rhs(row2)=uv_rhs(row2)+vrhs_temp
              endif
           end do
#ifdef use_cavity
        else
           uvel(1:3)=uf(nod3D_below_nod2D(1,elnodes2))      
           vvel(1:3)=uf(nod3D_below_nod2D(1,elnodes2)+n3)   
           usum=sum(uvel(1:3))
           vsum=sum(vvel(1:3))
           aux=-C_d*sqrt(usum**2+vsum**2)*inv3*inv12*val
           do q=1,3
              row=nod3D_below_nod2D(1,elnodes2(q))          
              row2=row+n3                                   
              urhs_temp=aux*(usum+uvel(q))
              vrhs_temp=aux*(vsum+vvel(q))
              if(use_cori_semi) then
                 uv_rhs(row)=uv_rhs(row)+beta*urhs_temp+gamma*vrhs_temp
                 uv_rhs(row2)=uv_rhs(row2)+beta*vrhs_temp-gamma*urhs_temp
              else
                 uv_rhs(row)=uv_rhs(row)+urhs_temp
                 uv_rhs(row2)=uv_rhs(row2)+vrhs_temp
              endif
           end do
        end if
#endif

        !bottom 
        !Bottom drag contribution  -C_d U|U| 
        elnodes2=bt_nds(elnodes2)        
        uvel(1:3)=uf(elnodes2)
        vvel(1:3)=uf(elnodes2+n3)                           
        usum=sum(uvel(1:3))
        vsum=sum(vvel(1:3))
        aux=-C_d*sqrt(usum**2+vsum**2)*inv3*inv12*val
        do q=1,3
           row=elnodes2(q)
           row2=row+n3                                     
           urhs_temp=aux*(usum+uvel(q))
           vrhs_temp=aux*(vsum+vvel(q))
           if(use_cori_semi) then
              uv_rhs(row)=uv_rhs(row)+beta*urhs_temp+gamma*vrhs_temp
              uv_rhs(row2)=uv_rhs(row2)+beta*vrhs_temp-gamma*urhs_temp
           else
              uv_rhs(row)=uv_rhs(row)+urhs_temp
              uv_rhs(row2)=uv_rhs(row2)+vrhs_temp
           endif
        end do
     end do
  end if


  if(biharmonic_visc) call biharmonic_viscosity


  !update A-B part (coriolis) and vertical boundary conditions
  do row=1, myDim_nod3d              
     row2=row+n3                           
     if(.not.use_cori_semi) then
        uv_rhs(row)=uv_rhs(row)+alpha_AB*ucori(row)+(1.0-alpha_AB)*ucori_back(row)
        uv_rhs(row2)=uv_rhs(row2)+alpha_AB*vcori(row)+(1.0-alpha_AB)*vcori_back(row)
        ucori_back(row)=ucori(row)
        vcori_back(row)=vcori(row)
        ucori(row)=0.
        vcori(row)=0.
     endif
     ! boundary conditions
     ind=index_nod3d(row)
     if((ind==11).or.(ind==21).or.(ind==31)) then
        uv_rhs(row)=0.0
        uv_rhs(row2)=0.0
#ifdef use_non_hydrostatic
        uv_rhs(row2+n3)=0.0              
#endif 
     end if
  end do
end subroutine velocity_rhs
!
!----------------------------------------------------------------------------
!
subroutine velocity_rhs_update
  ! rhs for updating velocity after computing ssh
  !
  ! Coded by Qiang Wang
  ! Reviewed by ??
  !-------------------------------------------------------------
  
  use o_param
  use o_array
  use o_mesh
  use o_elements
  use g_config
  use g_parfe
  implicit none
  !
  integer                        :: m, elem, el2, elnodes(4), elnodes23(4)
  integer                        :: row, row2, ind, n3
  real(kind=8)                   :: aux, vc(4), dx(4), dy(4), npx, npy
  real(kind=8)                   :: inv4
  real(kind=8)                   :: cori_p, dparam, beta, gamma
#ifdef use_non_hydrostatic
  real(kind=8)                   :: dz(4), npz
#endif
  !
  inv4=0.25_8
  n3=myDim_nod3D+eDim_nod3D          
  do row=1, myDim_nod3d                 
     row2=row+n3                     
     uv_rhs(row)=0.
     uv_rhs(row2)=0.
#ifdef use_non_hydrostatic
     uv_rhs(row2+n3)=0.               
#endif
  enddo
  !
  do elem=1, myDim_elem3d             
     elnodes=elem3D_nodes(:,elem)
     elnodes23=nod2D_corresp_to_nod3D(elnodes)
     dx=bafux_3D(:,elem)
     dy=bafuy_3D(:,elem)
#ifdef use_non_hydrostatic
     dz=bafuz_3D(:,elem)
#endif
     if(use_cori_semi) then
        aux=voltetra(elem)*inv4
        !coriolis parameter
        el2=elem2d_corresp_to_elem3d(elem)
        cori_p=coriolis_param_elem2d(el2)
        !terms for semi-implicit schemes
        dparam=dt_inv**2 + alpha_trapez**2*cori_p**2
        beta=dt_inv/dparam
        gamma=cori_p*alpha_trapez/dparam
     else
        aux=voltetra(elem)*inv4*dt
     endif
     !ssh/non-hydrostatic pressure:
#ifndef use_non_hydrostatic
#ifdef use_semiimplicit_scheme
     vc=ssh0(elnodes23)*(gamma_stab-1.0+theta_ssh)-ssh(elnodes23)*theta_ssh
#else
     vc=ssh0(elnodes23)*gamma_stab-ssh(elnodes23)
#endif
     npx=g*sum(dx*vc)
     npy=g*sum(dy*vc)
#else
     vc=g*ssh0(elnodes23)*gamma_stab+nhp0(elnodes)*gamma_stab_nh - &
          (g*ssh(elnodes23)+nhp(elnodes))
     npx=sum(dx*vc)
     npy=sum(dy*vc)
     npz=sum(dz*vc) 
#endif
     !
     !assembling over nodes
     if(use_cori_semi) then
        uv_rhs(elnodes)=uv_rhs(elnodes)+(beta*npx+gamma*npy)*aux
        uv_rhs(elnodes+n3)=uv_rhs(elnodes+n3)+(beta*npy-gamma*npx)*aux 
#ifdef use_non_hydrostatic   
        uv_rhs(elnodes+2*n3)=uv_rhs(elnodes+2*n3)+npz*aux*dt            
#endif
     else
        uv_rhs(elnodes)=uv_rhs(elnodes)+npx*aux
        uv_rhs(elnodes+n3)=uv_rhs(elnodes+n3)+npy*aux                 
#ifdef use_non_hydrostatic   
        uv_rhs(elnodes+2*n3)=uv_rhs(elnodes+2*n3)+npz*aux              
#endif  
     endif
  end do
  !
  ! boundary conditions
  do row=1,myDim_nod3d                      
     row2=row+n3                           
     ind=index_nod3d(row)
     if((ind==11).or.(ind==21).or.(ind==31)) then
        uv_rhs(row)=0.0
        uv_rhs(row2)=0.0 
#ifdef use_non_hydrostatic
        uv_rhs(row2+n3)=0.0                
#endif 
     end if
  end do
end subroutine velocity_rhs_update
!
!----------------------------------------------------------------------------
!
subroutine impl_vertvisc
  ! apply separated implicit vertical viscosity; lumped mass matrix is used
  !
  ! Coded by Qiang Wang and Sergey Danilov
  ! Reviewed by ??
  !-------------------------------------------------------------
  
  use o_param
  use o_array
  use o_mesh
  use o_elements
  use g_config
  use g_parfe
  implicit none

  integer            :: m, n2, nlay, k, i, q, p, cnt
  integer            :: row, elem, elnodes(4), n_el
  integer            :: nodup, nodlo
  real(kind=8)       :: vol, vol_visc, inv4, Av_el
  real(kind=8)       :: a(max_num_layers), b(max_num_layers), c(max_num_layers)
  real(kind=8)       :: ru(max_num_layers), rv(max_num_layers)

  inv4=0.25_8

  do n2=1,myDim_nod2d  !only myDim_nod2d required

     !if(mod(index_nod3d(n2),10)==1) cycle  !no-slip bc.	

     nlay=num_layers_below_nod2d(n2)+1

     ! assemble three diagonal matrix and rhs
     do k=1,nlay
        a(k)=0.
        b(k)=0.
        c(k)=0.
        ru(k)=0.
        rv(k)=0.

        row=nod3d_below_nod2d(k,n2)

        if(mod(index_nod3d(row),10)==1) then
           b(k)=1.0
           cycle !no-slip bc.
        end if

	nodup=0
	nodlo=0
        if(k>1) nodup=nod3d_below_nod2d(k-1,n2)
        if(k<nlay) nodlo=nod3d_below_nod2d(k+1,n2)
        do i=1,nod_in_elem3D(row)%nmb
           elem=nod_in_elem3D(row)%addresses(i)
           elnodes=elem3D_nodes(:,elem)
           vol_visc=voltetra(elem)
#ifdef use_fullfreesurf
           n_el=map_elem(elem)
           if(n_el==0) then
              vol=voltetra(elem)
           else
              vol=voltetra_new(n_el)
           end if
#else
           vol=voltetra(elem)
#endif

           do q=1,4
              if(elnodes(q)==row) then
                 p=q
                 exit
              end if
           end do

           cnt=0
           Av_el=0.0

           !first entry
           if(k>1) then
              do q=1,4
                 if(elnodes(q)==nodup) then
                    Av_el=Av(nodup)    !Av(nodup) is the Av for this segment
                    if(trim(mix_scheme)=='MY2p5') then
                       Av_el=Av_el+Av0
                    end if
                    a(k)=a(k) + bafuz_3d(p,elem) * bafuz_3d(q,elem) * Av_el * vol_visc
                    cnt=1
                    exit
                 end if
              end do
           end if

           !third entry
           if(k<nlay) then
              do q=1,4
                 if(elnodes(q)==nodlo) then
                    Av_el=Av(row)    !Av(row) is the Av for this segment
                    if(trim(mix_scheme)=='MY2p5') then
                       Av_el=Av_el+Av0
                    end if
                    c(k)=c(k) + bafuz_3d(p,elem) * bafuz_3d(q,elem) * Av_el * vol_visc
                    cnt=1
                    exit
                 end if
              end do
           end if

           !second entry
           b(k)=b(k) + bafuz_3d(p,elem) * bafuz_3d(p,elem) * Av_el * vol_visc

           !second entry, mass matrix part
           b(k)=b(k) + vol*inv4*dt_inv
           !rhs
           ru(k)=ru(k) + uf(row)*vol*inv4*dt_inv
           rv(k)=rv(k) + uf(myDim_nod3d+eDim_nod3D+row)*vol*inv4*dt_inv 

        end do  !i
     end do  !k

     ! surface and bottom boundary conditions
     ru(1)=ru(1) + uv_sfc_force(n2,1)
     rv(1)=rv(1) + uv_sfc_force(n2,2)
     ru(nlay)=ru(nlay) + uv_bott_force(n2,1)
     rv(nlay)=rv(nlay) + uv_bott_force(n2,2)

     ! the sweep algorithm
     ! forward step
     ! prepare coefficients 
     c(1)=-c(1)/b(1)
     ru(1)=ru(1)/b(1)
     rv(1)=rv(1)/b(1)
     do k=2,nlay
        b(k)=1.0/(a(k)*c(k-1)+b(k))
        c(k)=-c(k)*b(k)
        ru(k)=(ru(k)-a(k)*ru(k-1))*b(k)
        rv(k)=(rv(k)-a(k)*rv(k-1))*b(k)
     end do
     ! backward sweep
     do k=nlay-1,1,-1
        ru(k)=ru(k)+c(k)*ru(k+1)
        rv(k)=rv(k)+c(k)*rv(k+1)
     end do

     ! update the velocity
     do k=1,nlay
        row=nod3d_below_nod2d(k,n2)
        uf(row)=ru(k)
        uf(row+myDim_nod3d+eDim_nod3D)=rv(k)   
     end do

  end do  !m

  call com_3D(uf(1:myDim_nod3d+eDim_nod3D))     
  call com_3D(uf(1+myDim_nod3d+eDim_nod3D:2*(myDim_nod3d+eDim_nod3D))) 

end subroutine impl_vertvisc
!
!============================================================================
!
subroutine compute_ssh_rhs
  ! rhs for ssh
  !
  ! Coded by Qiang Wang
  ! Reviewed by ??
  !-------------------------------------------------------------
  
  use o_param
  use o_array
  use o_mesh
  use o_elements
  use g_config
  use g_parfe
  implicit none
  !
  integer             :: elem, el2, elnodes(4), elnodes2(3)
  integer             :: q, row, m, mn(3), n3
  real(kind=8)        :: vol, dx(3), dy(3), us, vs
  real(kind=8)        :: inv4, inv12, aux, aux1, aux2
  real(kind=8)        :: u_el(4), v_el(4), tri_u(3), tri_v(3)
  real(kind=8)        :: cori_p, dparam, beta, gamma
#ifdef use_fullfreesurf
  integer             :: n_el
#endif
  !
  inv4=0.25_8
  inv12=1.0_8/12.0_8
  n3=myDim_nod3D+eDim_nod3D       
  do row=1,myDim_nod2d            
     ssh_rhs(row)=0.
  enddo

  ! divergence contribution
  do elem=1,myDim_elem3d        
     el2=elem2d_corresp_to_elem3d(elem)
     elnodes2=elem2d_nodes(:,el2)
     elnodes=elem3D_nodes(:,elem)
     !elnodes23=nod2D_corresp_to_nod3D(elnodes)
     dx=bafux_2D(:,el2)
     dy=bafuy_2D(:,el2)
#ifdef use_fullfreesurf
     n_el=map_elem(elem)
     if(n_el==0) then
        vol=voltetra(elem)
     else
        vol=voltetra_new(n_el)
     end if
#else
     vol=voltetra(elem)
#endif    
#ifdef use_semiimplicit_scheme
     u_el=uf(elnodes)*theta_vel+uf0(elnodes)*(1.0-theta_vel)
     v_el=uf(elnodes+n3)*theta_vel+uf0(elnodes+n3)*(1.0-theta_vel)   
#else
     u_el=uf(elnodes)
     v_el=uf(elnodes+n3)                                           
#endif
     !
     aux=g*(gamma_stab-1.0)
#ifdef use_semiimplicit_scheme
     aux=aux*theta_vel
#endif
     aux1=sum(dx*ssh0(elnodes2))*aux
     aux2=sum(dy*ssh0(elnodes2))*aux
     if(use_cori_semi) then
        cori_p=coriolis_param_elem2d(el2)
        dparam=dt_inv**2 + alpha_trapez**2*cori_p**2
        beta=dt_inv/dparam
        gamma=cori_p*alpha_trapez/dparam
        us=sum(u_el)*inv4+beta*aux1+gamma*aux2
        vs=sum(v_el)*inv4+beta*aux2-gamma*aux1 
     else
        us=sum(u_el)*inv4+dt*aux1
        vs=sum(v_el)*inv4+dt*aux2 
     endif
     !
     do q=1,3
	row=elnodes2(q)
	ssh_rhs(row)=ssh_rhs(row)+(dx(q)*us+dy(q)*vs)*vol
     end do
  end do

#ifdef use_fullfreesurf
  ! P-E contribution
  do el2=1,myDim_elem2d                
     elnodes2=elem2D_nodes(:,el2)
     vs=sum(water_flux(elnodes2))  ! '+' is defined as upwards 
     vol=voltriangle(el2)*inv12
     do q=1, 3
        row=elnodes2(q)
        ssh_rhs(row)=ssh_rhs(row) - (vs+water_flux(row))*vol
     end do
  end do
#endif

  ! open boundary contribution
#ifdef use_opbnd_tide
  if(trim(tide_opbnd_type)=='ssh') then
     do q=1, nmbr_opbnd_t2d
        row=opbnd_n2d(q)
        ssh_rhs(row)=opbnd_z_tide(q)-opbnd_z0_tide(q)  ! increments
     end do
  else
     do el2=1, nmbr_opbnd_tri
        elnodes2=opbnd_tri(el2,1:3)
        elnodes2=nod2d_corresp_to_nod3d(elnodes2)
        mn=mapping_opbnd_n2d(elnodes2)
        vol=opbnd_nv(el2,4)*inv12
	tri_u=opbnd_u_tide(mn)
 	tri_v=opbnd_v_tide(mn)
        tri_v=tri_u*opbnd_nv(el2,1)+tri_v*opbnd_nv(el2,2)
        if(trim(tide_opbnd_type)=='Flather') then
           tri_v=tri_v + sqrt(g/opbnd_dep(mn))*(ssh(elnodes2)-opbnd_z_tide(mn)) 
        end if
        vs=sum(tri_v)
        do q=1, 3
	   row=elnodes2(q)
           ssh_rhs(row)=ssh_rhs(row) - (vs+tri_v(q))*vol 
        end do
     end do
  end if
#endif
#ifdef use_opbnd_restoring
  do q=1, nmbr_opbnd_t2d
     row=opbnd_n2D(q)
     m=mapping_opbnd_n2d(row)
     ssh_rhs(row)=ssh_rhs(row)+opbnd_ssh_rhs(m)
  end do
#endif

end subroutine compute_ssh_rhs
!===================================================================	      

#ifdef use_non_hydrostatic
subroutine compute_nhp_rhs
  ! rhs for non-hydro.
  ! this is only a test version!
  !
  ! Coded by Qiang Wang
  ! Reviewed by ??
  !--------------------------------------------------------
  
  use o_param
  use o_array
  use o_mesh
  use o_elements
  use g_config
  use g_parfe
  implicit none
  !
  integer             :: elem, elem2, q, row, m, n3
  integer             :: elnodes(4), elnodes23(4), elnodes2(3)
  real(kind=8)        :: vol, dx(4), dy(4), dz(4), us, vs, ws
  real(kind=8)        :: inv4, inv12, vc(4), aux1, aux2
  real(kind=8)        :: cori_p, dparam, beta, gamma
#ifdef use_fullfreesurf
  integer             :: n_el
#endif
  !
  inv4=0.25_8
  inv12=1.0_8/12.0_8
  n3=myDim_nod3D+eDim_nod3D
  do row=1,myDim_nod3d              
     nhp_rhs(row)=0.
  enddo
  !
  do elem=1,myDim_elem3d          
     elnodes=elem3D_nodes(:,elem)
     elnodes23=nod2d_corresp_to_nod3d(elnodes)
#ifdef use_fullfreesurf
     n_el=map_elem(elem)
     if(n_el==0) then
        dx=bafux_3D(:,elem)
        dy=bafuy_3D(:,elem)
        dz=bafuz_3D(:,elem)
        vol=voltetra(elem)
     else
        dx=bafux_3D_new(:,n_el)
        dy=bafuy_3D_new(:,n_el)
        dz=bafuz_3D_new(:,n_el)
        vol=voltetra_new(n_el)
     end if
#else
     dx=bafux_3D(:,elem)
     dy=bafuy_3D(:,elem)
     dz=bafuz_3D(:,elem)
     vol=voltetra(elem)
#endif    
     !
     vc=nhp(elnodes)*gamma_stab_nh+g*ssh0(elnodes23)*gamma_stab-g*ssh(elnodes23)
     aux1=sum(dx*vc)
     aux2=sum(dy*vc)
     if(use_cori_semi) then
        elem2=elem2d_corresp_to_elem3d(elem)
        cori_p=coriolis_param_elem2d(elem2)
        dparam=dt_inv**2 + alpha_trapez**2*cori_p**2
        beta=dt_inv/dparam
        gamma=cori_p*alpha_trapez/dparam
        us=sum(uf(elnodes))*inv4 + beta*aux1+gamma*aux2
        vs=sum(uf(elnodes+n3))*inv4 + beta*aux2-gamma*aux1    
     else
        us=sum(uf(elnodes))*inv4 + dt*aux1
        vs=sum(uf(elnodes+n3))*inv4 + dt*aux2                   
     endif
     ws=sum(uf(elnodes+2*n3))*inv4 + dt*sum(dz*nhp(elnodes))*gamma_stab_nh 
     do q=1,4
	row=elnodes(q)
	nhp_rhs(row)=nhp_rhs(row)+(dx(q)*us+dy(q)*vs+dz(q)*ws)*vol
     end do
  end do
  !
  ! surface integral 
  do elem2=1, myDim_elem2d               
     vol=voltriangle(elem2)*inv12
     elnodes2=elem2D_nodes(:,elem2)
     us=sum(dssh(elnodes2))
     vs=sum(water_flux(elnodes2))
     do q=1,3
        row=elnodes2(q)
        nhp_rhs(row)  = nhp_rhs(row) - &
             vol*((us+dssh(row))*dt_inv+vs+water_flux(row))
     end do
  end do
end subroutine compute_nhp_rhs
#endif
!===================================================================
!
subroutine uv_solve
  ! Tricky solver for velocity
  ! Getting better efficiency by loosing accuracy!
  ! Better solution method in the future required.
  !
  ! Coded by Qiang Wang
  ! Reviewed by Sergey Danilov
  !---------------------------------------------------------

  use o_matrices
  use o_mesh
  use o_elements
  use o_param
  use o_array
  use g_config
  use g_parfe
  implicit none
  !
  integer                         :: n, m, cn, clo, clo2, location(100)
  integer                         :: row, row2, n3
  real(kind=8)                    :: u_rhs_new
#ifdef use_non_hydrostatic
  integer                         :: row3
  real(kind=8), allocatable       :: temp_array(:)
  allocate(temp_array(myDim_nod3d))
#endif
  !
  n3=myDim_nod3D+eDim_nod3D
  !the first approximation
  do row=1,myDim_nod3d                   
     row2=row+n3    	                  
     duf(row)=uv_rhs(row)/uv_lump(row)
     duf(row2)=uv_rhs(row2)/uv_lump(row)
#ifdef use_non_hydrostatic
     row3=row2+n3                         
     duf(row3)=uv_rhs(row3)/uv_lump(row)
#endif
  end do
  call com_3d(duf(1:n3))                   
  call com_3d(duf(1+n3:2*n3))             
#ifdef use_non_hydrostatic
  call com_3d(duf(1+2*n3:3*n3))           
#endif
  !

  !iterate

  do n=1, num_iter_solve-1
     do row=1,myDim_nod3d                 
 	row2=row+n3                     
        clo=uvstiff%rowptr(row)-uvstiff%rowptr(1)+1     
        clo2=uvstiff%rowptr(row+1)-uvstiff%rowptr(1)    
        cn=clo2-clo+1
        location(1:cn)=nghbr_nod3D(row)%addresses      
        u_rhs_new=uv_rhs(row) - sum(uvstiff%values(clo:clo2)*duf(location(1:cn)))
        dtracer(row,1)=duf(row)+u_rhs_new/uv_lump(row)
	u_rhs_new=uv_rhs(row2) - sum(uvstiff%values(clo:clo2)*duf(location(1:cn)+n3)) 
        dtracer(row,2)=duf(row2)+u_rhs_new/uv_lump(row)
#ifdef use_non_hydrostatic
        row3=row2+n3                                    
        u_rhs_new=uv_rhs(row3) - sum(uvstiff%values(clo:clo2)*duf(location(1:cn)+2*n3))
      	temp_array(m)=duf(row3)+u_rhs_new/uv_lump(row)
#endif
     end do
     do row=1,myDim_nod3d                  
	row2=row+n3                     
        duf(row)=dtracer(row,1)
        duf(row2)=dtracer(row,2)
#ifdef use_non_hydrostatic
        row3=row2+n3                     
        duf(row3)=temp_array(m)
#endif
     end do
     call com_3d(duf(1:n3))              
     call com_3d(duf(1+n3:2*n3))         
#ifdef use_non_hydrostatic
     call com_3d(duf(1+2*n3:3*n3))       
#endif

  end do
  !
#ifdef use_non_hydrostatic
  deallocate(temp_array)
#endif  
  !
end subroutine uv_solve
!
!========================================================================
!
#ifndef use_non_hydrostatic
subroutine compute_vvel_rhs
  ! Compute RHS in vertical velocity equation
  !
  ! Coded by Qiang Wang
  ! Reviewed by ??
  !-------------------------------------------------------------
  
  use o_mesh
  use o_elements
  use o_param
  use o_array
  use g_config
  use g_parfe
  implicit none
  !
  integer                 :: i, m, n, q, elem, elem2, row, n3
  integer                 :: elnodes(4), elnodes2(3), mn(3)
  real(kind=8)            :: val, ux, vy, inv4, inv12, dx(4), dy(4)
  real(kind=8)            :: vc(3), aux1, aux2, tri_u(3), tri_v(3)
  real(kind=8)            :: cori_p, dparam, beta, gamma
#ifdef use_fullfreesurf
  integer                 :: n_el
#endif  

  inv4=0.25_8
  inv12=1.0_8/12.0_8
  n3=myDim_nod3D+eDim_nod3D
  do row=1,myDim_nod3D              
     wrhs(row)=0.0
  enddo

  ! compute -\nabla\tilde\phi v* d\Omega   
  do elem=1, myDim_elem3d             
     elem2=elem2d_corresp_to_elem3d(elem)
     elnodes=elem3D_nodes(:,elem)
     elnodes2=elem2D_nodes(:,elem2)
     ux=inv4*sum(uf0(elnodes))
     vy=inv4*sum(uf0(elnodes+n3))    
#ifdef use_fullfreesurf
     n_el=map_elem(elem)
     if(n_el==0) then
        dx=bafux_3D(:,elem)
        dy=bafuy_3D(:,elem)
        val=voltetra(elem)
     else
        dx=bafux_3D_new(:,n_el)
        dy=bafuy_3D_new(:,n_el)
        val=voltetra_new(n_el)
     endif
#else
     dx=bafux_3D(:,elem)
     dy=bafuy_3D(:,elem)
     val=voltetra(elem)
#endif
     !
#ifdef use_semiimplicit_scheme
     vc=-g*(theta_ssh*ssh(elnodes2)-(gamma_stab-1.0+theta_ssh)*ssh0(elnodes2))
#else
     vc=-g*(ssh(elnodes2)-gamma_stab*ssh0(elnodes2))
#endif
     aux1=sum(vc*bafux_2D(:,elem2))
     aux2=sum(vc*bafuy_2D(:,elem2))
     if(use_cori_semi) then
        cori_p=coriolis_param_elem2d(elem2)
        dparam=dt_inv**2 + alpha_trapez**2*cori_p**2
        beta=dt_inv/dparam
        gamma=cori_p*alpha_trapez/dparam
        ux=ux+beta*aux1+gamma*aux2
        vy=vy+beta*aux2-gamma*aux1
     else
        ux=ux+dt*aux1
        vy=vy+dt*aux2
     endif
     !
     do q=1,4
        row=elnodes(q)
        wrhs(row) = wrhs(row)  - val*(dx(q)*ux+dy(q)*vy)
     end do
  end do

  ! sea surface boundary condition/integral
  do elem2=1, myDim_elem2d                   
     val=voltriangle(elem2)*inv12
     elnodes2=elem2D_nodes(:,elem2)
     ux=sum(dssh(elnodes2))
#ifdef use_fullfreesurf
     vy=sum(water_flux(elnodes2))
#endif
     do q=1,3
        row=elnodes2(q)
        n=nod3d_below_nod2d(1,row)
        wrhs(n) = wrhs(n) + &
             val*(ux+dssh(row))*dt_inv
#ifdef use_fullfreesurf
        wrhs(n) = wrhs(n) + val*(vy+water_flux(row))
#endif   
     end do
  end do


!!$  ! open boundary contribution
!!$#if defined(use_opbnd_tide) || defined(use_opbnd_restoring)
!!$  do elem2=1, nmbr_opbnd_tri
!!$     elnodes2=opbnd_tri(elem2,1:3)  !contains the 3d nodes on the opbnd.
!!$     val=opbnd_nv(elem2,4)*inv12
!!$     tri_u=uf(elnodes2)
!!$     tri_v=uf(elnodes2+nod3d)
!!$     tri_v=tri_u*opbnd_nv(elem2,1)+tri_v*opbnd_nv(elem2,2)
!!$     vy=sum(tri_v)
!!$     do q=1,3
!!$        row=elnodes2(q)
!!$        wrhs(row)=wrhs(row) + (vy+tri_v(q))*val 
!!$     end do
!!$  end do
!!$#endif


  ! Check solvability
  do row=1, myDim_nod2d     
     val=0.0
     do n=1,num_layers_below_nod2D(row)+1
        m = nod3D_below_nod2D(n,row)
        val=val+wrhs(m) 
     end do
     ! As a rule, val is a factor 10^(-4) smaller
     ! than wrhs(row), and is even smaller compared to 
     ! wrhs at nodes located deeper. Yet this accuracy
     ! is insufficient and thus solvability is enforced.
     val=val/real(num_layers_below_nod2D(row)+1)
     do n=1,num_layers_below_nod2D(row)+1
        m = nod3D_below_nod2D(n,row)
        wrhs(m)= wrhs(m)-val
     end do
  end do

end subroutine compute_vvel_rhs
#endif
!
!========================================================================
!
subroutine solve_wpot
  ! Solve vertical velocity potential using sweep algorithm
  ! The tetra. discretization allow to solve the equ. column by column.
  !
  ! Coded by Sergey Danilov and Qiang Wang
  ! Reviewed by ??
  !-------------------------------------------------------------
  
  use o_param
  use o_array
  use o_mesh
  use o_matrices
  use g_config
  use g_parfe
  implicit none

  integer            :: n2, nlay, k, row
  real(kind=8)       :: a(max_num_layers), b(max_num_layers)
  real(kind=8)       :: c(max_num_layers), rw(max_num_layers)

  do n2=1,myDim_nod2d                      
     nlay=num_layers_below_nod2d(n2)+1

     ! matrix entries
     a(1:nlay)=wpot_matrix(1,1:nlay,n2) 
     b(1:nlay)=wpot_matrix(2,1:nlay,n2)  
     c(1:nlay)=wpot_matrix(3,1:nlay,n2)  

     !rhs
     do k=1,nlay
        row=nod3d_below_nod2d(k,n2)
        rw(k)=wrhs(row)
     end do

     ! due to increased order of operator in this equation, we should/can 
     ! specify the boundary (surface) value. Otherwise, the full matrix is 
     ! singular and not solvable.
     rw(1)=0.0

     if(nlay>2) then
        ! the sweep algorithm
        ! forward step
        ! prepare coefficients 
        c(2)=-c(2)/b(2)
        rw(2)=rw(2)/b(2)
        do k=3,nlay
           b(k)=a(k)*c(k-1)+b(k)
           c(k)=-c(k)/b(k)
           rw(k)=(rw(k)-a(k)*rw(k-1))/b(k)
        end do
        ! backward sweep
        do k=nlay-1,2,-1
           rw(k)=rw(k)+c(k)*rw(k+1)
        end do
     else
        rw(2)=rw(2)/b(2)
     end if

     ! put results to w array
     do k=1,nlay
        row=nod3d_below_nod2d(k,n2)
        w(row)=rw(k)
     end do

  end do  !m

  call com_3D(w)

end subroutine solve_wpot
!
!========================================================================
!
#ifndef use_non_hydrostatic
subroutine vvel_nodes
  ! 
  ! Coded by Qiang Wang
  ! Reviewed by ??
  !---------------------------------------------------------
  
  use o_param
  use o_array
  use o_mesh
  use o_elements
  use g_config
  use g_parfe
  implicit none
  !
  integer   :: col, elem, n, m, elnodes(4)
  real      :: vol, dz(4), elemw
  !
  !  Computes approximate vertical velocities at nodes as mean over neighbour
  !  elements

  do col=1,myDim_nod3d              
     wrhs(col)=0.
     vol=0.0
     do n=1,nod_in_elem3D(col)%nmb
        elem = nod_in_elem3D(col)%addresses(n)
        elnodes=elem3D_nodes(:,elem)
        dz=bafuz_3D(:,elem)
        vol=vol+voltetra(elem)
        elemw=sum(dz*w(elnodes))
        wrhs(col)=wrhs(col)+elemw*voltetra(elem) 
     end do
     wrhs(col)=wrhs(col)/vol
  end do
end subroutine vvel_nodes
#endif        
!============================================================================   
!
subroutine biharmonic_viscosity
  ! biharmonic viscosity
  !
  ! Coded by Sergey Danilov, 2006
  ! Reviewed by Qiang Wang, 2006
  ! Adjusted by Qiang Wang into the tetra. version, 2008
  ! Modified by Qiang Wang to allow for Smagorinsky viscosity, 2009
  !--------------------------------------------------------------

  use o_mesh
  use o_elements
  use o_param 
  use o_array
  use o_matrices
  use g_config
  use g_parfe
  implicit none
  integer                      :: m, q, row, row2, elnodes(4)
  integer                      :: elem, elem2, ind, lay, n3, elem_type
  integer                      :: n_el, bott_ind
  real(kind=8)                 :: dx(4), dy(4), dz(4), vol
  real(kind=8)                 :: vtr, Ah, Ah_sma, d1, d2
  real(kind=8)                 :: u_el(4), v_el(4), inv4 
  real(kind=8)                 :: udx, udy, udz, aux1
  real(kind=8)                 :: vdx, vdy, vdz, aux2
  real(kind=8)                 :: S(3), rotate_coe, temp_coe, temp_coe2
#ifdef use_non_hydrostatic
  integer                      :: row3
  real(kind=8)                 :: w_el(4), wdx, wdy, wdz, aux3
#endif
  real(kind=8)                 :: cori_p, dparam, beta, gamma
  !
  inv4=0.25_8
  n3=myDim_nod3D+eDim_nod3D
  do row=1, myDim_nod3d        
     row2=row+n3                
     duf(row)=0.
     duf(row2)=0.
#ifdef use_non_hydrostatic
     row3=row2+n3               
     duf(row3)=0.
#endif
  enddo
  !
  ! ==== store laplacian of velocity, with terms due to
  !      differentiation of metrics ignored
  do elem=1,myDim_elem3d        
     elnodes=elem3D_nodes(:,elem)
     elem2=elem2d_corresp_to_elem3d(elem)
     elem_type=grid_type_elem2d(elem2)
     vol=voltetra(elem)
     dx=bafux_3d(:, elem)
     dy=bafuy_3d(:, elem)
     u_el=uf(elnodes)
     v_el=uf(n3+elnodes)       
     udx=sum(dx*u_el)
     udy=sum(dy*u_el)
     vdx=sum(dx*v_el)
     vdy=sum(dy*v_el)
#ifdef use_non_hydrostatic
     w_el=uf(n3*2+elnodes)    
     wdx=sum(dx*w_el)
     wdy=sum(dy*w_el)
#endif
     if (elem_type==1) then
        dz=bafuz_3d(:,elem)
        udz=sum(dz*u_el)
        vdz=sum(dz*v_el)
#ifdef use_non_hydrostatic
        wdz=sum(dz*w_el)
#endif
        lay=elem3d_layer(elem)
        S(1)=grid_slope(1,lay,elem2)
        S(2)=grid_slope(2,lay,elem2)
        aux1=S(1)**2+S(2)**2
        S(3)=sqrt(aux1)
        rotate_coe=1.0/(1.0+aux1) 
     end if
     !
     do q=1,4 
        row=elnodes(q)
        row2=row+n3          
	ind=index_nod3D(row)
        if((ind==11).or.(ind==21).or.(ind==31)) cycle
        !
        if(elem_type==1) then  !along-sigma viscosity
           !diagonal part1 (lateral)
	   temp_coe=rotate_coe*vol
           aux1=-(dx(q)*udx*(1.0+S(2)*S(2))+ &
                dy(q)*udy*(1.0+S(1)*S(1)))*temp_coe
           aux2=-(dx(q)*vdx*(1.0+S(2)*S(2))+ &
                dy(q)*vdy*(1.0+S(1)*S(1)))*temp_coe
#ifdef use_non_hydrostatic
           aux3=-(dx(q)*wdx*(1.0+S(2)*S(2))+ &
                dy(q)*wdy*(1.0+S(1)*S(1)))*temp_coe
#endif
           !diagonal part2 (cross slope)
	   temp_coe2=S(3)*S(3)*temp_coe
           aux1=aux1-dz(q)*udz*temp_coe2
           aux2=aux2-dz(q)*vdz*temp_coe2
#ifdef use_non_hydrostatic
           aux3=aux3-dz(q)*wdz*temp_coe2
#endif
           !off diagonal part1 (lateral) --> (1,3),(2,3)
           aux1=aux1-(S(1)*dx(q)+S(2)*dy(q))*udz*temp_coe
           aux2=aux2-(S(1)*dx(q)+S(2)*dy(q))*vdz*temp_coe
#ifdef use_non_hydrostatic
           aux3=aux3-(S(1)*dx(q)+S(2)*dy(q))*wdz*temp_coe
#endif
           !off diagonal part1 (lateral)--> (1,2),(2,1), '- -' -> '+'
	   temp_coe2=S(1)*S(2)*temp_coe
           aux1=aux1 + (dx(q)*udy+dy(q)*udx)*temp_coe2
           aux2=aux2 + (dx(q)*vdy+dy(q)*vdx)*temp_coe2
#ifdef use_non_hydrostatic
           aux3=aux3 + (dx(q)*wdy+dy(q)*wdx)*temp_coe2
#endif
           !off diagonal part2 (cross slope) --> (3,1),(3,2)
           aux1=aux1-(S(1)*udx+S(2)*udy)*dz(q)*temp_coe
           aux2=aux2-(S(1)*vdx+S(2)*vdy)*dz(q)*temp_coe 
#ifdef use_non_hydrostatic
           aux3=aux3-(S(1)*wdx+S(2)*wdy)*dz(q)*temp_coe 
#endif
        else
           aux1=-(dx(q)*udx + dy(q)*udy)*vol
           aux2=-(dx(q)*vdx + dy(q)*vdy)*vol
#ifdef use_non_hydrostatic
           aux3=-(dx(q)*wdx + dy(q)*wdy)*vol
#endif
        endif
        !
        duf(row)=duf(row) + aux1
        duf(row2)=duf(row2) +aux2
#ifdef use_non_hydrostatic
        row3=row2+n3                         
        duf(row3)=duf(row3) + aux3
#endif
     end do
  end do
  !
  ! ==== Return to physical space ===============
  do row=1,myDim_nod3d                           
     row2=row+n3                                 
     ind=index_nod3D(row)
     if((ind==11).or.(ind==21).or.(ind==31)) cycle
     duf(row)=duf(row)/uv_lump(row)  
     duf(row2)=duf(row2)/uv_lump(row)
#ifdef use_non_hydrostatic
     row3=row2+n3                                 
     duf(row3)=duf(row3)/uv_lump(row)
#endif
  end do
  !
  ! ==== Communication ==========================
  call com_3d(duf(1:n3))                          
  call com_3d(duf(1+n3:2*n3))                    
#ifdef use_non_hydrostatic
  call com_3d(duf(1+2*n3:3*n3))                  
#endif
  !
  ! ==== Compute Laplacian of Laplacian =========
  do elem=1,myDim_elem3d                     
     elnodes=elem3D_nodes(:,elem)
     elem2=elem2d_corresp_to_elem3d(elem)
     elem_type=grid_type_elem2d(elem2)
     if(use_cori_semi) then
        cori_p=coriolis_param_elem2D(elem2)
        dparam=dt_inv**2 + alpha_trapez**2*cori_p**2
        beta=dt_inv/dparam
        gamma=cori_p*alpha_trapez/dparam
        vol=voltetra(elem)
     else
        vol=voltetra(elem)*dt
     endif
     dx=bafux_3d(:, elem)
     dy=bafuy_3d(:, elem)
     u_el=duf(elnodes)
     v_el=duf(elnodes+n3)                       
     udx=sum(dx*u_el)
     udy=sum(dy*u_el)
     vdx=sum(dx*v_el)
     vdy=sum(dy*v_el)
#ifdef use_non_hydrostatic
     w_el=duf(n3*2+elnodes)                    
     wdx=sum(dx*w_el)
     wdy=sum(dy*w_el)
#endif

     !biharmonic viscosity 
     vtr=voltriangle(elem2)
     Ah=-Ahb0   
     if(scale_mixing_h) Ah=Ah*(sqrt(vtr/scalevol))**3

     !special: test increasing boundary Ah, Qiang
     !if(any(mod(index_nod3d(elnodes),10)==1)) then
     !   Ah=Ah*3.0
     !elseif(increase_equ_zonal_visc .and. any(abs(geolat(elnodes))<10.*rad)) then
     !   Ah=Ah*fac_visc_increase
     !end if
     
     if(increase_equ_zonal_visc .and. any(abs(geolat(elnodes))<10.*rad)) then
        Ah=Ah*fac_visc_increase
     end if

     if(smagorinsky_visc) then
	d1=sum(dx*uf(elnodes))-sum(dy*uf(elnodes+n3))  
        d2=sum(dy*uf(elnodes))+sum(dx*uf(elnodes+n3))  
        Ah_sma=sqrt(d1*d1+d2*d2)*vtr                   !laplacian viscosity; set C/pi to 1
        Ah_sma=max(Ah_sma, 40.0/1.0e8*vtr)             !limit from below: 20m^2/s on 10km
	if(Ah_sma*dt/vtr > 0.01) Ah_sma=0.01*vtr/dt    !limit from above
        Ah_sma=-Ah_sma*vtr/8.0                         !biharmonic viscosity
        Ah=Ah+Ah_sma
     endif

     if (elem_type==1) then
        dz=bafuz_3d(:,elem)
        udz=sum(dz*u_el)
        vdz=sum(dz*v_el)
#ifdef use_non_hydrostatic
        wdz=sum(dz*w_el)
#endif
        lay=elem3d_layer(elem)
        S(1)=grid_slope(1,lay,elem2)
        S(2)=grid_slope(2,lay,elem2)
        aux1=S(1)**2+S(2)**2
        S(3)=sqrt(aux1)
        rotate_coe=1.0/(1.0+aux1) 
     end if
     !
     !
     do q=1, 4 
        row=elnodes(q)
        row2=row+n3                          
        !
        if(elem_type==1) then  !along-sigma viscosity
           !diagonal part1 (lateral)
	   temp_coe=Ah*rotate_coe*vol
           aux1=-(dx(q)*udx*(1.0+S(2)*S(2))+ &
                dy(q)*udy*(1.0+S(1)*S(1)))*temp_coe
           aux2=-(dx(q)*vdx*(1.0+S(2)*S(2))+ &
                dy(q)*vdy*(1.0+S(1)*S(1)))*temp_coe
#ifdef use_non_hydrostatic
           aux3=-(dx(q)*wdx*(1.0+S(2)*S(2))+ &
                dy(q)*wdy*(1.0+S(1)*S(1)))*temp_coe
#endif
           !diagonal part2 (cross slope)
	   temp_coe2=S(3)*S(3)*temp_coe
           aux1=aux1-dz(q)*udz*temp_coe2
           aux2=aux2-dz(q)*vdz*temp_coe2
#ifdef use_non_hydrostatic
           aux3=aux3-dz(q)*wdz*temp_coe2
#endif
           !off diagonal part1 (lateral) --> (1,3),(2,3)
           aux1=aux1-(S(1)*dx(q)+S(2)*dy(q))*udz*temp_coe
           aux2=aux2-(S(1)*dx(q)+S(2)*dy(q))*vdz*temp_coe
#ifdef use_non_hydrostatic
           aux3=aux3-(S(1)*dx(q)+S(2)*dy(q))*wdz*temp_coe
#endif
           !off diagonal part1 (lateral)--> (1,2),(2,1), '- -' -> '+'
	   temp_coe2=S(1)*S(2)*temp_coe
           aux1=aux1 + (dx(q)*udy+dy(q)*udx)*temp_coe2
           aux2=aux2 + (dx(q)*vdy+dy(q)*vdx)*temp_coe2
#ifdef use_non_hydrostatic
           aux3=aux3 + (dx(q)*wdy+dy(q)*wdx)*temp_coe2
#endif
           !off diagonal part2 (cross slope) --> (3,1),(3,2)
           aux1=aux1-(S(1)*udx+S(2)*udy)*dz(q)*temp_coe
           aux2=aux2-(S(1)*vdx+S(2)*vdy)*dz(q)*temp_coe 
#ifdef use_non_hydrostatic
           aux3=aux3-(S(1)*wdx+S(2)*wdy)*dz(q)*temp_coe 
#endif
        else
           aux1=-Ah*(dx(q)*udx + dy(q)*udy)*vol
           aux2=-Ah*(dx(q)*vdx + dy(q)*vdy)*vol 
#ifdef use_non_hydrostatic
           aux3=-Ah*(dx(q)*wdx + dy(q)*wdy)*vol
#endif
        endif
        !
        if(use_cori_semi) then
           uv_rhs(row)=uv_rhs(row)+aux1*beta+aux2*gamma
           uv_rhs(row2)=uv_rhs(row2)+aux2*beta-aux1*gamma
#ifdef use_non_hydrostatic
           row3=row2+n3                            
           uv_rhs(row3)=uv_rhs(row3)+aux3*dt
#endif
        else
           uv_rhs(row)=uv_rhs(row)+aux1
           uv_rhs(row2)=uv_rhs(row2)+aux2
#ifdef use_non_hydrostatic
           row3=row2+n3                            
           uv_rhs(row3)=uv_rhs(row3)+aux3
#endif
        endif
     end do
  end do
  !
end subroutine biharmonic_viscosity
!
!--------------------------------------------------------------------------
