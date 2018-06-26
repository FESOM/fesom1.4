!1) Diagnose mixed layer depth 
!2) Diagnose SGS parameterization

subroutine compute_mixlay
  ! Diagnose mixed layer depth
  !
  ! Coded by Qiang Wang
  ! Reviewed by ??
  !-------------------------------------------------

  use o_MESH
  use o_PARAM
  use o_array
  use g_meanarrays
  use g_PARFE
  implicit none

  integer         :: m, n2, n3, k, row
  real(kind=8)    :: dens_surf, dens, dens_up, md, bf, bf_up 
  real(kind=8)    :: buoyancy_crit, sigma_theta_crit, smallvalue

  smallvalue=1.0e-20
  buoyancy_crit=0.0003
  sigma_theta_crit=0.125   !kg/m3, Levitus method

  !mixed layer depth  !Levitus method
  do n2=1,ToDim_nod2d
     row=nod3D_below_nod2D(2,n2) !assume the second layer is 10m. 
     call fcn_dens0(tracer(row,1),tracer(row,2),dens_surf)
     md=-10.0
     dens_up=dens_surf
     do k=3,num_layers_below_nod2d(n2)+1
        n3=nod3d_below_nod2d(k,n2)
        call fcn_dens0(tracer(n3,1),tracer(n3,2),dens)
        if(dens>=dens_surf+sigma_theta_crit) then
           md=md+(coord_nod3d(3,n3)-md)/(dens-dens_up+smallvalue)*(dens_surf+sigma_theta_crit-dens_up)
           exit
        else
           md=coord_nod3d(3,n3)
           dens_up=dens
        end if
     end do
     mlotst_mean%local_values(n2) = mlotst_mean%local_values(n2) + abs(md)
  end do


!!$  !mixed layer depth
!!$  do n2=1,ToDim_nod2d
!!$     row=nod3D_below_nod2D(1,n2)  
!!$     call fcn_dens0(tracer(row,1),tracer(row,2),dens_surf)
!!$     md=0.0
!!$     bf_up=0.0
!!$     do k=2,num_layers_below_nod2d(n2)+1
!!$        n3=nod3d_below_nod2d(k,n2)
!!$        call fcn_dens0(tracer(n3,1),tracer(n3,2),dens)
!!$        bf=g*(dens-dens_surf)/dens
!!$        if(bf>=buoyancy_crit) then
!!$           md=md+(coord_nod3d(3,n3)-md)/(bf-bf_up+smallvalue)*(buoyancy_crit-bf_up)
!!$           exit
!!$        else
!!$           md=coord_nod3d(3,n3)
!!$           bf_up=bf
!!$        end if
!!$     end do
!!$     mixlay_dep_mean(n2)=mixlay_dep_mean(n2) + abs(md)     
!!$  end do


  !!mixed layer depth
  !do n2=1,ToDim_nod2d
  !   row=nod3D_below_nod2D(1,n2)           
  !   md=0.0
  !   bf_up=0.0
  !   do k=2,num_layers_below_nod2d(n2)+1
  !      n3=nod3d_below_nod2d(k,n2)
  !      call fcn_density(tracer(row,1),tracer(row,2),coord_nod3d(3,n3),dens_surf)
  !      bf=g*(density_insitu(n3)-dens_surf)/density_insitu(n3)
  !      if(bf>=buoyancy_crit) then
  !         md=md+(coord_nod3d(3,n3)-md)/(bf-bf_up+smallvalue)*(buoyancy_crit-bf_up)
  !         exit
  !      else
  !         md=coord_nod3d(3,n3)
  !         bf_up=bf
  !      end if
  !   end do
  !   mixlay_dep_mean(n2)=mixlay_dep_mean(n2) + abs(md)     
  !end do

end subroutine compute_mixlay
!
!--------------------------------------------------------------------------------------------
!
subroutine process_elem2node(id,array_3d)
  ! Average elementwise variables to nodes
  ! Used for output ocean SGS parameterization
  ! Called in write_means_part2
  !
  ! Coded by Qiang Wang
  ! Reviewed by ??
  !-------------------------------------------------

  use o_mesh
  use o_elements
  use o_array
  use g_meanarrays
  use g_parfe
  implicit none

  integer       :: id,  elem, elnodes(4), row
  real(kind=8)  :: array_3d(nod3D)
  real(kind=8), allocatable  :: aux(:)  

  allocate(aux(myDim_nod3D+eDim_nod3D)) 
  aux=0.0 

  if(id==1) then
     do elem=1,myDim_elem3d       
        elnodes=elem3d_nodes(:,elem)
        aux(elnodes)=aux(elnodes)+sgs_u(elem)*voltetra(elem) 
     end do
     do row=1,myDim_nod3d         
        aux(row)=aux(row)/wrhs(row)  !wrhs contains cluster volume  
     end do
     call broadcast3d(aux,array_3d) 
  elseif(id==2) then
     do elem=1,myDim_elem3d          
        elnodes=elem3d_nodes(:,elem)
        aux(elnodes)=aux(elnodes)+sgs_v(elem)*voltetra(elem) 
     end do
     do row=1,myDim_nod3d           
        aux(row)=aux(row)/wrhs(row)  !wrhs contains cluster volume  
     end do
     call broadcast3d(aux,array_3d) 
  elseif(id==3) then
     do elem=1,myDim_elem3d          

        elnodes=elem3d_nodes(:,elem)
        aux(elnodes)=aux(elnodes)+sgs_ut(elem)*voltetra(elem) 
     end do
     do row=1,myDim_nod3d           
        aux(row)=aux(row)/wrhs(row)  !wrhs contains cluster volume  
     end do
     call broadcast3d(aux,array_3d)  
  elseif(id==4) then
     do elem=1,myDim_elem3d         
        elnodes=elem3d_nodes(:,elem)
        aux(elnodes)=aux(elnodes)+sgs_vt(elem)*voltetra(elem)
     end do
     do row=1,myDim_nod3d          
        aux(row)=aux(row)/wrhs(row)  !wrhs contains cluster volume 
     end do
     call broadcast3d(aux,array_3d)  
  elseif(id==5) then
     do elem=1,myDim_elem3d         
        elnodes=elem3d_nodes(:,elem)
        aux(elnodes)=aux(elnodes)+sgs_us(elem)*voltetra(elem) 
     end do
     do row=1,myDim_nod3d          
        aux(row)=aux(row)/wrhs(row)  !wrhs contains cluster volume 
     end do
     call broadcast3d(aux,array_3d)  
  elseif(id==6) then
     do elem=1,myDim_elem3d         
        elnodes=elem3d_nodes(:,elem)
        aux(elnodes)=aux(elnodes)+sgs_vs(elem)*voltetra(elem) 
     end do
     do row=1,myDim_nod3d          
        aux(row)=aux(row)/wrhs(row)  !wrhs contains cluster volume
     end do
     call broadcast3d(aux,array_3d)  
  end if
  deallocate(aux)
end subroutine process_elem2node
