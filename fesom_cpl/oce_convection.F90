subroutine convect_adjust
  ! Increase vertical mixing in case of static instability.
  ! It is already considered (included) by PP scheme.
  ! Only required for other cases
  !
  ! Coded by Qiang Wang
  ! Reviewed by ??
  !---------------------------------------------------
  
  use o_MESH
  use o_ELEMENTS
  use o_PARAM
  use o_array
  use g_config
  use g_PARFE
  implicit none
  !
  integer     	:: col,n,j,node_low,node_up,node
  real(kind=8)	:: density_up, drho

  do col=1,toDim_nod2d    

     n=num_layers_below_nod2D(col)

     do j=1,n       
        node=nod3D_below_nod2D(j,col)
        node_up = node
        node_low = nod3D_below_nod2D(j+1,col)

        call fcn_density(tracer(node_up,1),tracer(node_up,2), &
             coord_nod3d(3,node_low),density_up)

        drho = density_up - density_insitu(node_low)

        if(drho>=0.0) then
           if(allow_convect_global) then
              Kv(node,1)=diff_conv_limit
              Av(node)=visc_conv_limit
           elseif(coord_nod2d(2,col)>0.0) then     
              Kv(node,1)=diff_conv_limit
              Av(node)=visc_conv_limit
           endif
        end if

     end do

  end do
end subroutine convect_adjust

