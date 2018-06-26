
subroutine landice_water_init
  ! init land ice melting rate
  use o_mesh
  use o_elements
  use o_param
  use o_array
  use g_config
  use g_forcing_param
  use g_forcing_arrays
  use g_clock
  use g_parfe
  implicit none

  integer                     :: n, i, j, num_reg, num_nod, num_days
  integer                     :: n_loc, fileID
  integer, allocatable        :: temp_arr2d(:), nodes_in_region(:)
  real(kind=8)                :: vol, vol_sum, aux
  real(kind=8), allocatable   :: totalflux(:)
  character*300               :: file_name
  character                   :: c_reg_ind

  ! read the number of region and the total yearly discharge in each region:
  file_name=trim(meshpath)//'landice_yearly_mass_loss.out' 
  fileID=160
  open(fileID, file=file_name)
  read(fileID,*) num_reg
  allocate(totalflux(num_reg))
  read(fileID,*) totalflux      !unit: Gt/year
  close(fileID)

  allocate(temp_arr2d(nod2d))
  temp_arr2d=0
  do n=1, myDim_nod2D  !note: eDim_nod2D should not be included in this case
     temp_arr2d(myList_nod2D(n))=n
  end do

  ! yearly mean runoff
  do i=1,num_reg
     
     !read in nodes in the region
     write(c_reg_ind,'(i1)') i   !assume that num_reg is less than 10
     file_name=trim(meshpath)//'landice_nodes_in_region_'//c_reg_ind//'.out' 
     fileID=160
     open(fileID, file=file_name)
     read(fileID,*) num_nod
     allocate(nodes_in_region(num_nod))
     read(fileID,*) nodes_in_region
     close(fileID)

     vol=0.0
     do n=1,num_nod
        n_loc=temp_arr2d(nodes_in_region(n))
        if(n_loc>0) then
           do j=1,nod_in_elem2d(n_loc)%nmb
              vol=vol+voltriangle(nod_in_elem2d(n_loc)%addresses(j))
           end do
        end if
     end do

     vol_sum=0.0
     call MPI_Barrier(MPI_COMM_FESOM,MPIerr)
     call MPI_AllREDUCE(vol, vol_sum, 1, MPI_DOUBLE_PRECISION,MPI_SUM, &
          MPI_COMM_FESOM, MPIerr)
 
     vol_sum=vol_sum/3.0
     
     aux=totalflux(i)*1.0e9/real(ndpyr)/86400.0_8/vol_sum  !m/s

     do n=1,num_nod
        n_loc=temp_arr2d(nodes_in_region(n))
        if(n_loc>0) then
           runoff_landice(n_loc)=aux
        end if
     end do

     deallocate(nodes_in_region)
  end do

  call com_2D(runoff_landice)

  ! seasonality
  num_days=sum(num_day_in_month(fleapyear,landice_start_mon:landice_end_mon))
  aux=real(ndpyr)/real(num_days)
  landice_season=0.0
  landice_season(landice_start_mon:landice_end_mon)=aux
    
  deallocate(temp_arr2d, totalflux)

  if(mype==0) write(*,*) 'Land-ice melt water fluxes prepared.'

end subroutine landice_water_init
!
!----------------------------------------------------------------------------------
!
subroutine add_landice_water
  use o_array
  use g_forcing_arrays
  use g_clock
  use g_parfe
  implicit none

  water_flux=water_flux-runoff_landice*landice_season(month)
  !water_flux is positive for upward

end subroutine add_landice_water
