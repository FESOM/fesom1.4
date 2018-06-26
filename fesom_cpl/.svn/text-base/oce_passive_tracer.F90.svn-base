module o_passive_tracer_mod
  use o_param
  use o_array
  use o_mesh
  use o_elements
  use g_config
  use g_forcing_arrays
  use g_clock
  use g_parfe
  implicit none

  integer, allocatable, dimension(:)        :: index_passive_tracer
  integer, allocatable, dimension(:,:)      :: passive_tracer_loc_index
  real(kind=8), allocatable, dimension(:,:) :: ptr_sfc_force

contains


  subroutine passive_tracer_init

    integer              :: i, j, k, n, n3, row, fileID
    integer              :: n_loc, num_nod
    integer, allocatable :: temp_arr2d(:), nodes_release(:)
    character(1)         :: cptrind
    character(4)         :: tr_name
    character(2000)       :: file_name

    !--------------------------------------------------------------  
    ! find index
    allocate(index_passive_tracer(num_passive_tracer))
    do j=1, num_passive_tracer
       write(cptrind,'(i1)') j
       tr_name='ptr'//cptrind
       do i=1, num_tracer
          if(prog_tracer_name(i) == tr_name) index_passive_tracer(j)=i
       end do
    end do

    !--------------------------------------------------------------  
    ! initial values
    do j=1, num_age_tracer
       tracer(:,index_passive_tracer(j))=ptr_background_value
    end do

    !--------------------------------------------------------------
    ! in case that p.tr is restored in a region
    if(passive_tracer_restore) then

       ! set passive tracer location index: 1 at release, 0 otherwise

       allocate(passive_tracer_loc_index(ToDim_nod3d,num_passive_tracer))
       passive_tracer_loc_index=0

       allocate(temp_arr2d(nod2d))
       temp_arr2d=0
       do n=1, ToDim_nod2D
          temp_arr2d(myList_nod2D(n))=n
       end do

       do j=1, num_passive_tracer
          write(cptrind,'(i1)') j
          tr_name='ptr'//cptrind
          file_name=trim(meshpath)//'passive_tracer_restore_nodes_'//tr_name//'.out'
          fileID=160
          open(fileID, file=file_name)
          read(fileID,*) num_nod
          allocate(nodes_release(num_nod))
          read(fileID,*) nodes_release
          close(fileID)
          do n=1,num_nod
             n_loc=temp_arr2d(nodes_release(n))
             if(n_loc>0) then
                n_loc=nod3d_below_nod2d(1,n_loc)
                passive_tracer_loc_index(n_loc,j)=1
                tracer(n_loc,index_passive_tracer(j))=ptr_restore_value
             end if
          end do
          deallocate(nodes_release)
       end do

       deallocate(temp_arr2d)

       !--------------------------------------------------------------
       ! in case restore volume
       if(ptr_restore_in_volume) then
          do i=1,ToDim_nod2d
             row=nod3d_below_nod2d(1,i)
             do j=1, num_passive_tracer
                if(passive_tracer_loc_index(row,j)==1) then
                   do k=2,num_layers_below_nod2d(i)+1
                      n3=nod3d_below_nod2d(k,i)
                      passive_tracer_loc_index(n3,j)=1
                      tracer(n3,index_passive_tracer(j))=ptr_restore_value
                   end do
                end if
             end do
          end do
       end if

    end if

    !--------------------------------------------------------------
    ! in case that passive tracers enter through surface fluxes
    if(passive_tracer_flux) then
       allocate(ptr_sfc_force(ToDim_nod2d,num_passive_tracer))
    end if

  end subroutine passive_tracer_init
  !
  !-------------------------------------------------------------------------
  !
  subroutine ptr_sfc_bc

    integer      :: elem, j, elnodes(3), elnodes2(3)
    real(kind=8) :: auxf, entries(3)

    if(passive_tracer_flux) then

       ptr_sfc_force=0.0

       do elem=1,myDim_elem2d             
          elnodes2=elem2D_nodes(:,elem)
          elnodes=nod3D_below_nod2D(1,elnodes2)   
          auxf=voltriangle(elem)/12.0_8    
          do j=1,num_passive_tracer
             entries=-auxf*(tracer(elnodes,2)+tracer(elnodes,index_passive_tracer(j))) &
                  * runoff_landice(elnodes2)*landice_season(month)
             ptr_sfc_force(elnodes2,j)=ptr_sfc_force(elnodes2,j)+sum(entries)+entries
          end do
       end do

    end if

  end subroutine ptr_sfc_bc
  !
  !-------------------------------------------------------------------------
  !
  subroutine ptr_cutoff_restore

    integer   :: j, row

    if(passive_tracer_restore) then
       do j=1, num_passive_tracer
          do row=1,ToDim_nod3d
             if(passive_tracer_loc_index(row,j)==1) then
                tracer(row,index_passive_tracer(j))=ptr_restore_value
             end if
          end do
       end do
    endif

  end subroutine ptr_cutoff_restore


end module o_passive_tracer_mod
