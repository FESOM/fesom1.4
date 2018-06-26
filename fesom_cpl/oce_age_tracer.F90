module o_age_tracer_mod
  use o_param
  use o_array
  use o_mesh
  use o_elements
  use g_config
  use g_clock
  use g_parfe
  implicit none

  integer, allocatable, dimension(:)   :: index_age_tracer  
  integer, allocatable, dimension(:,:) :: age_tracer_loc_index

contains


  subroutine age_tracer_init

    integer              :: i, j, k, n, n3, row, fileID
    integer              :: n_loc, num_nod
    integer, allocatable :: temp_arr2d(:), nodes_release(:)
    character(1)         :: cageind
    character(4)         :: tr_name
    character(2000)       :: file_name

    !--------------------------------------------------------------  
    ! find index
    allocate(index_age_tracer(num_age_tracer))
    do j=1, num_age_tracer
       write(cageind,'(i1)') j
       tr_name='age'//cageind
       do i=1, num_tracer
          if(prog_tracer_name(i) == tr_name) index_age_tracer(j)=i
       end do
    end do

    !--------------------------------------------------------------  
    ! set initial values
    do j=1, num_age_tracer
       tracer(:,index_age_tracer(j))=0.0
    end do

    !-------------------------------------------------------------- 
    ! restore time scale at the release region
    if(zero_age_at_release) age_tracer_restore_time=dt

    !--------------------------------------------------------------    	
    ! set age tracer location index: 1 at release, 0 otherwise

    allocate(age_tracer_loc_index(ToDim_nod3d,num_age_tracer))
    age_tracer_loc_index=0

    allocate(temp_arr2d(nod2d))
    temp_arr2d=0
    do n=1, ToDim_nod2D
       temp_arr2d(myList_nod2D(n))=n
    end do

    do j=1, num_age_tracer
       write(cageind,'(i1)') j
       tr_name='age'//cageind
       file_name=trim(meshpath)//'age_tracer_release_nodes_'//tr_name//'.out'
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
             age_tracer_loc_index(n_loc,j)=1
          end if
       end do
       deallocate(nodes_release)
    end do

    deallocate(temp_arr2d)

    !--------------------------------------------------------------
    ! in case release in volume

    if(age_release_in_volume) then
       do i=1,ToDim_nod2d
          row=nod3d_below_nod2d(1,i)
          do j=1, num_age_tracer
             if(age_tracer_loc_index(row,j)==1) then
                do k=2,num_layers_below_nod2d(i)+1
                   n3=nod3d_below_nod2d(k,i)
                   age_tracer_loc_index(n3,j)=1
                end do
             end if
          end do
       end do
    end if

  end subroutine age_tracer_init
  !
  !-------------------------------------------------------------------------
  !
  subroutine age_tracer_tendency

    integer      :: j, elem, elnodes(4), res_ind_elem(4)
    real(kind=8) :: inv20, vol, source_elem(4), sum_source
    real(kind=8) :: inside_val(4), outside_val

    inv20=1./20.
    outside_val=1.0/(86400.0*(365+fleapyear))

    do elem=1,myDim_elem3d              
       elnodes=elem3D_nodes(:,elem)
       vol=voltetra(elem)*inv20   
       do j=1, num_age_tracer
          inside_val=-tracer(elnodes,index_age_tracer(j))/age_tracer_restore_time
          res_ind_elem=age_tracer_loc_index(elnodes,j) ! 1-inside, 0-outside
          source_elem=inside_val*res_ind_elem + outside_val*(1.0-res_ind_elem)
          sum_source=sum(source_elem)
          tracer_rhs(elnodes,index_age_tracer(j))=tracer_rhs(elnodes,index_age_tracer(j)) &
               + (sum_source+source_elem(elnodes))*vol
       end do
    end do

  end subroutine age_tracer_tendency
  !
  !-------------------------------------------------------------------------
  !
  subroutine age_tracer_cutoff_restore

    integer   :: j, row

    do j=1, num_age_tracer
       do row=1,ToDim_nod3d
          tracer(row,index_age_tracer(j))= &
               max(tracer(row,index_age_tracer(j)),0.)
          if(zero_age_at_release .and. age_tracer_loc_index(row,j)==1) then
             tracer(row,index_age_tracer(j))=0.0
          end if
       end do
    end do

  end subroutine age_tracer_cutoff_restore

end module o_age_tracer_mod
