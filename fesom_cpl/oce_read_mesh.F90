!--------------------------------------------------------------
! Reads mesh and communication information in a distributed way
! some extra info. (nodal flag for cavity, region type
! of 2d elements, sigma grid slope) is also read in here.

subroutine read_mesh
  use o_PARAM
  use o_MESH
  use o_ELEMENTS
  use o_MATRICES
  use o_ARRAY
  use g_config
  use g_PARfe 
  use g_rotate_grid
  use ocean_mesh_module, only : ocean_mesh__add_depth, ocean_mesh__add_wetrange, ocean_mesh__trim_empty_wetranges, ocean_mesh__depths, ocean_mesh__global_wetranges
  use mpi_topology_module, only : mpi_topology
  use fesom_env_module, only : fesom_env_mesh_partition_path, fesom_env__readable
  implicit none

  integer       :: n, m, k, nn, fileID, ind
  integer       :: vert_nodes(100)
  integer       :: nchunk, chunk_size, ipos, iofs, mesh_check
  real(kind=8)  :: x, y, z, rx, ry
  real(kind=8)  :: t0, t1
  character*10  :: mype_string
  character*300 :: file_name
  character(:), allocatable :: dist_mesh_dir
  integer       :: ierror              ! return error code
  integer, allocatable, dimension(:)        :: mapping
  integer, allocatable, dimension(:,:)      :: ibuff
  real(kind=8), allocatable, dimension(:,:) :: rbuff
  integer, allocatable, dimension(:,:)      :: auxbuff ! will be used for reading aux3d.out
  integer last_wet_2Dindex
  integer levelindex, wetcolindex, wetcolsize, totalcolpos
  real(kind=8) lastdepth
  integer l, r
  logical am_i_output_writer ! save some RAM and do not store wetranges on nodes which do not write data
  
  !mesh related files will be read in chunks of chunk_size
  chunk_size=100000
  !==============================
  ! Allocate mapping array (chunk_size)
  ! It will be used for several purposes 
  !==============================
  allocate(mapping(chunk_size))
  allocate(ibuff(chunk_size,4), rbuff(chunk_size,3))

  mapping=0 
  !==============================
  t0=MPI_Wtime()
  write(mype_string,'(i4.4)') mype  
  dist_mesh_dir = fesom_env_mesh_partition_path(trim(meshpath), npes)

  !=======================
  ! rank partitioning vectors
  !=======================
  file_name=trim(dist_mesh_dir)//'/rpart.out' 
  fileID=200+mype
  open(fileID, file=trim(file_name)) 
  allocate(part2D(npes+1), part3D(npes+1))

  read(fileID,*) n
  if (n.ne.npes) then
     write(*,*) 'current NPES does not coincide with that used for mesh pre-partition'
     call par_ex
     stop
  end if
  part2D(1)=1
  read(fileID,*) part2D(2:npes+1)
  do n=2, npes+1
     part2D(n)=part2D(n-1)+part2D(n)
  end do

  part3D(1)=1
  read(fileID,*) part3D(2:npes+1)
  do n=2, npes+1
     part3D(n)=part3D(n-1)+part3D(n)
  end do
  close(fileID)
  !write(*,*) 'rpart is read'

  !===========================
  ! Lists of nodes and elements 
  ! in global indexing. Not everything
  ! is needed
  !===========================

  file_name=trim(dist_mesh_dir)//'/my_list'//trim(mype_string)//'.out'  
  fileID=200+mype  

  open(fileID, file=trim(file_name))
  read(fileID,*) n

  read(fileID,*) myDim_nod2D
  read(fileID,*) eDim_nod2D
  ToDim_nod2d=myDim_nod2d+eDim_nod2d
  allocate(myList_nod2D(myDim_nod2D+eDim_nod2D)) 	 
  read(fileID,*) myList_nod2D

  read(fileID,*) myDim_nod3D
  read(fileID,*) eDim_nod3D 	 
  ToDim_nod3d=myDim_nod3d+eDim_nod3d
  allocate(myList_nod3D(myDim_nod3D+eDim_nod3D)) 	 
  read(fileID,*) myList_nod3D

  read(fileID,*) myDim_elem2D
  allocate(myList_elem2D(myDim_elem2D))
  read(fileID,*) myList_elem2D

  read(fileID,*) myDim_elem3D
  allocate(myList_elem3D(myDim_elem3D))
  read(fileID,*) myList_elem3D ! m

  close(fileID)

  nod3D=part3D(npes+1)-1
  nod2D=part2D(npes+1)-1

  !==============================
  ! read 2d node data
  !==============================

  allocate(coord_nod2D(2,myDim_nod2D+eDim_nod2D))
  allocate(index_nod2D(myDim_nod2D+eDim_nod2D))	 			

  if (mype==0) then
    file_name=trim(meshpath)//'nod2d.out'
    open(fileID, file=file_name)
    read(fileID,*) n      ! nod2D, we know it already
    write(*,*) 'reading '// trim(file_name)   
  end if
  mesh_check=0
  do nchunk=0, (nod2D-1)/chunk_size
     mapping(1:chunk_size)=0
     do n=1, myDim_nod2D+eDim_nod2D
        ipos=(myList_nod2D(n)-1)/chunk_size
        if (ipos==nchunk) then
           iofs=myList_nod2D(n)-nchunk*chunk_size
           mapping(iofs)=n
        end if
     end do

     k=min(chunk_size, nod2D-nchunk*chunk_size)
     if (mype==0) then
        do n=1, k
           read(fileID,*) ibuff(n,1), rbuff(n,1:2), ibuff(n,2)
        end do
     end if
     call MPI_BCast(rbuff(1:k,1), k, MPI_DOUBLE_PRECISION, 0, MPI_COMM_FESOM, ierror)
     call MPI_BCast(rbuff(1:k,2), k, MPI_DOUBLE_PRECISION, 0, MPI_COMM_FESOM, ierror)
     call MPI_BCast(ibuff(1:k,2), k, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)

     do n=1, k
        x=rbuff(n,1)
        y=rbuff(n,2)
        if (force_rotation) then
           rx=x
           ry=y
           call g2r(rx*rad, ry*rad, x, y)
           x=x/rad
           y=y/rad
        end if        
        if (mapping(n)>0) then
           mesh_check=mesh_check+1
           coord_nod2D(1,mapping(n))=x
           coord_nod2D(2,mapping(n))=y
           index_nod2D(mapping(n))=ibuff(n,2)
        end if
     end do
  end do
  if (mype==0) close(fileID)

  if (mesh_check/=myDim_nod2D+eDim_nod2D) then
     write(*,*) 'ERROR while reading nod2d.out on mype=', mype
     write(*,*) mesh_check, ' values have been read in according to partitioning'
     write(*,*) 'it does not equal to myDim_nod2D+eDim_nod2D = ', myDim_nod2D+eDim_nod2D
  end if

  !==============================
  ! read 2d elem data
  !==============================
  if (mype==0)  then 
     file_name=trim(meshpath)//'elem2d.out'
     open(fileID, file=file_name)
     read(fileID,*) elem2d
     write(*,*) 'reading '// trim(file_name)   
  end if
  call MPI_BCast(elem2d, 1, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)

  allocate(elem2D_nodes(3, myDim_elem2D))
  do nchunk=0, (elem2D-1)/chunk_size
     mapping(1:chunk_size)=0
     do n=1, myDim_elem2D
        ipos=(myList_elem2D(n)-1)/chunk_size
        if (ipos==nchunk) then
           iofs=myList_elem2D(n)-nchunk*chunk_size
           mapping(iofs)=n
        end if
     end do

     k=min(chunk_size, elem2D-nchunk*chunk_size)
     if (mype==0) then
        do n=1, k
           read(fileID,*) ibuff(n, 1:3)
        end do
     end if

     call MPI_BCast(ibuff(1:k,1), k, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)
     call MPI_BCast(ibuff(1:k,2), k, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)
     call MPI_BCast(ibuff(1:k,3), k, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)

     do n=1, k
        if (mapping(n)>0) then
           elem2D_nodes(1,mapping(n))=ibuff(n,1)
           elem2D_nodes(2,mapping(n))=ibuff(n,2)
           elem2D_nodes(3,mapping(n))=ibuff(n,3)
        end if
     end do
  end do
  if (mype==0) close(fileID)
  ! nodes in elem2d are in global numbering. convert to local:
  do nchunk=0, (nod2D-1)/chunk_size
     mapping(1:chunk_size)=0
     do n=1, myDim_nod2D+eDim_nod2D
        ipos=(myList_nod2D(n)-1)/chunk_size
        if (ipos==nchunk) then
           iofs=myList_nod2D(n)-nchunk*chunk_size
           mapping(iofs)=n
        end if
     end do
     do n=1, myDim_elem2D
        do m=1,3
           nn=elem2D_nodes(m, n)
           ipos=(nn-1)/chunk_size
           if (ipos==nchunk) then
              iofs=nn-nchunk*chunk_size
              ! minus sign is required to avoid modified entry being modified in another chunk
              ! will be changed to plus at the end
              elem2D_nodes(m,n)=-mapping(iofs) 
           end if
        end do
     end do
  end do
  elem2D_nodes=-elem2D_nodes

  !==============================
  ! Ice shelf variables
  !==============================
  allocate(cavity_flag_nod2d(myDim_nod2d+eDim_nod2d))
  cavity_flag_nod2d=0
#ifdef use_cavity
  if (mype==0) then 
     file_name=trim(meshpath)//'cavity_flag_nod2d.out'
     open(fileID, file=file_name)
     write(*,*) 'reading '// trim(file_name)   
  end if
  do nchunk=0, (nod2D-1)/chunk_size
     mapping(1:chunk_size)=0
     do n=1, myDim_nod2D+eDim_nod2D
        ipos=(myList_nod2D(n)-1)/chunk_size
        if (ipos==nchunk) then
           iofs=myList_nod2D(n)-nchunk*chunk_size
           mapping(iofs)=n
        end if
     end do

     k=min(chunk_size, nod2D-nchunk*chunk_size)
     if (mype==0) then
        do n=1, k
           read(fileID,*) ibuff(n, 1)
        end do
     end if
     call MPI_BCast(ibuff(1:k,1), k, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)
     do n=1, k
        if(mapping(n)>0) cavity_flag_nod2d(mapping(n))=ibuff(n, 1)
     end do
  end do
  if (mype==0) close(fileID)
#endif


  !=============================
  ! Region type of 2d elements: z-level or sigma
  !=============================
  allocate(grid_type_elem2d(myDim_elem2d))
  if(grid_type==1) then
     grid_type_elem2d=0
  elseif(grid_type==2) then
     grid_type_elem2d=1
  else
     if (mype==0) then 
        file_name=trim(meshpath)//'grid_type_elem2d.out'
        open(fileID, file=file_name)
        write(*,*) 'reading '// trim(file_name)   
     end if
     do nchunk=0, (elem2D-1)/chunk_size
        mapping(1:chunk_size)=0
        do n=1, myDim_elem2D
           ipos=(myList_elem2D(n)-1)/chunk_size
           if (ipos==nchunk) then
              iofs=myList_elem2D(n)-nchunk*chunk_size
              mapping(iofs)=n
           end if
        end do

        k=min(chunk_size, elem2D-nchunk*chunk_size)
        if (mype==0) then
           do n=1, k
              read(fileID,*) ibuff(n, 1)
           end do
        end if

        call MPI_BCast(ibuff(1:k,1), k, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)

        do n=1, k
           read(fileID,*) nn
           if(mapping(n)>0) grid_type_elem2d(mapping(n))=ibuff(n, 1)
        end do
     end do
     if (mype==0) close(fileID)
  end if

  !==============================
  ! read 3d node data
  !==============================
  am_i_output_writer = mpi_topology%am_i_host_head_rank(MPI_COMM_FESOM)
  lastdepth = 1
  allocate(coord_nod3D(3,myDim_nod3D+eDim_nod3D))
  allocate(index_nod3D(myDim_nod3D+eDim_nod3D))	 			
  if (mype==0) then 
     file_name=trim(meshpath)//'nod3d.out' 
     open(fileID, file=file_name)
     read(fileID,*) n      ! nod3D, we know it already
     write(*,*) 'reading '// trim(file_name)   
  end if
  do nchunk=0, (nod3D-1)/chunk_size
     mapping(1:chunk_size)=0
     do n=1, myDim_nod3D+eDim_nod3D
        ipos=(myList_nod3D(n)-1)/chunk_size
        if (ipos==nchunk) then
           iofs=myList_nod3D(n)-nchunk*chunk_size
           mapping(iofs)=n
        end if
     end do

     k=min(chunk_size, nod3D-nchunk*chunk_size)
     if (mype==0) then
        do n=1, k
           read(fileID,*) ibuff(n,1), rbuff(n,1:3), ibuff(n,2)
        end do
     end if

     call MPI_BCast(rbuff(1:k,1), k, MPI_DOUBLE_PRECISION, 0, MPI_COMM_FESOM, ierror)
     call MPI_BCast(rbuff(1:k,2), k, MPI_DOUBLE_PRECISION, 0, MPI_COMM_FESOM, ierror)
     call MPI_BCast(rbuff(1:k,3), k, MPI_DOUBLE_PRECISION, 0, MPI_COMM_FESOM, ierror)
     call MPI_BCast(ibuff(1:k,2), k, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)

     do n=1, k
        ! grab the depth for each level
        if( abs(lastdepth-rbuff(n,3)) > epsilon(lastdepth) ) then ! i.e. lastdepth > newdepth
          if(am_i_output_writer) call ocean_mesh__add_depth(-rbuff(n,3)) ! flip sign to have depth do positive downwards
          lastdepth = rbuff(n,3)
        end if
     
        x=rbuff(n,1)
        y=rbuff(n,2)
        if (mapping(n)>0) then
        if (force_rotation) then
           rx=x
           ry=y
           call g2r(rx*rad, ry*rad, x, y)
           x=x/rad
           y=y/rad
        end if
        coord_nod3D(1,mapping(n))=x
        coord_nod3D(2,mapping(n))=y
        coord_nod3D(3,mapping(n))=rbuff(n,3)
        index_nod3D(mapping(n))=ibuff(n,2)
        end if 
     end do
  end do
  if (mype==0) close(fileID)

  !==============================
  ! read 3d elem data
  !==============================

  if (mype==0) then
     file_name=trim(meshpath)//'elem3d.out' 
     open(fileID, file=file_name)
     read(fileID,*) elem3D
     write(*,*) 'reading '// trim(file_name)   
  end if
  call MPI_BCast(elem3d, 1, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)
  allocate(elem3D_nodes(4, myDim_elem3D))
  do nchunk=0, (elem3D-1)/chunk_size
     mapping(1:chunk_size)=0
     do n=1, myDim_elem3D
        ipos=(myList_elem3D(n)-1)/chunk_size
        if (ipos==nchunk) then
           iofs=myList_elem3D(n)-nchunk*chunk_size
           mapping(iofs)=n
        end if
     end do

     k=min(chunk_size, elem3D-nchunk*chunk_size)
     if (mype==0) then
        do n=1, k
           read(fileID,*) ibuff(n,1:4)
        end do
     end if
     call MPI_BCast(ibuff(1:k,1), k, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)
     call MPI_BCast(ibuff(1:k,2), k, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)
     call MPI_BCast(ibuff(1:k,3), k, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)
     call MPI_BCast(ibuff(1:k,4), k, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)

     do n=1, min(chunk_size, elem3D-nchunk*chunk_size)
        if (mapping(n)>0) then
           elem3D_nodes(1,mapping(n))=ibuff(n,1)
           elem3D_nodes(2,mapping(n))=ibuff(n,2)
           elem3D_nodes(3,mapping(n))=ibuff(n,3)
           elem3D_nodes(4,mapping(n))=ibuff(n,4)
        end if
      end do
  end do
  if (mype==0) close(fileID)

  ! nodes in elem3d are in natural numbering. convert to local:
  do nchunk=0, (nod3D-1)/chunk_size
     mapping(1:chunk_size)=0
     do n=1, myDim_nod3D+eDim_nod3D
        ipos=(myList_nod3D(n)-1)/chunk_size
        if (ipos==nchunk) then
           iofs=myList_nod3D(n)-nchunk*chunk_size
           mapping(iofs)=n
        end if
     end do
     do n=1, myDim_elem3D
        do m=1, 4
           nn=elem3D_nodes(m, n)
           ipos=(nn-1)/chunk_size
           if (ipos==nchunk) then
              iofs=nn-nchunk*chunk_size
              ! minus sign is required to avoid modified entry being modified in another chunk
              ! will be changed to plus at the end
              elem3D_nodes(m,n)=-mapping(iofs)
           end if
        end do
     end do
  end do
  elem3D_nodes=-elem3D_nodes

  !==============================
  ! read aux. arrays 
  !==============================

  if (mype==0) then 
     file_name=trim(meshpath)//'aux3d.out'  
     open(fileID, file=file_name)
     read(fileID, *) max_num_layers
     write(*,*) 'reading '// trim(file_name)   
  end if
  call MPI_BCast(max_num_layers, 1, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)
  allocate(auxbuff(chunk_size, max_num_layers)) ! will be deallocated after having read aux3d.out
  !=============================
  ! nod3D_below_nod2D
  !============================= 
  ! ATTENTION: the array is to be stored in slices of max_num_layers
  !            or as a column
  allocate(nod3D_below_nod2D(max_num_layers,myDim_nod2D+eDim_nod2D))
  do nchunk=0, (nod2D-1)/chunk_size
     mapping(1:chunk_size)=0
     do n=1, myDim_nod2D+eDim_nod2D
        ipos=(myList_nod2D(n)-1)/chunk_size
        if (ipos==nchunk) then
           iofs=myList_nod2D(n)-nchunk*chunk_size
           mapping(iofs)=n
        end if
     end do

     k=min(chunk_size, nod2D-nchunk*chunk_size)
     if (mype==0) then
        do n=1, k
           read(fileID,*) auxbuff(n, 1:max_num_layers) ! i.e. each 2D-node
        end do
     end if
     ! broadcast for each layer
     do n=1, max_num_layers
        call MPI_BCast(auxbuff(1:k,n), k, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)
     end do
    ! here we have the -999s on each process available in auxbuff
    ! auxbuf(each 2D-node,48)
    ! nod3D_below_nod2D(48:local2Dnodes+ghostnodes)
    do n=1, k ! each 2D-node        
      if (mapping(n)>0)  then ! select nodes belonging to this process
           nod3D_below_nod2D(:,mapping(n)) = auxbuff(n,:) ! for node n, all layers from auxbuf to nod3D_below_nod2D
      end if
    end do
  end do
  deallocate(auxbuff) ! reading finished

  where (nod3D_below_nod2D <= 0) ! -999 gets replaced with 0
        nod3D_below_nod2D=0
  end where
  do nchunk=0, (nod3D-1)/chunk_size
     mapping(1:chunk_size)=0
     do n=1, myDim_nod3D+eDim_nod3D
        ipos=(myList_nod3D(n)-1)/chunk_size
        if (ipos==nchunk) then
           iofs=myList_nod3D(n)-nchunk*chunk_size
           mapping(iofs)=n
        end if
     end do
     do n=1, myDim_nod2D+eDim_nod2D
        do m=1, max_num_layers
           nn=nod3D_below_nod2D(m,n)
           ipos=(nn-1)/chunk_size
           if (nn > 0 .and. ipos==nchunk) then
              iofs=nn-nchunk*chunk_size
              ! minus sign is required to avoid modified entry being modified in another chunk
              ! will be changed to plus at the end
              nod3D_below_nod2D(m,n)=-mapping(iofs)
           end if           
        end do
      end do
  end do
  nod3D_below_nod2D=-nod3D_below_nod2D
  where (nod3D_below_nod2D == 0)
        nod3D_below_nod2D=-9999
  end where  
  !=============================
  ! nod2D_corresp_to_nod3D
  !============================= 

  last_wet_2Dindex = nod2D
  levelindex = 0
  wetcolsize = 0
  allocate(nod2D_corresp_to_nod3D(myDim_nod3D+eDim_nod3D)) 
  do nchunk=0, (nod3D-1)/chunk_size
     mapping(1:chunk_size)=0
     do n=1, myDim_nod3D+eDim_nod3D
        ipos=(myList_nod3D(n)-1)/chunk_size
        if (ipos==nchunk) then
           iofs=myList_nod3D(n)-nchunk*chunk_size
           mapping(iofs)=n
        end if
     end do

     k=min(chunk_size, nod3D-nchunk*chunk_size)
     if (mype==0) then
        do n=1, k
           read(fileID,*) ibuff(n,1) ! each wet 3D node
        end do
     end if
     call MPI_BCast(ibuff(1:k,1), k, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)

    do n=1, k ! each 3D wet node
      
      if(last_wet_2Dindex >= ibuff(n,1)) then
        ! we entered a new level, ibuff(n,1) is the first wet node
        if(wetcolsize > 0) then ! store last wetrange
          if(am_i_output_writer) call ocean_mesh__add_wetrange(levelindex, wetcolindex, wetcolsize)
        end if
        levelindex = levelindex +1
        totalcolpos = ibuff(n,1)
        wetcolindex = ibuff(n,1)
        wetcolsize = 1
      else
        ! same level as before
        if(last_wet_2Dindex+1 .eq. ibuff(n,1)) then
          ! the previous node was also wet, no new wet range
          totalcolpos = totalcolpos +1
          wetcolsize = wetcolsize +1
        else
          ! the previous node was a dry node, so a new wet range starts here
          if(wetcolsize > 0) then ! store last wetrange
            if(am_i_output_writer) call ocean_mesh__add_wetrange(levelindex, wetcolindex, wetcolsize)
          end if
          totalcolpos = totalcolpos + ibuff(n,1) - last_wet_2Dindex ! skip dry nodes
          wetcolindex = totalcolpos
          wetcolsize = 1
        end if
      end if
      
      last_wet_2Dindex = ibuff(n,1)
            
      if (mapping(n)>0)  then ! select nodes belonging to this process
        nod2D_corresp_to_nod3D(mapping(n))=ibuff(n,1)
      end if
     end do
  end do
  if(wetcolsize > 0) then ! store last wetrange
    if(am_i_output_writer) call ocean_mesh__add_wetrange(levelindex, wetcolindex, wetcolsize)
  end if

  ! remove any empty levels from the end of wetranges (i.e. levels without any ocean nodes)
  call ocean_mesh__trim_empty_wetranges

  ! put to wetranges and depth c++ data storage
  if(am_i_output_writer) then
    do l=1,size(ocean_mesh__global_wetranges)
      do r=1,size(ocean_mesh__global_wetranges(l)%ranges)
        call add_ocean_levels_wetrange(l, ocean_mesh__global_wetranges(l)%ranges(r)%index, ocean_mesh__global_wetranges(l)%ranges(r)%size)
      end do
    end do
    do l=1,size(ocean_mesh__depths)
      call add_ocean_depth(ocean_mesh__depths(l))
    end do
  end if

  do nchunk=0, (nod2D-1)/chunk_size
     mapping(1:chunk_size)=0
     do n=1, myDim_nod2D+eDim_nod2D
        ipos=(myList_nod2D(n)-1)/chunk_size
        if (ipos==nchunk) then
           iofs=myList_nod2D(n)-nchunk*chunk_size
           mapping(iofs)=n
        end if
     end do
     do n=1, myDim_nod3D+eDim_nod3D
        m=nod2D_corresp_to_nod3D(n)
        ipos=(m-1)/chunk_size
        if (ipos==nchunk) then
           iofs=m-nchunk*chunk_size
           ! minus sign is required to avoid modified entry being modified in another chunk
           ! will be changed to plus at the end
           nod2D_corresp_to_nod3D(n)=-mapping(iofs)
        end if
     end do
  end do
  nod2D_corresp_to_nod3D=-nod2D_corresp_to_nod3D
  !=============================
  ! elem2D_corresp_to_elem3D
  !=============================
  allocate(elem2D_corresp_to_elem3D(myDim_elem3D)) 

  do nchunk=0, (elem3D-1)/chunk_size
     mapping(1:chunk_size)=0
     do n=1, myDim_elem3D
        ipos=(myList_elem3D(n)-1)/chunk_size
        if (ipos==nchunk) then
           iofs=myList_elem3D(n)-nchunk*chunk_size
           mapping(iofs)=n
        end if
     end do

     k=min(chunk_size, elem3D-nchunk*chunk_size)
     if (mype==0) then
        do n=1, k
           read(fileID,*) ibuff(n,1)
        end do
     end if
     call MPI_BCast(ibuff(1:k,1), k, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)

     do n=1, k
        if (mapping(n)>0) then
           elem2D_corresp_to_elem3D(mapping(n))=ibuff(n,1)
        end if
     end do
  end do

  do nchunk=0, (elem2D-1)/chunk_size
     mapping(1:chunk_size)=0
     do n=1, myDim_elem2D
        ipos=(myList_elem2D(n)-1)/chunk_size
        if (ipos==nchunk) then
           iofs=myList_elem2D(n)-nchunk*chunk_size
           mapping(iofs)=n
        end if
     end do
     do n=1, myDim_elem3D
        m=elem2D_corresp_to_elem3D(n)
        ipos=(m-1)/chunk_size
        if (ipos==nchunk) then
           iofs=m-nchunk*chunk_size
           ! minus sign is required to avoid modified entry being modified in another chunk
           ! will be changed to plus at the end
           elem2D_corresp_to_elem3D(n)=-mapping(iofs)
        end if
     end do
  end do
  elem2D_corresp_to_elem3D=-elem2D_corresp_to_elem3D
  if (mype==0) close(fileID)

  !correcting 3d nodal indices
#ifndef use_opbnd_restoring
#ifndef use_opbnd_tide
  do n=1,myDim_nod3d+eDim_nod3D
     if (index_nod3D(n)==12) index_nod3D(n)=11
     if (index_nod3D(n)==22) index_nod3D(n)=21
     if (index_nod3D(n)==32) index_nod3D(n)=31
  end do
#endif
#endif

  ! ==============================
  ! read sigma slope and define layer of elem.
  ! ==============================
  if(grid_type/=1) call read_grid_slope	
   

  if(mype==0) then
     write(*,*) 'mesh (according to pre-partition) is read in'	
     write(*,*) 'configured with nod2D=',nod2D,' nod3D=',nod3D
  end if

  ! ==============================
  ! Communication information
  ! ==============================
  file_name=trim(dist_mesh_dir)//'/com_info'//trim(mype_string)//'.out'  
  fileID=200+mype  
  open(fileID, file=file_name)
  read(fileID,*)  n
  read(fileID,*) com_nod2D%rPEnum
  allocate(com_nod2D%rPE(com_nod2D%rPEnum))
  read(fileID,*) com_nod2D%rPE
  allocate(com_nod2D%rptr(com_nod2D%rPEnum+1))
  read(fileID,*) com_nod2D%rptr
  allocate(com_nod2D%rlist(eDim_nod2D))
  read(fileID,*) com_nod2D%rlist

  read(fileID,*) com_nod2D%sPEnum
  allocate(com_nod2D%sPE(com_nod2D%sPEnum))
  read(fileID,*) com_nod2D%sPE
  allocate(com_nod2D%sptr(com_nod2D%sPEnum+1))
  read(fileID,*) com_nod2D%sptr
  n=com_nod2D%sptr(com_nod2D%sPEnum+1)-1
  allocate(com_nod2D%slist(n))
  read(fileID,*) com_nod2D%slist

  read(fileID,*) com_nod3D%rPEnum
  allocate(com_nod3D%rPE(com_nod3D%rPEnum))
  read(fileID,*) com_nod3D%rPE
  allocate(com_nod3D%rptr(com_nod3D%rPEnum+1))
  read(fileID,*) com_nod3D%rptr
  allocate(com_nod3D%rlist(eDim_nod3D))
  read(fileID,*) com_nod3D%rlist

  read(fileID,*) com_nod3D%sPEnum
  allocate(com_nod3D%sPE(com_nod3D%sPEnum))
  read(fileID,*) com_nod3D%sPE
  allocate(com_nod3D%sptr(com_nod3D%sPEnum+1))
  read(fileID,*) com_nod3D%sptr
  n=com_nod3D%sptr(com_nod3D%sPEnum+1)-1
  allocate(com_nod3D%slist(n))
  read(fileID,*) com_nod3D%slist

  !contiguous listing will be read
  allocate(myCList_nod2D(myDim_nod2D+eDim_nod2D))
  allocate(myCList_nod3D(myDim_nod3D+eDim_nod3D))

  ! two versions of partitioning are available currenlty and differ
  ! by storing mapping data in a separate file instead of in every com_info file
  file_name=trim(dist_mesh_dir)//'/mapping'//'.out'
  if(fesom_env__readable(file_name)) then ! mapping.out file is present, so we can speed up reading the partitioned mesh and do not read everything on every process
  close(fileID) ! close com_info files and open mapping.out

  if (mype==0) then
     fileID=201
     open(fileID, file=file_name, action='read')
  end if

  !3D part
  do nchunk=0, (nod3D-1)/chunk_size
     mapping(1:chunk_size)=0
     do n=1, myDim_nod3D+eDim_nod3D
        ipos=(myList_nod3D(n)-1)/chunk_size
        if (ipos==nchunk) then
           iofs=myList_nod3D(n)-nchunk*chunk_size
           mapping(iofs)=n
        end if
     end do


     k=min(chunk_size, nod3D-nchunk*chunk_size)
     if (mype==0) then
        do n=1, k
           read(fileID,*) ibuff(n,1)
        end do
     end if
     call MPI_BCast(ibuff(1:k,1), k, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)

     do n=1, k
        if (mapping(n)>0) then
           myCList_nod3D(mapping(n))=ibuff(n,1)
        end if 
     end do
  end do

  !2D part
  do nchunk=0, (nod2D-1)/chunk_size
     mapping(1:chunk_size)=0
     do n=1, myDim_nod2D+eDim_nod2D
        ipos=(myList_nod2D(n)-1)/chunk_size
        if (ipos==nchunk) then
           iofs=myList_nod2D(n)-nchunk*chunk_size
           mapping(iofs)=n
        end if
     end do

     k=min(chunk_size, nod2D-nchunk*chunk_size)
     if (mype==0) then
        do n=1, k
           read(fileID,*) ibuff(n,1)
        end do
     end if
     call MPI_BCast(ibuff(1:k,1), k, MPI_INTEGER, 0, MPI_COMM_FESOM, ierror)

     do n=1, k
        if (mapping(n)>0) then
           myCList_nod2D(mapping(n))=ibuff(n,1)
        end if 
     end do
  end do

  else ! continue reading mapping array from com_info file on each core
  !3D part
  do nchunk=0, (nod3D-1)/chunk_size
     mapping(1:chunk_size)=0
     do n=1, myDim_nod3D+eDim_nod3D
        ipos=(myList_nod3D(n)-1)/chunk_size
        if (ipos==nchunk) then
           iofs=myList_nod3D(n)-nchunk*chunk_size
           mapping(iofs)=n
        end if
     end do
     do n=1, min(chunk_size, nod3D-nchunk*chunk_size)
        read(fileID,*) nn
        if (mapping(n)>0) then
           myCList_nod3D(mapping(n))=nn
        end if 
     end do
  end do

  !2D part
  do nchunk=0, (nod2D-1)/chunk_size
     mapping(1:chunk_size)=0
     do n=1, myDim_nod2D+eDim_nod2D
        ipos=(myList_nod2D(n)-1)/chunk_size
        if (ipos==nchunk) then
           iofs=myList_nod2D(n)-nchunk*chunk_size
           mapping(iofs)=n
        end if
     end do
     do n=1, min(chunk_size, nod2D-nchunk*chunk_size)
        read(fileID,*) nn
        if (mapping(n)>0) then
           myCList_nod2D(mapping(n))=nn
        end if 
     end do
  end do
  end if
  close(fileID)
  deallocate(rbuff, ibuff)
  deallocate(mapping)
  allocate(col_pos(myDim_nod3D+eDim_nod3D))

  col_pos=0
  t1=MPI_Wtime()

  if(mype==0) write(*,*) 'comm and mapping is read in ', t1-t0, ' seconds'
end subroutine  read_mesh
!
!=========================================================================
!
subroutine read_grid_slope  
  ! read sigma grid slope
  use o_param
  use o_MESH
  use o_ELEMENTS
  use g_config
  use g_parfe

  implicit none

  integer                               :: i, j, k, n, fileID
  real(kind=8)                          :: temp(2)
  integer, allocatable, dimension(:)    :: mapping
  integer                               :: nchunk, chunk_size, ipos, iofs
  integer                               :: ierror ! return error code
  real(kind=8), allocatable, dimension(:,:,:) :: rbuff
  ! sigma grid slope
  allocate(grid_slope(2, max_num_layers-1,myDim_elem2d))
  grid_slope=0.0

  !mesh related files will be read in chunks of chunk_size
  chunk_size=100000
  allocate(rbuff(chunk_size, max_num_layers-1, 2))
  !==============================
  ! Allocate mapping array (chunk_size)
  ! It will be used for several purposes 
  !==============================
  allocate(mapping(chunk_size))
  mapping=0 

  if (mype==0) fileID=101
  if (mype==0) open(fileID,file=trim(MeshPath)//'sigma_grid_slope_elem.out', status='old')

  do nchunk=0, (elem2D-1)/chunk_size
     mapping(1:chunk_size)=0
     do n=1, myDim_elem2D
        ipos=(myList_elem2D(n)-1)/chunk_size
        if (ipos==nchunk) then
           iofs=myList_elem2D(n)-nchunk*chunk_size
           mapping(iofs)=n
        end if
     end do

     k=min(chunk_size, elem2D-nchunk*chunk_size)
     if (mype==0) then
        do n=1, k
           do j=1,max_num_layers-1
              read(fileID,*) rbuff(n,j,:)
           end do
        end do
     end if

     do j=1,max_num_layers-1
        call MPI_BCast(rbuff(1:k, j, 1), k, MPI_DOUBLE_PRECISION, 0, MPI_COMM_FESOM, ierror)
        call MPI_BCast(rbuff(1:k, j, 2), k, MPI_DOUBLE_PRECISION, 0, MPI_COMM_FESOM, ierror)
     end do

     do n=1, k
        do j=1,max_num_layers-1
           if (mapping(n)>0) then
              grid_slope(:,j,mapping(n))=rbuff(n, j, :)
           end if
           end do
      end do
  end do
  if (mype==0) close(fileID)
  deallocate(rbuff)
  deallocate(mapping)
end subroutine read_grid_slope
!
!=========================================================================
