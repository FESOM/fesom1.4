! 
! P1-P1 tetrahedral discretization
! 
! Distributed memory setup
!
! Initialization routines: they are performed in global memory model, 
! but only a small part of arrays is allocated. 2Gb should be sufficient for 
! meshes with about 10 M nodes. 
!
!=============================================================================

program MAIN
  use o_PARAM
  use o_ELEMENTS
  use o_MESH
  use g_PARFE
  use o_MATRICES
  use g_config
  implicit none

  character(len=100)   :: nmlfile

  call par_init

  nmlfile ='namelist.config'    ! name of config. namelist file
  open (20,file=nmlfile)
  read (20,NML=paths)
  close (20)

  call read_2Dmesh
  call read_3Dmesh           ! nlayer part is added there      

  call build_nghbr_arrays    ! Is only used for computing 2D neighbourhood 
  call ssh_stiff_construct
  call set_par_support_ini
  call save_dist_mesh
  call par_ex
end program MAIN
!
!--------------------------------------------------------------------------
!
subroutine read_2Dmesh
  use g_config
  use o_PARAM
  use o_MESH
  use o_ELEMENTS
  implicit none

  integer           :: i, n, ind

  open (20,file=trim(MeshPath)//'nod2d.out',  status='old')
  open (21,file=trim(MeshPath)//'elem2d.out', status='old')
  write(*,*) '2D mesh is opened'

  read(20,*) nod2D 
  close(20)
  !
  read(21,*)  elem2D      
  allocate(elem2D_nodes(3,elem2D))  
  do n=1, elem2D
     read(21,*) elem2D_nodes(:,n)
  end do
  close(21)
end subroutine read_2Dmesh
!
!---------------------------------------------------------------------------
!
subroutine read_3Dmesh
  use g_config
  use o_PARAM
  use o_MESH
  use o_ELEMENTS
  implicit none
  
  integer           :: i, j, n, node
   
  open(10,file=trim(meshpath)//'nod3d.out', status='old')
  open(11,file=trim(meshpath)//'elem3d.out',status='old')
  open(12,file=trim(meshpath)//'aux3d.out', status='old')

  write(*,*) '3D mesh is opened'

  ! Read node data
  read(10,*) nod3D   
  close(10)

  ! Read the element data
  read(11, *)  elem3D   
  allocate(elem3D_nodes(4,elem3D)) 
  do i=1,elem3D 
     read(11,*) elem3D_nodes(:,i)
  end do
  close(11)

  ! Read auxilliary data
  read(12,*) max_num_layers
  allocate(nod3D_below_nod2D(max_num_layers,nod2D))
  read(12,*) nod3D_below_nod2D
  close(12)

  allocate(num_layers_below_nod2D(nod2D))
  num_layers_below_nod2D=-1
  do n=1,nod2D
     do j=1,max_num_layers
        node=nod3D_below_nod2D(j,n)
        if (node > 0) then
           num_layers_below_nod2D(n)=num_layers_below_nod2D(n) + 1
        else
           exit
        end if
     end do
  end do
end subroutine read_3Dmesh
!
!---------------------------------------------------------------------------
!
subroutine build_nghbr_arrays
  use o_MESH
  use o_ELEMENTS
  use o_PARAM
  use g_config
  implicit none
  
  integer                            :: j, k, m, a, tr(3), flag, el, ml
  integer                            :: i, pos
  integer, dimension(100)            :: AUX=0
  integer, allocatable, dimension(:) :: ind
  
  !--------------- 2D mesh-------------------------------
  ! Builds nod_in_elem2D  
  allocate(ind(nod2D))
  ind=0
  do j=1,elem2D
     tr=elem2D_nodes(:,j)
     ind(tr)=ind(tr)+1
  end do
  allocate(nod_in_elem2D(nod2D))
  nod_in_elem2D%nmb=ind(1:nod2D)    
  do j=1,nod2D   
     allocate(nod_in_elem2D(j)%addresses(ind(j)))

  end do
  ind=0
  do i=1,elem2D   
     tr=elem2D_nodes(:,i)
     ind(tr)=ind(tr)+1
     do j=1, 3
        nod_in_elem2D(tr(j))%addresses(ind(tr(j)))=i
     end do
  end do
  ! The list of elements is ordered, and no sorting is needed
  
  ! Builds nghbr_nod2D  
  allocate(nghbr_nod2D(nod2D))
  ind=0
  do j=1, nod2D
     flag=0
     do m=1,nod_in_elem2D(j)%nmb
        el=nod_in_elem2D(j)%addresses(m)
        do k=1, 3
           a=elem2D_nodes(k,el)       
           if (ind(a)==0) then  
              ind(a)=1 
              flag=flag+1         
              aux(flag)=a
           end if
        end do
     end do
     nghbr_nod2D(j)%nmb=flag
     allocate(nghbr_nod2D(j)%addresses(flag))
     ! We need to sort array aux(1:flag)
     do m=flag,1,-1
        ml=maxloc(aux(1:flag), 1)
        nghbr_nod2D(j)%addresses(m)=aux(ml)
        ind(aux(ml))=0
        aux(ml)=-999
     end do
  end do
  
  deallocate(ind)
end subroutine build_nghbr_arrays
!
!----------------------------------------------------------------------------
!
subroutine ssh_stiff_construct
  use o_PARAM
  use o_MESH
  use o_ELEMENTS
  use o_MATRICES
  use g_config
  !
  implicit none
  integer                             :: i, j, k, m, q, row, ipos, is, ie
  integer                             :: elem, elem2, elnodes2(3)
  real(kind=8)                        :: dx(3), dy(3), dz
  
  ! a)
  sshstiff%dim=nod2D
  allocate(sshstiff%rowptr(sshstiff%dim+1))
  sshstiff%rowptr(1)=1
  do k=1,nod2D
     sshstiff%rowptr(k+1)=sshstiff%rowptr(k)+nghbr_nod2D(k)%nmb
  end do
  sshstiff%nza=sshstiff%rowptr(nod2D+1)-1
  ! b)
  allocate(sshstiff%colind(sshstiff%nza))
  allocate(sshstiff%values(sshstiff%nza))
  sshstiff%values=0.0
  ! c)
  do k=1,nod2D
     sshstiff%colind(sshstiff%rowptr(k):sshstiff%rowptr(k+1)-1)= &
          nghbr_nod2D(k)%addresses
  end do

end subroutine ssh_stiff_construct
!
!---------------------------------------------------------------------------
!
subroutine com_global2local
  use g_parfe
  use o_MESH
  implicit none

  integer   n
  integer, allocatable, dimension(:) :: temp

  allocate(temp(nod3D)) 

  ! =========
  ! 2D nodes
  ! =========
  ! Replace global numbering with a local one
  temp(1:nod2D)=0
  do n=1, myDim_nod2D+eDim_nod2D
     temp(myList_nod2D(n))=n
  end do
  do n=1, com_nod2D%sptr(com_nod2D%sPEnum+1)-1
     com_nod2D%slist(n)=temp(com_nod2D%slist(n))
  end do

  do n=1, com_nod2D%rptr(com_nod2D%rPEnum+1)-1
     com_nod2D%rlist(n)=temp(com_nod2D%rlist(n))
  end do
  ! com_nod2D%rlist should be myDim_nod2D+1:myDim_nod2D+eDim_nod2D

  ! =========
  ! 3D nodes
  ! =========
  ! Replace global numbering with a local one
  temp(1:nod3D)=0
  do n=1, myDim_nod3D+eDim_nod3D
     temp(myList_nod3D(n))=n
  end do
  do n=1, com_nod3D%sptr(com_nod3D%sPEnum+1)-1
     com_nod3D%slist(n)=temp(com_nod3D%slist(n))
  end do

  do n=1, com_nod3D%rptr(com_nod3D%rPEnum+1)-1
     com_nod3D%rlist(n)=temp(com_nod3D%rlist(n))
  end do
  ! com_nod3D%rlist should be myDim_nod3D+1:myDim_nod3D+eDim_nod3D

  deallocate(temp)
end subroutine com_global2local
!
!=============================================================================
!
subroutine save_dist_mesh
  use o_PARAM
  use o_MESH
  use o_ELEMENTS
  use o_MATRICES
  use o_ARRAY
  use g_config
  use g_parfe 
  implicit none

  integer        n, m, fileID, nend, nini, ed(2)
  character*10   mype_string
  character*300   file_name
  character*300   dist_mesh_dir
  integer, allocatable, dimension(:)  :: temp, ztemp, count

  allocate(temp(elem3D), ztemp(max_num_layers))  ! serves for mapping
  allocate(count(npes+1))
  write(mype_string,'(i4.4)') mype  
  !dist_mesh_dir=trim(meshpath)//'dist/'
  write(dist_mesh_dir,'(a,i0,a1)') trim(meshpath)//'dist', npes, '/'
  ! ==============================
  ! rank partitioning:
  ! ==============================

  if(mype==0) then
     file_name=trim(dist_mesh_dir)//'rpart.out'  
     fileID=200+mype  
     open(fileID, file=file_name)
     count=0;
     do n=1, nod2D
        m=part2D(n);
        count(m+1)=count(m+1)+1
     end do
     write(fileID,*) npes
     write(fileID,*) count(1:npes)
     count=0;
     do n=1, nod3D
        m=part3D(n);
        count(m+1)=count(m+1)+1
     end do
     write(fileID,*) count(1:npes)
     close(fileID)
  end if


  file_name=trim(dist_mesh_dir)//'my_list'//trim(mype_string)//'.out'  
  fileID=200+mype  
  ! =============================   
  ! lists of owned nodes and elements
  ! =============================
  open(fileID, file=file_name)
  write(fileID,*) mype
  write(fileID,*) myDim_nod2D
  write(fileID,*) eDim_nod2D 	 
  write(fileID,*) myList_nod2D

  write(fileID,*) myDim_nod3D
  write(fileID,*) eDim_nod3D 	 
  write(fileID,*) myList_nod3D

  write(fileID,*) myDim_elem2D
  write(fileID,*) myList_elem2D

  write(fileID,*) myDim_elem3D
  write(fileID,*) myList_elem3D
  close(fileID)       

  ! =========================  
  ! communication information
  ! ========================= 
  call com_global2local   
  file_name=trim(dist_mesh_dir)//'com_info'//trim(mype_string)//'.out' 
  fileID=200+mype  
  open(fileID, file=file_name)
  write(fileID,*) mype
  write(fileID,*) com_nod2D%rPEnum
  write(fileID,*) com_nod2D%rPE
  write(fileID,*) com_nod2D%rptr
  write(fileID,*) com_nod2D%rlist
  write(fileID,*) com_nod2D%sPEnum
  write(fileID,*) com_nod2D%sPE
  write(fileID,*) com_nod2D%sptr
  write(fileID,*) com_nod2D%slist

  write(fileID,*) com_nod3D%rPEnum
  write(fileID,*) com_nod3D%rPE
  write(fileID,*) com_nod3D%rptr
  write(fileID,*) com_nod3D%rlist
  write(fileID,*) com_nod3D%sPEnum
  write(fileID,*) com_nod3D%sPE
  write(fileID,*) com_nod3D%sptr
  write(fileID,*) com_nod3D%slist

  ! ================================
  ! mapping ( PE contiguous 3D numbering) 	 
  ! ================================  
  ! matrices use contiguous numbering
  ! so the correspondence array is required 
  ! to assembly them correctly 
  count=0;
  do n=1, nod3D
     m=part3D(n);
     count(m+2)=count(m+2)+1
     temp(n)=count(m+2)
  end do
  count(1)=1
  do n=2,npes+1	  
     count(n)=count(n)+count(n-1)
  end do
  ! Now count == part in range partitioning   

  do n=1,nod3D
     temp(n)=temp(n)+count(part3D(n)+1)-1
     write(fileID,*) temp(n)  
  end do

  ! ================================
  ! mapping ( PE contiguous 2D numbering) 	 
  ! ================================  

  count=0;
  do n=1, nod2D
     m=part2D(n);
     count(m+2)=count(m+2)+1
     temp(n)=count(m+2)
  end do
  count(1)=1
  do n=2,npes+1	  
     count(n)=count(n)+count(n-1)
  end do
  ! Now count == part in range partitioning   

  do n=1,nod2D
     temp(n)=temp(n)+count(part2D(n)+1)-1
     write(fileID,*) temp(n)  
  end do
  close(fileID)

  deallocate(count, ztemp, temp)

  write(*,*) 'Distributed mesh is saved'
end subroutine  save_dist_mesh














