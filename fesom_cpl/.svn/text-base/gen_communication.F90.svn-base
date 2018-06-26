subroutine com_2D(arr2d)
  use o_MESH
  use i_ARRAY
  use g_PARFE 
  implicit none

  real(kind=8), intent(inout)	:: arr2d(ToDim_nod2D)

  integer       :: sreq(com_nod2D%sPEnum)
  integer       :: rreq(com_nod2D%rPEnum)
  integer       :: n, sn, rn, nend
  integer       :: count, source
  real(kind=8)  :: s_buff(com_nod2D%sptr(com_nod2D%sPEnum+1))
  real(kind=8)  :: r_buff(com_nod2D%rptr(com_nod2D%rPEnum+1))
 
! Issueing MPI_IRECV first may allow the MPI_SEND to write directly into the
! receiving buffer, without a copy to an intermediate internal buffer. 
  rn = com_nod2D%rPEnum
  do n=1,rn
     source=com_nod2D%rPE(n)

     count=com_nod2D%rptr(n+1) - com_nod2D%rptr(n)

     call MPI_IRECV(r_buff(com_nod2D%rptr(n)), count, MPI_DOUBLE_PRECISION, source, source, &
         MPI_COMM_FESOM, rreq(n), MPIerr) 
  end do

  ! Put data to be communicated into send buffer
  sn = com_nod2D%sPEnum
  nend = com_nod2D%sptr(sn+1) - 1
  s_buff(1:nend) = arr2d(com_nod2D%slist(1:nend))

  do n=1, sn
     count=com_nod2D%sptr(n+1) - com_nod2D%sptr(n)

     call MPI_ISEND(s_buff(com_nod2D%sptr(n)), count, MPI_DOUBLE_PRECISION, com_nod2D%sPE(n), &
          mype, MPI_COMM_FESOM, sreq(n), MPIerr)
  end do

  call MPI_WAITALL(sn,rreq,MPI_STATUSES_IGNORE, MPIerr)

  ! Put received data to their destination

  nend = com_nod2D%rptr(rn+1) - 1
  arr2d(com_nod2D%rlist(1:nend)) = r_buff(1:nend)

  call MPI_WAITALL(sn,sreq,MPI_STATUSES_IGNORE, MPIerr)

end subroutine com_2D
!
!===================================================================
!
subroutine com_3D(arr3d)
  use o_MESH
  use o_ARRAY
  use g_PARFE 
  implicit none
  !
  real(kind=8), intent(inout)	:: arr3d(ToDim_nod3D)

  integer  	:: sreq(com_nod3D%sPEnum)
  integer  	:: rreq(com_nod3D%rPEnum)
  integer  	:: n, sn, rn, nend
  integer	:: count, source
  real(kind=8)  :: s_buff(com_nod3D%sptr(com_nod3D%sPEnum+1))
  real(kind=8)  :: r_buff(com_nod3D%rptr(com_nod3D%rPEnum+1))


  rn = com_nod3D%rPEnum

 do n=1,rn
     source = com_nod3D%rPE(n)

     count = com_nod3D%rptr(n+1) - com_nod3D%rptr(n)

     call MPI_IRECV(r_buff(com_nod3D%rptr(n)), count, MPI_DOUBLE_PRECISION, source, &
          source, MPI_COMM_FESOM, rreq(n), MPIerr) 
  end do

  ! Put data to be communicated into send buffer 
  sn = com_nod3D%sPEnum   
  nend=com_nod3D%sptr(sn+1) - 1
  s_buff(1:nend) = arr3d(com_nod3D%slist(1:nend))

  do n=1, sn
     count = com_nod3D%sptr(n+1) - com_nod3D%sptr(n)

     call MPI_ISEND(s_buff(com_nod3D%sptr(n)), count, MPI_DOUBLE_PRECISION, com_nod3D%sPE(n), mype, & 
          MPI_COMM_FESOM, sreq(n), MPIerr)
  enddo
 

  call MPI_WAITALL(rn,rreq,MPI_STATUSES_IGNORE, MPIerr)

  ! Put received data to their destination

  nend = com_nod3D%rptr(rn+1) - 1
  arr3d(com_nod3D%rlist(1:nend)) = r_buff(1:nend)

  call MPI_WAITALL(sn,sreq,MPI_STATUSES_IGNORE, MPIerr)

end subroutine com_3D
!===================================================================
!
subroutine broadcast3D(arr3D, arr3Dglobal)
  ! Makes nodal information available to all PE
  ! arr3d is any array like TF or SF of local size.
  ! arr3Dglobal is an array of nod3D size which 
  ! should be allocated before calling this routine.
  ! It will be filled with information on other PE in 
  ! natural numbering. The routine can be used to organize
  ! output in the same way as in global memory setup  
  !===================================================================

  use g_PARFE
  use o_MESH
  use o_ELEMENTS

  implicit none

  integer      ::  i, n, nTS, sender, status(MPI_STATUS_SIZE)
  integer, allocatable, dimension(:) ::  isendbuf, irecvbuf

  real(kind=8) ::  arr3D(myDim_nod3D+eDim_nod3D)
  real(kind=8) ::  arr3Dglobal(nod3D)
  real(kind=8), allocatable, dimension(:) ::  sendbuf, recvbuf

  call MPI_Barrier(MPI_COMM_FESOM, MPIERR)

  if ( mype == 0 ) then
     if (npes>1) then
        arr3Dglobal(myList_nod3D(1:myDim_nod3D))=arr3D(1:myDim_nod3D)
     end if
     do  n = 1, npes-1

        call MPI_RECV( nTS, 1, MPI_INTEGER, MPI_ANY_SOURCE, &
             0, MPI_COMM_FESOM, status, MPIerr )
        sender = status(MPI_SOURCE)
        allocate( recvbuf(1:nTS), irecvbuf(1:nTS) )
        call MPI_RECV( irecvbuf(1), nTS, MPI_INTEGER, sender, &
             1, MPI_COMM_FESOM, status, MPIerr )
        call MPI_RECV( recvbuf(1), nTS, MPI_DOUBLE_PRECISION, sender, &
             2, MPI_COMM_FESOM, status, MPIerr )

        do i = 1, nTS
           arr3Dglobal(irecvbuf(i)) = recvbuf(i)
        enddo
        deallocate( recvbuf, irecvbuf )

     enddo

  else

     allocate( sendbuf(1:myDim_nod3D), isendbuf(1:myDim_nod3D) )
     do n = 1, myDim_nod3D
        isendbuf(n) = myList_nod3D(n)
        sendbuf(n)  = arr3D(n)
     enddo
     call MPI_SEND( myDim_nod3D, 1, MPI_INTEGER, 0, 0, MPI_COMM_FESOM, MPIerr )
     call MPI_SEND( isendbuf(1), myDim_nod3D, MPI_INTEGER, 0, 1, &
          MPI_COMM_FESOM, MPIerr )
     call MPI_SEND( sendbuf(1), myDim_nod3D, MPI_DOUBLE_PRECISION, &
          0, 2, MPI_COMM_FESOM, MPIerr )
     deallocate( sendbuf, isendbuf )

  endif

  call MPI_BCAST( arr3Dglobal, nod3d, MPI_DOUBLE_PRECISION, 0, &
       MPI_COMM_FESOM, MPIerr)

end subroutine broadcast3D
!
!===================================================================
!
subroutine broadcast2D(arr2D, arr2Dglobal)
  ! Makes nodal information available to all PE 
  ! As the preceeding routine, but for 2D arrays
  !===================================================================

  use g_PARFE
  use o_MESH
  use o_ELEMENTS

  implicit none

  integer      ::  i, n, nTS, sender, status(MPI_STATUS_SIZE)
  integer, allocatable, dimension(:) ::  isendbuf, irecvbuf

  real(kind=8) ::  arr2D(myDim_nod2D+eDim_nod2D)
  real(kind=8) ::  arr2Dglobal(nod2D)
  real(kind=8), allocatable, dimension(:) ::  sendbuf, recvbuf

  call MPI_Barrier(MPI_COMM_FESOM, MPIERR)
  if ( mype == 0 ) then
     if (npes>1) then
        arr2Dglobal(myList_nod2D(1:myDim_nod2D))=arr2D(1:myDim_nod2D)
     end if
     do  n = 1, npes-1

        call MPI_RECV( nTS, 1, MPI_INTEGER, MPI_ANY_SOURCE, &
             0, MPI_COMM_FESOM, status, MPIerr )
        sender = status(MPI_SOURCE)
        allocate( recvbuf(1:nTS), irecvbuf(1:nTS) )
        call MPI_RECV( irecvbuf(1), nTS, MPI_INTEGER, sender, &
             1, MPI_COMM_FESOM, status, MPIerr )
        call MPI_RECV( recvbuf(1), nTS, MPI_DOUBLE_PRECISION, sender, &
             2, MPI_COMM_FESOM, status, MPIerr )

        do i = 1, nTS
           arr2Dglobal(irecvbuf(i)) = recvbuf(i)
        enddo
        deallocate( recvbuf, irecvbuf )

     enddo

  else

     allocate( sendbuf(1:myDim_nod2D), isendbuf(1:myDim_nod2D) )
     do n = 1, myDim_nod2D
        isendbuf(n) = myList_nod2D(n)
        sendbuf(n)  = arr2D(n)
     enddo
     call MPI_SEND( myDim_nod2D, 1, MPI_INTEGER, 0, 0, MPI_COMM_FESOM, MPIerr )
     call MPI_SEND( isendbuf(1), myDim_nod2D, MPI_INTEGER, 0, 1, &
          MPI_COMM_FESOM, MPIerr )
     call MPI_SEND( sendbuf(1), myDim_nod2D, MPI_DOUBLE_PRECISION, &
          0, 2, MPI_COMM_FESOM, MPIerr )
     deallocate( sendbuf, isendbuf )

  endif

  call MPI_BCAST( arr2Dglobal, nod2d, MPI_DOUBLE_PRECISION, 0, &
       MPI_COMM_FESOM, MPIerr)

end subroutine broadcast2D
!===================================================================
!===================================================================
!
!===================================================================
subroutine gather3D(arr3D, arr3Dglobal)
  ! Makes nodal information available to _MASTER_ PE 
  !===================================================================

  use g_PARFE
  use o_MESH
  use o_ELEMENTS

  implicit none

  integer  ::  i, n, n3D

  real(kind=8), intent(in)  :: arr3D(myDim_nod3D+eDim_nod3D)
  real(kind=8), intent(out) :: arr3Dglobal(nod3D)
  real(kind=8)              :: recvbuf(nod3D)
  integer                   :: req(npes-1)


  call MPI_Barrier(MPI_COMM_FESOM, MPIERR)
  if ( mype == 0 ) then
     
     do  n = 1, npes-1
        n3D    = remPtr_nod3D(n+1) - remPtr_nod3D(n) 
        call MPI_IRECV(recvbuf(remPtr_nod3D(n)), n3D, MPI_DOUBLE_PRECISION, n, 2, MPI_COMM_FESOM, req(n), MPIerr)
     enddo

     arr3Dglobal(myList_nod3D(1:myDim_nod3D)) = arr3D(1:myDim_nod3D)

    call MPI_WAITALL(npes-1, req, MPI_STATUSES_IGNORE, MPIerr)

    arr3Dglobal(remList_nod3D(1 : remPtr_nod3D(npes)-1)) &
                    = recvbuf(1 : remPtr_nod3D(npes)-1)

  else

     call MPI_SEND( arr3D, myDim_nod3D, MPI_DOUBLE_PRECISION, 0, 2, MPI_COMM_FESOM, MPIerr )

  endif


end subroutine gather3D
!===================================================================
!
subroutine gather2D(arr2D, arr2Dglobal)
  ! Makes nodal information available to _MASTER_ PE 
  ! As the preceeding routine, but for 2D arrays
  !===================================================================

  use g_PARFE
  use o_MESH
  use o_ELEMENTS

  implicit none

  integer  ::  i, n, n2D

  real(kind=8), intent(in)  :: arr2D(myDim_nod2D+eDim_nod2D)
  real(kind=8), intent(out) :: arr2Dglobal(nod2D)
  real(kind=8)              :: recvbuf(nod2D)
  integer                   :: req(npes-1)


  call MPI_Barrier(MPI_COMM_FESOM, MPIERR)
  if ( mype == 0 ) then
     
     do  n = 1, npes-1
        n2D    = remPtr_nod2D(n+1) - remPtr_nod2D(n) 
        call MPI_IRECV(recvbuf(remPtr_nod2D(n)), n2D, MPI_DOUBLE_PRECISION, n, 2, MPI_COMM_FESOM, req(n), MPIerr )
     enddo

     arr2Dglobal(myList_nod2D(1:myDim_nod2D)) = arr2D(1:myDim_nod2D)

     call MPI_WAITALL(npes-1, req,MPI_STATUSES_IGNORE, MPIerr)

     arr2Dglobal(remList_nod2D(1 : remPtr_nod2D(npes)-1)) &
                     = recvbuf(1 : remPtr_nod2D(npes)-1)

  else

     call MPI_SEND( arr2D, myDim_nod2D, MPI_DOUBLE_PRECISION, 0, 2, MPI_COMM_FESOM, MPIerr )

  endif


end subroutine gather2D
!===================================================================
subroutine gather3D_real4(arr3D, arr3Dglobal)
  ! Makes nodal information available to _MASTER_ PE 
  ! As the preceeding routine, but for 3D arrays again, and a local
  ! conversion from 8-Byte to 4-Byte float for netcdf output.
  !===================================================================

  use g_PARFE
  use o_MESH
  use o_ELEMENTS

  implicit none

  integer  ::  i, n,  n3D

  real(kind=8), intent(in)  :: arr3D(myDim_nod3D+eDim_nod3D)
  real(kind=4), intent(out) :: arr3Dglobal(nod3D)
  real(kind=4)              :: sendbuf(myDim_nod3D)
  real(kind=4)              :: recvbuf(nod3D)
  integer                   :: req(npes-1)


  call MPI_Barrier(MPI_COMM_FESOM, MPIerr)

  if ( mype == 0 ) then
     
     do  n = 1, npes-1
        n3D    = remPtr_nod3D(n+1) - remPtr_nod3D(n) 
        call MPI_IRECV(recvbuf(remPtr_nod3D(n)), n3D, MPI_REAL, n, 2, MPI_COMM_FESOM, req(n), MPIerr )
     enddo

     arr3Dglobal(myList_nod3D(1:myDim_nod3D)) = real(arr3D(1:myDim_nod3D),4)

    call MPI_WAITALL(npes-1, req, MPI_STATUSES_IGNORE, MPIerr)

    arr3Dglobal(remList_nod3D( 1: remPtr_nod3D(npes)-1)) &
                    = recvbuf( 1: remPtr_nod3D(npes)-1)

  else
     sendbuf(1:myDim_nod3D) = real(arr3D(1:myDim_nod3D),4)

     call MPI_SEND( sendbuf, myDim_nod3D, MPI_REAL, 0, 2, MPI_COMM_FESOM, MPIerr )

  endif


end subroutine gather3D_real4


subroutine gather3D_real4_at_rank(arr3D, arr3Dglobal, root)
  use g_PARFE, only : myList_nod3D, remList_nod3D, eDim_nod3D, myDim_nod3D, MPI_COMM_FESOM, MPI_REAL, MPI_STATUSES_IGNORE, MPIerr, mype, npes, remPtr_nod3D, MPI_REQUEST_NULL, rank0List_nod3D
  use o_MESH, only : nod3D
  implicit none
  integer  ::  remote_rank, remote_node_count
  real(kind=8), intent(in)  :: arr3D(myDim_nod3D+eDim_nod3D)
  real(kind=4), intent(out) :: arr3Dglobal(nod3D)
  real(kind=4)              :: sendbuf(myDim_nod3D)
  real(kind=4)              :: recvbuf(nod3D)
  integer                   :: req(npes-1)
  integer :: request_index
  integer, intent(in) :: root ! rank of receiving process
  integer :: i
  integer :: rank0Dim_nod3D

  !call MPI_Barrier(MPI_COMM_FESOM, MPIERR)

  if ( mype == root ) then
    request_index = 1
    rank0Dim_nod3D = nod3D - remPtr_nod3D(npes) +1
    do remote_rank = 0, npes-1
      if(remote_rank == root) cycle
      if(remote_rank == 0) then
        remote_node_count = rank0Dim_nod3D
      else
        remote_node_count = remPtr_nod3D(remote_rank+1) - remPtr_nod3D(remote_rank)
      endif
      if(remote_rank == 0) then
        call MPI_IRECV(recvbuf(remPtr_nod3D(npes)), remote_node_count, MPI_REAL, remote_rank, 2, MPI_COMM_FESOM, req(request_index), MPIerr)
      else
        call MPI_IRECV(recvbuf(remPtr_nod3D(remote_rank)), remote_node_count, MPI_REAL, remote_rank, 2, MPI_COMM_FESOM, req(request_index), MPIerr)
      endif
      request_index = request_index + 1
    enddo
    
    call MPI_WAITALL(size(req), req, MPI_STATUSES_IGNORE, MPIerr)    
    
    do remote_rank = 0, npes-1
      if(remote_rank == root) then
         arr3Dglobal(myList_nod3D(1:myDim_nod3D)) = real(arr3D(1:myDim_nod3D),4) ! local data
      else if(remote_rank == 0) then
        arr3Dglobal(rank0List_nod3D(1:rank0Dim_nod3D)) = recvbuf(remPtr_nod3D(npes):remPtr_nod3D(npes)+rank0Dim_nod3D-1) ! rank 0 data
      else
        arr3Dglobal(remList_nod3D(remPtr_nod3D(remote_rank):remPtr_nod3D(remote_rank+1)-1)) = recvbuf(remPtr_nod3D(remote_rank):remPtr_nod3D(remote_rank+1)-1) ! data of any rank!=0
      endif
    enddo
  else
     sendbuf(1:myDim_nod3D) = real(arr3D(1:myDim_nod3D),4)
     call MPI_SEND(sendbuf, myDim_nod3D, MPI_REAL, root, 2, MPI_COMM_FESOM, MPIerr)
  endif
end subroutine


subroutine gather2D_real4_at_rank(arr2D, arr2Dglobal, root)
  use g_PARFE, only : myList_nod2D, remList_nod2D, eDim_nod2D, myDim_nod2D, MPI_COMM_FESOM, MPI_REAL, MPI_STATUSES_IGNORE, MPIerr, mype, npes, remPtr_nod2D, MPI_REQUEST_NULL, rank0List_nod2D
  use o_MESH, only : nod2d
  implicit none
  integer  ::  remote_rank, remote_node_count
  real(kind=8), intent(in)  :: arr2D(myDim_nod2D+eDim_nod2D)
  real(kind=4), intent(out) :: arr2Dglobal(nod2D)
  real(kind=4) :: arr2Dglobalx(nod2D)
  real(kind=4)              :: sendbuf(myDim_nod2D)
  real(kind=4)              :: recvbuf(nod2D)
  integer                   :: req(npes-1)
  integer :: request_index
  integer, intent(in) :: root ! rank of receiving process
  integer :: i
  integer :: rank0Dim_nod2D

  !call MPI_Barrier(MPI_COMM_FESOM, MPIERR)

  if ( mype == root ) then
    request_index = 1
    rank0Dim_nod2D = nod2D - remPtr_nod2D(npes) +1
    do remote_rank = 0, npes-1
      if(remote_rank == root) cycle
      if(remote_rank == 0) then
        remote_node_count = rank0Dim_nod2D
      else
        remote_node_count = remPtr_nod2D(remote_rank+1) - remPtr_nod2D(remote_rank)
      endif
      if(remote_rank == 0) then
        call MPI_IRECV(recvbuf(remPtr_nod2D(npes)), remote_node_count, MPI_REAL, remote_rank, 2, MPI_COMM_FESOM, req(request_index), MPIerr)
      else
        call MPI_IRECV(recvbuf(remPtr_nod2D(remote_rank)), remote_node_count, MPI_REAL, remote_rank, 2, MPI_COMM_FESOM, req(request_index), MPIerr)
      endif
      request_index = request_index + 1
    enddo
    
    call MPI_WAITALL(size(req), req, MPI_STATUSES_IGNORE, MPIerr)    
    
    do remote_rank = 0, npes-1
      if(remote_rank == root) then
         arr2Dglobal(myList_nod2D(1:myDim_nod2D)) = real(arr2D(1:myDim_nod2D),4) ! local data
      else if(remote_rank == 0) then
        arr2Dglobal(rank0List_nod2D(1:rank0Dim_nod2D)) = recvbuf(remPtr_nod2D(npes):remPtr_nod2D(npes)+rank0Dim_nod2D-1) ! rank 0 data
      else
        arr2Dglobal(remList_nod2D(remPtr_nod2D(remote_rank):remPtr_nod2D(remote_rank+1)-1)) = recvbuf(remPtr_nod2D(remote_rank):remPtr_nod2D(remote_rank+1)-1) ! data of any rank!=0
      endif
    enddo
  else
     sendbuf(1:myDim_nod2D) = real(arr2D(1:myDim_nod2D),4)
     call MPI_SEND(sendbuf, myDim_nod2D, MPI_REAL, root, 2, MPI_COMM_FESOM, MPIerr)
  endif
end subroutine
