module ocean_mesh_module
  
  implicit none
  public ocean_mesh__depths, ocean_mesh__add_depth, ocean_mesh__add_wetrange, ocean_mesh__trim_empty_wetranges, ocean_mesh__purge, ocean_mesh__global_wetranges, ocean_mesh__wetlevelsize
  private

  real(kind=8), protected, save, allocatable, dimension(:) :: ocean_mesh__depths
  integer, protected, save :: ocean_mesh__wetlevelsize = 0

  type Range
    integer :: index, size
  end type Range
  type Ranges
    type(Range), allocatable :: ranges(:)
  end type Ranges
  type(Ranges), protected, allocatable :: ocean_mesh__global_wetranges(:) ! wet node index/size pairs per level
  
contains


  subroutine ocean_mesh__purge
    if(allocated(ocean_mesh__depths)) deallocate(ocean_mesh__depths)
    if(allocated(ocean_mesh__global_wetranges)) deallocate(ocean_mesh__global_wetranges)
  end subroutine


  subroutine ocean_mesh__add_depth(d)
    real(kind=8) d
    real(kind=8), allocatable, dimension(:) :: tmparr
    
    if( .not. allocated(ocean_mesh__depths)) then
      allocate(ocean_mesh__depths(1))
    else
      allocate( tmparr(size(ocean_mesh__depths)+1) )
      tmparr(1:size(ocean_mesh__depths)) = ocean_mesh__depths
      deallocate(ocean_mesh__depths)
      call move_alloc(tmparr, ocean_mesh__depths)
    end if
    
    ocean_mesh__depths(size(ocean_mesh__depths)) = d
  end subroutine
  
  
  subroutine ocean_mesh__add_wetrange(level, index, size_)
    integer level
    integer index
    integer size_
    type(Range), allocatable :: tmpranges(:)
        
    call ensure_level_allocated(level)

    ! add another range for this level    
    allocate( tmpranges(size(ocean_mesh__global_wetranges(level)%ranges)+1) )
    tmpranges(1:size(ocean_mesh__global_wetranges(level)%ranges)) = ocean_mesh__global_wetranges(level)%ranges
    deallocate(ocean_mesh__global_wetranges(level)%ranges)
    call move_alloc(tmpranges, ocean_mesh__global_wetranges(level)%ranges)
    
    ocean_mesh__global_wetranges(level)%ranges(size( ocean_mesh__global_wetranges(level)%ranges ))%index = index
    ocean_mesh__global_wetranges(level)%ranges(size( ocean_mesh__global_wetranges(level)%ranges ))%size = size_
  end subroutine


  ! remove any empty levels from the end of wetranges (i.e. levels without any ocean nodes)
  subroutine ocean_mesh__trim_empty_wetranges
    type(Ranges), allocatable :: tmpranges(:)
    
    if(allocated(ocean_mesh__global_wetranges)) then
      do while(size(ocean_mesh__global_wetranges(size(ocean_mesh__global_wetranges))%ranges) == 0)
        allocate( tmpranges(size(ocean_mesh__global_wetranges)-1) )
        tmpranges = ocean_mesh__global_wetranges(1:size(ocean_mesh__global_wetranges)-1)
        deallocate(ocean_mesh__global_wetranges)
        call move_alloc(tmpranges, ocean_mesh__global_wetranges)
      end do

      ocean_mesh__wetlevelsize = size(ocean_mesh__global_wetranges)
    end if
  end subroutine
    
  
  subroutine ensure_level_allocated(level)
    type(Ranges), allocatable :: tmparr(:)
    integer level
    integer old_size
    integer i
    
    if( .not. allocated(ocean_mesh__global_wetranges)) then
      old_size = 0
      allocate(ocean_mesh__global_wetranges(level))
    else if( size(ocean_mesh__global_wetranges) < level ) then
      old_size = size(ocean_mesh__global_wetranges)
      allocate( tmparr(level) )
      tmparr(1:size(ocean_mesh__global_wetranges)) = ocean_mesh__global_wetranges
      deallocate(ocean_mesh__global_wetranges)
      call move_alloc(tmparr, ocean_mesh__global_wetranges)
    else
      old_size = size(ocean_mesh__global_wetranges)
    end if

    do i=old_size+1, size(ocean_mesh__global_wetranges)
      allocate( ocean_mesh__global_wetranges(i)%ranges(0) )
    end do

    ocean_mesh__wetlevelsize = size(ocean_mesh__global_wetranges)
  end subroutine
   
end module
