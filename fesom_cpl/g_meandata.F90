module g_meandata
  implicit none
  private
  public :: Meandata, Meandata2D, Meandata3D, LocalMeandata, Snapshotdata2D, all_meandata
  type Meandata
    private
    character(:), public, allocatable :: name
    character(:), allocatable :: description
    character(:), allocatable :: unit
    character(:), allocatable :: dim_name
    real(kind=8), public, allocatable, dimension(:)  :: local_values
    real(kind=4), allocatable, dimension(:)  :: global_values
    integer, public :: addcounter
    integer :: rootrank
    contains
      procedure, pass :: init
      procedure, pass :: init_w_dimname
      procedure, pass :: var_size
      procedure, pass :: gather_data
      procedure, pass :: set_mean_output_file
      procedure, pass :: append_mean_output
      procedure, pass :: finalize
  end type

  type MeandataPtr
    class(Meandata), pointer :: ptr
  end type
  type(MeandataPtr), allocatable :: all_meandata(:)

  type, extends(Meandata) :: LocalMeandata
    contains
      procedure, pass :: init => init_local
      procedure, pass :: var_size => var_size_local
      procedure, pass :: gather_data => gather_data_local
  end type

  type, extends(Meandata) :: Meandata2D
    contains
      procedure, pass :: init => init_2D
      procedure, pass :: var_size => var_size_2D
      procedure, pass :: gather_data => gather_data_2D
  end type

  type, extends(Meandata) :: Meandata3D
    contains
      procedure, pass :: init => init_3D
      procedure, pass :: var_size => var_size_3D
      procedure, pass :: gather_data => gather_data_3D
      procedure, pass :: set_mean_output_file_with_ocean_levels
  end type

  type, extends(Meandata2D) :: Snapshotdata2D
    contains
      procedure, pass :: append_mean_output => append_mean_output_snapshot_2D
  end type
  
  contains

    subroutine init_w_dimname(this, local_size, name, description, unit, dim_name)
      implicit none
      class(Meandata), target :: this
      integer :: local_size
      character(len=*), intent(in) :: name, description, unit, dim_name
      type(MeandataPtr), allocatable :: tmparr(:)
      
      allocate(this%local_values(local_size))
      this%name = name
      this%description = description
      this%unit = unit
      this%dim_name = dim_name
      
      ! clean_meanarrays
      this%local_values = 0. 
      this%addcounter = 0
      
      ! add this instance to all_meandata array
      if( .not. allocated(all_meandata)) then
        allocate(all_meandata(1))
      else
        allocate( tmparr(size(all_meandata)+1) )
        tmparr(1:size(all_meandata)) = all_meandata
        deallocate(all_meandata)
        call move_alloc(tmparr, all_meandata)
      end if
      all_meandata(size(all_meandata))%ptr => this
    end subroutine


    subroutine init(this, local_size, name, description, unit)
      implicit none
      class(Meandata), target :: this
      integer :: local_size
      character(len=*), intent(in) :: name, description, unit
      stop 'this base class is not intended to be used, how to implement virtual functions in Fortran?'
    end subroutine


    subroutine init_local(this, local_size, name, description, unit)
      use mpi_topology_module, only : mpi_topology
      use g_parfe, only : MPI_COMM_FESOM
      implicit none
      class(LocalMeandata), target :: this
      integer :: local_size
      character(len=*), intent(in) :: name, description, unit

      this%rootrank = mpi_topology%next_host_head_rank(MPI_COMM_FESOM)
    
      call init_w_dimname(this, local_size, name, description, unit, 'nodes')
      call this%set_mean_output_file
    end subroutine


    subroutine init_2D(this, local_size, name, description, unit)
      use mpi_topology_module, only : mpi_topology
      use g_parfe, only : MPI_COMM_FESOM
      implicit none
      class(Meandata2D), target :: this
      integer :: local_size
      character(len=*), intent(in) :: name, description, unit

      this%rootrank = mpi_topology%next_host_head_rank(MPI_COMM_FESOM)
      
      call init_w_dimname(this, local_size, name, description, unit, 'nodes_2d')
      call this%set_mean_output_file
    end subroutine


    subroutine init_3D(this, local_size, name, description, unit)
      use mpi_topology_module, only : mpi_topology
      use g_parfe, only : MPI_COMM_FESOM
      implicit none
      class(Meandata3D), target :: this
      integer :: local_size
      character(len=*), intent(in) :: name, description, unit

      this%rootrank = mpi_topology%next_host_head_rank(MPI_COMM_FESOM)

      call init_w_dimname(this, local_size, name, description, unit, 'nodes_3d')
      call this%set_mean_output_file_with_ocean_levels
    end subroutine


    function var_size(this) result(s)
      implicit none
      class(Meandata) :: this
      integer :: s
      stop 'this base class is not intended to be used, how to implement virtual functions in Fortran?'
    end function


    subroutine gather_data(this, arrglobal)
      implicit none
      class(Meandata) :: this
      real(kind=4), dimension(:), intent(out) :: arrglobal
      stop 'this base class is not intended to be used, how to implement virtual functions in Fortran?'
    end subroutine


    subroutine set_mean_output_file(this)
      use g_parfe, only : mype
      use g_clock, only : ResultPath, runid, yearnew, month, day_in_month, year_start, month_start, day_start, include_fleapyear
      implicit none
      class(Meandata) :: this
      character(2000) :: filename
      character(len=8) :: firstdate
      character(len=8) :: schedule_unit
      character(len=100) :: schedule_info_txt, units_txt, calendar_txt
      integer :: schedule_first, schedule_rate

      if(mype==this%rootrank) then
        write(firstdate,'(I4,I2.2,I2.2)') yearnew, month, day_in_month
        filename=trim(ResultPath)//this%name//'_'//runid//'_'//firstdate//'.nc' 
        ! query information about the active schedule for this variable
        call schedule_info(this%name//CHAR(0), schedule_unit, schedule_first, schedule_rate)
        write(schedule_info_txt,'(a,a,a,i0,a,i0)') 'unit: ', trim(schedule_unit), ' first: ', schedule_first, ' rate: ', schedule_rate
        write(units_txt, '(a,i4.4,a,i2.2,a,i2.2,a)') 'seconds since ', year_start, '-', month_start, '-', day_start, ' 0:0:0'
        if (include_fleapyear) then
          calendar_txt = 'standard'
        else
          calendar_txt = 'noleap'
        end if
        call set_output_file(trim(filename)//CHAR(0), this%name//CHAR(0), this%var_size(), this%description//CHAR(0), this%unit//CHAR(0), this%dim_name//CHAR(0), 'output_schedule'//CHAR(0), trim(schedule_info_txt)//CHAR(0), trim(units_txt)//CHAR(0), trim(calendar_txt)//CHAR(0))
      end if
    end subroutine


    subroutine set_mean_output_file_with_ocean_levels(this)
      use g_parfe, only : mype
      use g_clock, only : ResultPath, runid, yearnew, month, day_in_month, year_start, month_start, day_start, include_fleapyear
      use o_mesh, only : nod2D
      use ocean_mesh_module, only : ocean_mesh__wetlevelsize
      use g_config, only : levelwise_output
      implicit none
      class(Meandata3D) :: this
      character(2000) :: filename
      character(len=8) :: firstdate
      character(len=8) :: schedule_unit
      character(len=100) :: schedule_info_txt, units_txt, calendar_txt
      integer :: schedule_first, schedule_rate

      if(mype==this%rootrank) then
        write(firstdate,'(I4,I2.2,I2.2)') yearnew, month, day_in_month
        filename=trim(ResultPath)//this%name//'_'//runid//'_'//firstdate//'.nc' 
        ! query information about the active schedule for this variable
        call schedule_info(this%name//CHAR(0), schedule_unit, schedule_first, schedule_rate)
        write(schedule_info_txt,'(a,a,a,i0,a,i0)') 'unit: ', trim(schedule_unit), ' first: ', schedule_first, ' rate: ', schedule_rate
        write(units_txt, '(a,i4.4,a,i2.2,a,i2.2,a)') 'seconds since ', year_start, '-', month_start, '-', day_start, ' 0:0:0'
        if (include_fleapyear) then
          calendar_txt = 'standard'
        else
          calendar_txt = 'noleap'
        end if
        if(levelwise_output) then
          call set_ocean_levels_output_file(trim(filename)//CHAR(0), this%name//CHAR(0), nod2d, ocean_mesh__wetlevelsize, this%description//CHAR(0), this%unit//CHAR(0), this%dim_name//CHAR(0), 'output_schedule'//CHAR(0), trim(schedule_info_txt)//CHAR(0), trim(units_txt)//CHAR(0), trim(calendar_txt)//CHAR(0))
        else
          call set_output_file(trim(filename)//CHAR(0), this%name//CHAR(0), this%var_size(), this%description//CHAR(0), this%unit//CHAR(0), this%dim_name//CHAR(0), 'output_schedule'//CHAR(0), trim(schedule_info_txt)//CHAR(0), trim(units_txt)//CHAR(0), trim(calendar_txt)//CHAR(0))
        end if
      end if
    end subroutine
    
    
    subroutine append_mean_output(this)
      use g_clock, only : dt, istep
      use g_parfe, only : mype
      implicit none
      class(Meandata) :: this
    
      this%local_values = this%local_values /float(this%addcounter) ! compute_means

      ! gather mean data and dump from single rank
      if(mype==this%rootrank) then
        if(.not. allocated(this%global_values)) then
          allocate(this%global_values(this%var_size()))
        end if
        call end_append_output(this%name//CHAR(0))
      else
        if(.not. allocated(this%global_values)) allocate(this%global_values(0))
      end if
      
      call this%gather_data(this%global_values)

      if(mype==this%rootrank) then        
        call begin_append_output(this%name//CHAR(0), this%var_size(), this%global_values, dt*istep)
      end if
      
      this%local_values = 0. ! clean_meanarrays
      this%addcounter = 0 ! clean_meanarrays
    end subroutine    


    ! write output without dividing by number of accumulations (addcounter)
    subroutine append_mean_output_snapshot_2D(this)
      use g_clock, only : dt, istep
      use g_parfe, only : mype
      implicit none
      class(Snapshotdata2D) :: this

      ! gather mean data and dump from single rank
      if(mype==this%rootrank) then
        if(.not. allocated(this%global_values)) then
          allocate(this%global_values(this%var_size()))
        end if
        call end_append_output(this%name//CHAR(0))
      else
        if(.not. allocated(this%global_values)) allocate(this%global_values(0))
      end if
      
      call this%gather_data(this%global_values)

      if(mype==this%rootrank) then        
        call begin_append_output(this%name//CHAR(0), this%var_size(), this%global_values, dt*istep)
      end if
      
      this%local_values = 0. ! clean_meanarrays
      this%addcounter = 0 ! clean_meanarrays
    end subroutine    


    subroutine finalize(this)
      use g_parfe, only : mype
      implicit none
      class(Meandata) :: this
      
      if(mype==this%rootrank) then
        call end_append_output(this%name//CHAR(0))
        if(allocated(this%global_values)) deallocate(this%global_values)
      else
        if(allocated(this%global_values)) deallocate(this%global_values)
      end if
    end subroutine    


    function var_size_local(this) result(s)
      use o_mesh, only : nod2d
      implicit none
      class(LocalMeandata) :: this
      integer :: s
      s = size(this%local_values)
    end function


    function var_size_2D(this) result(s)
      use o_mesh, only : nod2d
      implicit none
      class(Meandata2D) :: this
      integer :: s
      s = nod2d
    end function


    function var_size_3D(this) result(s)
      use o_mesh, only : nod3d
      implicit none
      class(Meandata3D) :: this
      integer :: s
      s = nod3d
    end function

    
    subroutine gather_data_local(this, arrglobal)
      implicit none
      class(LocalMeandata) :: this
      real(kind=4), dimension(:), intent(out) :: arrglobal
      arrglobal = this%local_values
    end subroutine


    subroutine gather_data_2D(this, arrglobal)
      implicit none
      class(Meandata2D) :: this
      real(kind=4), dimension(:), intent(out) :: arrglobal
      call gather2D_real4_at_rank(this%local_values, arrglobal, this%rootrank)
    end subroutine


    subroutine gather_data_3D(this, arrglobal)
      implicit none
      class(Meandata3D) :: this
      real(kind=4), dimension(:), intent(out) :: arrglobal
      call gather3D_real4_at_rank(this%local_values, arrglobal, this%rootrank)
    end subroutine
    
end module g_meandata
