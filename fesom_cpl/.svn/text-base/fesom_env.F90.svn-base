! sanitize fesom config, settings and environment
module fesom_env_module
  use g_config
  implicit none
  public fesom_env__readable, fesom_env__writable, fesom_env_mesh_partition_path, fesom_env__check_namelist_config
  private
  
contains

  subroutine fesom_env__check_namelist_config(namelistpath)
    character(len=*), optional, intent(in) :: namelistpath
    
    if(present(namelistpath)) then
      call read_namelist_config(namelistpath)
    else
      call read_namelist_config('namelist.config')
    end if
        
    call ensure_path_readable(MeshPath, 'MeshPath')
    call ensure_path_readable(ClimateDataPath, 'ClimateDataPath')
#ifndef __oasis
    call ensure_path_readable(ForcingDataPath, 'ForcingDataPath')
#endif

    call ensure_path_writable(ResultPath, 'ResultPath')
  end subroutine

 
  logical function fesom_env__readable(path) result(result)
    character(len=*), intent(in) :: path
    integer exitval
    
    call execute_command_line ("test -r "//path, exitstat=exitval) ! can not use inquire here as it is horribly unreliable
    result = (exitval==0)
  end function


  logical function fesom_env__writable(path) result(result)
    character(len=*), intent(in) :: path
    integer exitval
    
    call execute_command_line ("test -w "//path, exitstat=exitval) ! can not use inquire here as it is horribly unreliable
    result = (exitval==0)
  end function
  
  
  function fesom_env_mesh_partition_path(meshdir, partition_count) result(resultdir)
    integer, intent(in) :: partition_count
    character(len=*), intent(in) :: meshdir
    character(len=:), allocatable :: resultdir
    character(len=int(log10(real(partition_count)))+1) :: partition_count_str

    write(partition_count_str,'(i0)') partition_count
    resultdir = trim(meshdir)//'/dist/'//partition_count_str//'p'
  end function


  subroutine ensure_path_readable(path, name)
    character(len=*), intent(in) :: path
    character(len=*), intent(in), optional :: name
    character(:), allocatable :: nametxt
    
    nametxt = ''
    if(present(name)) nametxt = ' '//name
    
    if(.not. fesom_env__readable(path) .eqv. .true.) then
      print '(a)',('readable? no'//nametxt//' <'//trim(path)//'>')
      error stop ! can not pass the string to error stop, as even the currently newest ifort (18.0.1) does not support non-const strings here
    else
      print '(a)',('readable? yes'//nametxt//' <'//trim(path))//'>'
    end if
  end subroutine


  subroutine ensure_path_writable(path, name)
    character(len=*), intent(in) :: path
    character(len=*), intent(in), optional :: name
    character(:), allocatable :: nametxt
    
    nametxt = ' '
    if(present(name)) nametxt = nametxt//name
    
    if(.not. fesom_env__writable(path) .eqv. .true.) then
      print '(a)',('writable? no'//nametxt//' <'//trim(path))//'>'
      error stop ! can not pass the string to error stop, as even the currently newest ifort (18.0.1) does not support non-const strings here
    else
      print '(a)',('writable? yes'//nametxt//' <'//trim(path))//'>'
    end if
  end subroutine


  subroutine read_namelist_config(namelistpath)
    character(len=*), intent(in) :: namelistpath
    integer fileunit

    call ensure_path_readable(namelistpath)
    open(newunit=fileunit,file=namelistpath,action='read')
    read(fileunit,nml=paths)
    close(fileunit)
  end subroutine
   
end module
