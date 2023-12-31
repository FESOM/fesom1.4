module g_clock
  !combining RT and Lars version
  !
  use g_config
  implicit none
  save
  real(kind=8)             :: timeold, timenew     !time in a day, unit: sec
  integer                  :: dayold, daynew       !day in a year
  integer                  :: yearold, yearnew     !year before and after time step
  integer                  :: month, day_in_month  !month and day in a month
  integer                  :: fleapyear            !1 fleapyear, 0 not 
  integer                  :: ndpyr                !number of days in yearnew 
  integer                  :: num_day_in_month(0:1,12)
  character(4)             :: cyearold, cyearnew   !year as character string      
  integer                  :: year_start, month_start, day_start ! will be used in the output files to adopt to NetCDF Climate and Forecast (CF) Metadata Conventions
  integer :: day_count = 1, month_count = 1 ! total number of days/months this clock has stepped

  data num_day_in_month(0,:) /31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
  data num_day_in_month(1,:) /31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/


contains
  !
  !--------------------------------------------------------------------------------
  !
  subroutine clock

    implicit none
    integer         :: i
    real(kind=8)    :: aux1, aux2
    !
    timeold=timenew 
    dayold=daynew
    yearold=yearnew

    ! update time
    timenew=timenew+dt          

    ! update day
    if (timenew>86400.) then  !assumed that time step is less than one day!
       daynew=daynew+1
       timenew=timenew-86400.
       day_count = day_count+1
    endif

    ! update year
    if (daynew>ndpyr) then
       daynew=1
       yearnew=yearnew+1
       call check_fleapyr(yearnew, fleapyear)
       ndpyr=365+fleapyear
       write(cyearold,'(i4)') yearold
       write(cyearnew,'(i4)') yearnew
    endif

    ! find month and dayinmonth at new time step
    aux1=0
    do i=1,12
       aux2=aux1+num_day_in_month(fleapyear,i)
       if(daynew>aux1 .and. daynew<=aux2) then
          if(month .ne. i) month_count = month_count + 1
          month=i
          day_in_month=daynew-aux1
          exit
       end if
       aux1=aux2
    end do

  end subroutine clock
  !
  !--------------------------------------------------------------------------------
  !
  subroutine clock_init
    use o_param
    use g_parfe
    implicit none
    integer         :: i, daystart, yearstart
    real(kind=8)    :: aux1, aux2, timestart

    ! the model inialized at
    timestart=timenew
    daystart=daynew
    yearstart=yearnew
    seconds_from_yearstart=0

    ! init clock for this run
    open(99,file=runid//'.clock',status='old')
    read(99,*) timeold, dayold, yearold
    read(99,*) timenew, daynew, yearnew
    close(99)
    if(daynew==0) daynew=1

    ! check if this is a restart or not
    if(yearnew==yearstart .and. daynew==daystart .and. timenew==timestart) then
       r_restart=.false.
       yearold=yearnew-1 !required for checking if create new output files
    else
       r_restart=.true.
    end if

    ! year as character string 
    write(cyearold,'(i4)') yearold
    write(cyearnew,'(i4)') yearnew

    ! if restart model at beginning of a day, set timenew to be zero
    if (timenew==86400.) then  
       timenew=0.0
       daynew=daynew+1
    endif

    ! set timeold to be timenew, ready for initializing forcing fields,
    ! yearold should not be updated here, which is requird to open input files.
    timeold=timenew 
    dayold=daynew

    ! check fleap year
    call check_fleapyr(yearnew, fleapyear)
    ndpyr=365+fleapyear

    ! find month and dayinmonth at the new time step
    aux1=0
    do i=1,12
       aux2=aux1+num_day_in_month(fleapyear,i)
       if(daynew>aux1 .and. daynew<=aux2) then
          month=i
          day_in_month=daynew-aux1
          exit
       end if
       aux1=aux2
    end do

    if(mype==0) then
       write(*,*)'clock initialized at time ', real(timenew,4), daynew, yearnew
       if(r_restart) then
          write(*,*) 'THIS RUN IS A RESTART RUN!'
       end if
    end if

  year_start	=yearnew
  month_start	=month
  day_start	=day_in_month
  end subroutine clock_init
  !
  !-------------------------------------------------------------------------------
  !
  subroutine clock_finish
    use o_param
    implicit none
    !
    if ((daynew==ndpyr) .and. (timenew==86400.)) then
       timenew=0.0
       daynew=1
       yearnew=yearold+1
    endif

    open(99,file=runid//'.clock',status='unknown')
    write(99,*) timeold, dayold, yearold
    write(99,*) timenew, daynew, yearnew
    close(99)
  end subroutine clock_finish
  !
  !----------------------------------------------------------------------------
  !
  subroutine clock_newyear
    implicit none
    !
    if ((daynew>=ndpyr).and.(timenew==86400.)) then
       timenew=0.0
       daynew=1
       yearnew=yearold+1
       write(cyearnew,'(i4)') yearnew
    endif
  end subroutine clock_newyear
  !
  !----------------------------------------------------------------------------
  !
  subroutine check_fleapyr(year, flag)
    use O_param
    implicit none
    integer, intent(in) :: year      
    integer, intent(out):: flag

    flag=0

    if(.not.include_fleapyear) return

    if ((mod(year,4)==0.and.mod(year,100)/=0) .or. mod(year,400)==0) then
       flag=1
    endif
  end subroutine check_fleapyr
  !
  !----------------------------------------------------------------------------
  !
end module g_clock
