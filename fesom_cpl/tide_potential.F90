SUBROUTINE potential

USE o_MESH
USE o_ARRAY
USE o_PARAM
USE g_ROTATE_grid
USE g_config
use g_clock
  IMPLICIT NONE

  REAL(kind=8) :: DY, TY, TY2, TY3
  REAL(kind=8) :: h0, s0, p0, ta, ax, ay, lon, lat
  REAL(kind=8) :: Apt(11), Fpt(11), sxi(11)
  REAL(kind=8) :: ssh_sd, ssh_lp, ssh_d, c,time_cor !year, day_number
  INTEGER :: iz, n, j, period_m2, i

!     amplitude

  DATA Apt/0.242334,0.112841,0.046398,0.030704 & ! M2, S2, N2, K2
       ,0.141565,0.100514,0.046843,0.019256 &          ! K1, O1, P1, Q1
       ,0.041742,0.022026,0.019446/                              ! Mf, Mm, Ssa

!     phase

  DATA Fpt/1.40519d-4,1.45444d-4,1.37880d-4,1.45842d-4  &
       ,0.72921d-4,0.67598d-4,0.72523d-4,0.64959d-4          &
       ,0.53234d-5,0.26392d-5,0.03982d-5/

  !     day_number ---> is the day number of the year (day_number = 1 for January 1)
  !     year >= 1975 ---> is the year number	

 !year = 1975.0
 !day_number=FLOOR(time/(24.0*3600.0)) + 1.0
  DY =  daynew + 365.0*(yearnew - 1975.0) + floor((yearnew - 1975.0)/4.0)

  TY = (27392.500528 + 1.0000000356*DY)/36525.0

  TY2 = TY*TY
  TY3 = TY2*TY

  h0 = 279.69668 + 36000.768930485*TY + 3.03d-4*TY2
  CALL control_360(iz,h0)

  s0 = 270.434358 + 481267.88314137*TY - 0.001133*TY2 + &
       1.9d-6*TY3
  CALL control_360(iz,s0)

  p0 = 334.329653 + 4069.0340329575*TY - 0.010325*TY2 -&
       1.2d-5*TY3
  CALL control_360(iz,p0)

  sxi(1) = (h0 - s0)*2.0*rad                        ! M2
  sxi(2) = 0.0                                            ! S2
  sxi(3) = (2.0*h0 - 3.0*s0 + p0)*rad     ! N2
  sxi(4) = 2.0*h0*rad                                 ! K2

  sxi(5) = (h0 + 90.0)*rad                          ! K1
  sxi(6) = (h0 - 2.0*s0 - 90.0)*rad         ! O1
  sxi(7) = (h0 - 90.0)*rad                           ! P1
  sxi(8) = (h0 - 3.0*s0 + p0 - 90.0)*rad  ! Q1

  sxi(9) = 2.0*s0*rad                                  ! Mf
  sxi(10)= (s0 - p0)*rad                                    ! Mm
  sxi(11)= 2.0*h0*rad                                 ! Ssa 

  ! calculation equilibrium tides
  ! ssh_sd ----> semidiurnal (M2, S2, N2, K2)
  ! ssh_d  ----> diurnal     (K1, O1, P1, Q1)
  ! ssh_lp ----> long-period (Mf, Mm, Ssa)

  period_m2 = nint(12.42*3600.0/dt)
  
  ! shift of the origin of the universal standard time by the individual substitution (29.11.2011)
  ! time_cor(j) = time - sxi(j)/Fpt(j)

  DO n=1, ToDim_nod2D
     ax = coord_nod2d(1,n)
     ay = coord_nod2d(2,n)
     call r2g(lon, lat, ax, ay)   ! reserved for rotated meshes, see 3D model

     ssh_sd = 0.0
     DO j=1,4
        time_cor = dt*(real(istep)-1) - sxi(j)/Fpt(j)
        ssh_sd = ssh_sd + Apt(j)*SIN(0.5*pi - lat)*SIN(0.5*pi - lat)&
             *COS(Fpt(j)*time_cor + 2.0*lon + sxi(j))
     ENDDO
     ssh_d = 0.0
     DO j=5,8
        time_cor = dt*(real(istep)-1) - sxi(j)/Fpt(j)
            ssh_d = ssh_d + Apt(j)*SIN(2.0*(0.5*pi - lat))&
              *COS(Fpt(j)*time_cor + lon + sxi(j))
     ENDDO
     ssh_lp = 0.0
     DO j=9,11
        time_cor = dt*(real(istep)-1) - sxi(j)/Fpt(j)
                 ssh_lp = ssh_lp + Apt(j)*(3.0*SIN(0.5*pi - lat)**2.0 - 2.0)&
                      *COS(Fpt(j)*time_cor + sxi(j))
     ENDDO
!++++++++++++++++++++++++++++++++++++++++
! equilibrium tide for 11 harmonics
!++++++++++++++++++++++++++++++++++++++++

     ssh_gp(n) = ssh_sd + ssh_d + ssh_lp    !!! now, for pressure gradient (ssh - ssh_gp)
     !      c = (time/(dt*period_m2*nint(final_period/3.0_8)))
     !      if (c > 1.0_8) c = 1.0_8
  !   ssh_gp(n) = ssh_gp(n)        !*c

  ENDDO

END SUBROUTINE potential
!===========================================================================
SUBROUTINE control_360(iz,fr)

USE o_PARAM
  implicit none

  REAL(KIND=8) :: fr
  integer :: iz

  iz = INT(fr/360.0)
  IF (iz >= 1) fr = fr - float(iz)*360.0

END SUBROUTINE control_360
!===========================================================================
