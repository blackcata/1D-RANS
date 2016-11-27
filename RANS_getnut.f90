!-----------------------------------------------------------------------------------!
!
!   PROGRAM : RANS_getnut.f90
!
!   PURPOSE : Get Nu_t(turbulent kinematic viscosity) using k-e model relation
!
!                                                             2016.11.19 K.Noh
!
!-----------------------------------------------------------------------------------!

        SUBROUTINE GETNUT

            USE RANS_module,                                                    &
              ONLY : Ny, Cm, k, dis, nu_T, fm

            IMPLICIT NONE
            INTEGER :: i,j
            REAL(KIND=8) :: Rt, Ry

            DO j = 0,Ny
              nu_T(j) = Cm *fm(j)* k(j)**2. / dis(j)
              ! print*,fm(j),nu_T(j)
            END DO

        END SUBROUTINE GETNUT
