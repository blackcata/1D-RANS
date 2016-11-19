!-----------------------------------------------------------------------------------!
!
!   PROGRAM : RANS_getnut.f90
!
!   PURPOSE : Get Nu_t(turbulent kinematic viscosity) using k-e model relation
!
!                                                             2016.11.17 K.Noh
!
!-----------------------------------------------------------------------------------!

        SUBROUTINE GETNUT

            USE RANS_module,                                                    &
              ONLY : Ny, Cm, k, dis, nu_T

            IMPLICIT NONE
            INTEGER :: i,j

            DO j = 0,Ny
              nu_T(j) = Cm * k(j)**2 * dis(j)
            END DO

        END SUBROUTINE GETNUT
