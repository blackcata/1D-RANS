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
              ONLY : Ny, Cm, k, dis, nu_T, fm, nu

            IMPLICIT NONE
            INTEGER :: i,j
            REAL(KIND=8) :: Rt

            DO j = 0,Ny
              Rt      = k(j)**2 /(nu * dis(j))
              fm(j)   = exp(-3.4/(1+Rt/50)**2)
              nu_T(j) = Cm * fm(j)* k(j)**2 / dis(j)
            END DO

        END SUBROUTINE GETNUT
