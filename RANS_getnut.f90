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
              ONLY : Ny, Cm, k, dis, nu_T, fm, nu, u_tau, Y, del

            IMPLICIT NONE
            INTEGER :: i,j
            REAL(KIND=8) :: Rt, Ry
            REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: Y_tmp

            ALLOCATE( Y_tmp(0:Ny) )

            DO j = 0,Ny
               IF (j<Ny/2) THEN
                  Y_tmp(j) = Y(j)
               ELSE
                  Y_tmp(j) = 2*del - Y(j)
               END IF

              fm(j)   = ( 1 - exp(-Y_tmp(j)/(nu/u_tau)/40))**2  ! van Driest(1954)
              nu_T(j) = Cm *fm(j)* k(j)**2 / dis(j)

            END DO

        END SUBROUTINE GETNUT
