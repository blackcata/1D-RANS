!-----------------------------------------------------------------------------------!
!
!   PROGRAM : RANS_getk.f90
!
!   PURPOSE : Get k(turbulent kinetic energy) using k-e model relation
!
!                                                             2016.11.19 K.Noh
!
!-----------------------------------------------------------------------------------!

        SUBROUTINE GETK

            USE RANS_module,                                                    &
              ONLY : Ny, dy, nu, Sk, alpha, beta

            USE RANS_module,                                                    &
              ONLY : k, k_new, dis, U_new, nu_T

            IMPLICIT NONE
            INTEGER :: i,j

            REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: a,b,c,r,x

            ALLOCATE ( a(0:Ny), b(0:Ny), c(0:Ny), r(0:Ny), x(0:Ny))

            DO j = 1,Ny-1
              a(j) =   nu_T(j-1) + nu_T(j)
              b(j) = -(nu_T(j+1) + 2*nu_T(j) + nu_T(j-1))/alpha
              c(j) =   nu_T(j+1) + nu_T(j)
              r(j) = 2*dy*2*Sk* (dis(j) - nu_T(j) * (U_new(j+1) - U_new(j-1))/(2*dy**2))&
                      + a(j) * k(j) * (1-alpha) /alpha
            END DO

            x(0:Ny) = k(0:Ny)

            b(0)  = 1
            b(Ny) = 1
            a(Ny) = 0
            c(0)  = 0
            r(0)  = 0
            r(Ny) = 0

            CALL TDMA_Solver(a,b,c,r,x,Ny)

            k_new(0:Ny) = beta * x(0:Ny) + (1-beta) * k(0:Ny)
            DEALLOCATE(a,b,c,r,x)

        END SUBROUTINE GETK
