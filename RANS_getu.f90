!-----------------------------------------------------------------------------------!
!
!   PROGRAM : RANS_getnut.f90
!
!   PURPOSE : Get Nu_t(turbulent kinematic viscosity) using k-e model relation
!
!                                                             2016.11.19 K.Noh
!
!-----------------------------------------------------------------------------------!

        SUBROUTINE GETU

            USE RANS_module,                                                    &
              ONLY : Ny, dy, nu, del, Re_tau, alpha, beta

            USE RANS_module,                                                    &
              ONLY : U, U_new, nu_T

            IMPLICIT NONE
            INTEGER :: i,j
            REAL(KIND=8) :: C0, u_tau

            REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: a,b,c,r,x

            ALLOCATE ( a(0:Ny), b(0:Ny), c(0:Ny), r(0:Ny), x(0:Ny))

            u_tau = Re_tau*nu / del
            C0 = u_tau**2 / del

            DO j = 1,Ny-1
              a(j) =   2*nu + nu_T(j-1) + nu_T(j)
              b(j) = -(4*nu + nu_T(j+1) + 2*nu_T(j) + nu_T(j-1))/alpha
              c(j) =   2*nu + nu_T(j+1) + nu_T(j)
              r(j) = -2*dy*dy*C0 + a(j) * U(j) * (1-alpha) /alpha
            END DO

            x(0:Ny) = U(0:Ny)

            b(0)  = 1
            b(Ny) = 1
            a(Ny) = 0
            c(0)  = 0
            r(0)  = 0
            r(Ny) = 0

            CALL TDMA_Solver(a,b,c,r,x,Ny)

            U_new(0:Ny) = beta * x(0:Ny) + (1-beta) * U(0:Ny)

        END SUBROUTINE GETU
