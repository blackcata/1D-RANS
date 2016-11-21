!-----------------------------------------------------------------------------------!
!
!   PROGRAM : RANS_getdis.f90
!
!   PURPOSE : Get dis(dissipation) using k-e model relation
!
!                                                             2016.11.20 K.Noh
!
!-----------------------------------------------------------------------------------!

        SUBROUTINE GETDIS

            USE RANS_module,                                                    &
              ONLY : Ny, dy, nu, Se, Ce1, Ce2, K0, u_tau, alpha, beta

            USE RANS_module,                                                    &
              ONLY : k_new, dis, dis_new, U_new, nu_T, prod

            IMPLICIT NONE
            INTEGER :: i,j

            REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: a,b,c,r,x

            ALLOCATE ( a(0:Ny), b(0:Ny), c(0:Ny), r(0:Ny), x(0:Ny))

            !-----------------------------------------------------------!
            !                     Set TDMA constants
            !-----------------------------------------------------------!
            DO j = 1,Ny-1
              a(j) =   nu_T(j-1) + nu_T(j)
              b(j) = -(nu_T(j+1) + 2*nu_T(j) + nu_T(j-1))
              c(j) =   nu_T(j+1) + nu_T(j)
              r(j) = 2*dy**2*Se* ( Ce2*dis(j)**2/k_new(j)                       &
                                  - Ce1*nu_T(j)*dis(j)/k_new(j)*prod(j) )
            END DO

            x(0:Ny) = dis(0:Ny)

            !-----------------------------------------------------------!
            !                     Boundary conditions
            !-----------------------------------------------------------!
            b(0)  = 1
            b(Ny) = 1
            a(Ny) = 0
            c(0)  = 0
            r(0)  = u_tau**2*nu/(K0*dy)
            r(Ny) = u_tau**2*nu/(K0*dy)

            r(0:Ny) = r(0:Ny) + b(0:Ny) * dis(0:Ny)*(1-alpha)/alpha
            b(0:Ny) = b(0:Ny) / alpha

            !-----------------------------------------------------------!
            !                       Calculate TDMA
            !-----------------------------------------------------------!
            CALL TDMA_Solver(a,b,c,r,x,Ny)

            !-----------------------------------------------------------!
            !                   Relaxation & Update
            !-----------------------------------------------------------!
            dis_new(0:Ny) = beta * x(0:Ny) + (1-beta) * dis(0:Ny)
            dis(0:Ny) = dis_new(0:Ny)
            DEALLOCATE(a,b,c,r,x)

        END SUBROUTINE GETDIS
