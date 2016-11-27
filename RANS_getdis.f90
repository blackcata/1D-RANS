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
              ONLY : Ny, dy, nu, Se, Ce1, Ce2, Ct, alpha, beta, mode

            USE RANS_module,                                                    &
              ONLY : k_new, dis, dis_new, U_new, nu_T, prod, Rt

            IMPLICIT NONE
            INTEGER :: j
            REAL(KIND=8) :: T, f2
            REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: a,b,c,r,x

            ALLOCATE ( a(0:Ny), b(0:Ny), c(0:Ny), r(0:Ny), x(0:Ny))

            !-----------------------------------------------------------!
            !                     Set TDMA constants
            !-----------------------------------------------------------!
            IF (mode == 4) THEN
              Do j = 1,Ny-1
                T    = sqrt( (k_new(j)/dis(j))**2. + Ct**2.*(nu/dis(j)) )
                f2   = 1 - 2./9.*exp(-0.33*Rt(j)**0.5)

                r(j) = 2.*dy**2.*Se* ( Ce2*f2*dis(j) - Ce1*prod(j) ) / T
              END DO
            ELSE
              DO j =1,Ny-1
                r(j) = 2.*dy**2.*Se* ( Ce2*dis(j)**2./k_new(j)                  &
                                     - Ce1*dis(j)/k_new(j)*prod(j) )
              END DO
            END IF

            DO j = 1,Ny-1
              a(j) =   Se*2.*nu + nu_T(j-1) + nu_T(j)
              b(j) = -(Se*4.*nu + nu_T(j+1) + 2.*nu_T(j) + nu_T(j-1))
              c(j) =   Se*2.*nu + nu_T(j+1) + nu_T(j)
            END DO

            r(0:Ny) = r(0:Ny) + b(0:Ny) * dis(0:Ny)*(1.-alpha)/alpha
            b(0:Ny) = b(0:Ny) / alpha
            x(0:Ny) = dis(0:Ny)

            !-----------------------------------------------------------!
            !                     Boundary conditions
            !-----------------------------------------------------------!
            b(0)  = 1.
            b(Ny) = 1.
            a(Ny) = -1.
            c(0)  = 0.
            r(0)  = 2*nu*k_new(1)/(dy**2.)
            r(Ny) = 0

            !-----------------------------------------------------------!
            !                       Calculate TDMA
            !-----------------------------------------------------------!
            CALL TDMA_Solver(a,b,c,r,x,Ny)

            !-----------------------------------------------------------!
            !                   Relaxation & Update
            !-----------------------------------------------------------!
            dis_new(0:Ny) = beta * x(0:Ny) + (1.-beta) * dis(0:Ny)
            dis(0:Ny) = dis_new(0:Ny)
            DEALLOCATE(a,b,c,r,x)

        END SUBROUTINE GETDIS
