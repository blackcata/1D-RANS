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
              ONLY : Ny, dy, nu, Se, B0, C1, Ce1, Ce2, Cm, Ct, alpha, beta, mode

            USE RANS_module,                                                    &
              ONLY : k_new, dis, dis_new, U_new, nu_T, prod, Rt, fw

            IMPLICIT NONE
            INTEGER :: j
            REAL(KIND=8) :: T, f2, S, eta, Ce1_s, eta0
            REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: a,b,c,r,x

            ALLOCATE ( a(0:Ny), b(0:Ny), c(0:Ny), r(0:Ny), x(0:Ny))

            !-----------------------------------------------------------!
            !                     Set TDMA constants
            !-----------------------------------------------------------!
            IF (mode == 4) THEN
              Do j = 1,Ny-1
                T    = sqrt( (k_new(j)/dis(j))**2. + Ct**2.*(nu/dis(j)) )
                f2   = 1 - 2./9.*exp(-0.33*Rt(j)**0.5)
                S    = sqrt(2*((U_new(j+1)-U_new(j-1))/dy)**2 )
                eta  = S*k_new(j)/dis(j)
                eta0 = sqrt((Ce2-1)/(Cm*(Ce1-1)))
                Ce1_s= Ce1 - (eta*(1-eta/eta0))/(1+B0*eta**3)

                r(j) = 2.*dy**2.*Se*( ( Ce2*f2*dis(j) - Ce1_s*prod(j) ) / T     &
                     - (C1*(1-fw(j))*nu*nu_T(j)*                                &
                       ((U_new(j+1)-2*U_new(j)+U_new(j-1))/(dy**2))) )
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
