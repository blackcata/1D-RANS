!-----------------------------------------------------------------------------------!
!
!   PROGRAM : RANS_getfm.f90
!
!   PURPOSE : Get f_m(damping function) using k-e model relation
!
!                                                             2016.11.25 K.Noh
!
!-----------------------------------------------------------------------------------!

        SUBROUTINE GETFM

            USE RANS_module,                                                    &
              ONLY : Ny, Cm, k, dis, nu_T, fm, prod, nu, u_tau, Y, del, mode,   &
                     dy, fw, A0, A1, Cd, Cp, Ce

            IMPLICIT NONE
            INTEGER :: i,j
            REAL(KIND=8) :: L, fm1, fm2, C1
            REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: Rt,Y_tmp,a,b,c,r,x

            ALLOCATE( Y_tmp(0:Ny), Rt(0:Ny))
            ALLOCATE( a(0:Ny), b(0:Ny), c(0:Ny), r(0:Ny), x(0:Ny) )

            IF (mode == 4) THEN
              !-----------------------------------------------------------!
              !                     Set TDMA constants
              !-----------------------------------------------------------!
              DO j = 1,Ny-1
                Rt(j) = k(j)**2 /(nu * dis(j))
                L     = Cp*sqrt(k(j)**3/dis(j)**2 + Ce**2*(nu**3/dis(j))**(1./4))
                C1    = Rt(j)**1.5/(A0*L)**2

                a(j)  =   1
                b(j)  = -(2 + C1*dy**2)
                c(j)  =   1
                r(j)  = - C1*dy**2
              END DO
              x(0:Ny) = fw(0:Ny)

              !-----------------------------------------------------------!
              !                     Boundary conditions
              !-----------------------------------------------------------!
              b(0)  = 1
              b(Ny) = 1
              a(Ny) = 0
              c(0)  = 0
              r(0)  = 0
              r(Ny) = 0

              Rt(0) = 0
              Rt(Ny)= 0

              CALL TDMA_Solver(a,b,c,r,x,Ny)
              fw(0:Ny) = x(0:Ny)
            END IF

            !-----------------------------------------------------------!
            !                  Damping funnction setting
            !-----------------------------------------------------------!
            DO j = 0,Ny
               IF (j<Ny/2) THEN
                  Y_tmp(j) = Y(j)
               ELSE
                  Y_tmp(j) = 2*del - Y(j)
               END IF
            END DO

             SELECT CASE (mode)
                CASE(0) ! No wall model
                  fm(0:Ny) = 1

                CASE(1) ! van Driest (1954)
                  fm(0:Ny) = ( 1 - exp(-Y_tmp(0:Ny)/(nu/u_tau)/A1))**2

                CASE(2) ! Launder and Sharma (1974)
                  Rt(0:Ny) = k(0:Ny)**2 /(nu * dis(0:Ny))
                  fm(0:Ny) = exp( -3.4/(1 + Rt(0:Ny)/50)**2 )

                CASE(4) ! Park et al (1997)
                  DO j = 0,Ny
                    fm1   = (1 + Cd*exp(-(Rt(j)/120)**2)*Rt(j)**(-3/4))*fw(j)**2
                    fm2   = 7.0*(4.5 + 0.3*prod(j)/dis(j))                      &
                                 /(4.5 + 1.3*prod(j)/dis(j))**2

                    fm(j) = fm1 * fm2
                    !print*,fm(j),fm1,fm2
                  END DO
             END SELECT

            DEALLOCATE(Y_tmp)

        END SUBROUTINE GETFM
