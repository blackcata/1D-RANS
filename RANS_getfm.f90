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
                     dy, fw, A0, Cd

            IMPLICIT NONE
            INTEGER :: i,j
            REAL(KIND=8) :: L, fm1, fm2, C1
            REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: Rt,Y_tmp,a,b,c,r,x

            ALLOCATE( Y_tmp(0:Ny), Rt(0:Ny))
            ALLOCATE( a(0:Ny), b(0:Ny), c(0:Ny), r(0:Ny), x(0:Ny) )

            IF (mode ==4) THEN
              !-----------------------------------------------------------!
              !                     Set TDMA constants
              !-----------------------------------------------------------!
              DO j = 1,Ny-1
                Rt(j) = k(j)**2 /(nu * dis(j))
                L     = k(j)**1.5/dis(j)
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
              a(Ny) = -1
              c(0)  = -1
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

               SELECT CASE (mode)
                  CASE(0) ! No wall model
                    fm(j) = 1

                  CASE(1) ! van Driest (1954)
                    fm(j) = ( 1 - exp(-Y_tmp(j)/(nu/u_tau)/A0))**2

                  CASE(2) ! Launder and Sharma (1974)
                    Rt(j) = k(j)**2 /(nu * dis(j))
                    fm(j) = exp( -3.4/(1 + Rt(j)/50)**2 )

                  CASE(4) ! Park et al (1997)
                    fm1   = (1 + Cd*exp(-(Rt(j)/120)**2)*Rt(j)**(-3/4))*fw(j)**2
                    fm2   = 7.0*(4.5 + 0.3*prod(j)/dis(j))                      &
                                 /(4.5 + 1.3*prod(j)/dis(j))**2

                    fm(j) = fm1 * fm2
                    print*,fm(j),fm1,fm2
               END SELECT

            END DO

            DEALLOCATE(Y_tmp)

        END SUBROUTINE GETFM