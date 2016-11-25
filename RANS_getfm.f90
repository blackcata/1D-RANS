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
              ONLY : Ny, Cm, k, dis, nu_T, fm, nu, u_tau, Y, del, mode

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

               Ry      = sqrt(k(j))*Y_tmp(j)
               Rt      = k(j)**2 /(nu * dis(j))

               SELECT CASE (mode)
                  CASE(1) ! van Driest (1954)
                    fm(j)   = ( 1 - exp(-Y_tmp(j)/(nu/u_tau)/40))**2

                  CASE(2) ! Launder and Sharma (1974)
                    fm(j)   = exp( -3.4/(1+Rt/50)**2 )

                  CASE(4) ! Park et al (1997)

               END SELECT


            END DO

            DEALLOCATE(Y_tmp)

        END SUBROUTINE GETFM
