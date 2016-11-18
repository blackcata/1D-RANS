!-----------------------------------------------------------------------------------!
!
!   PROGRAM : RANS_setup.f90
!
!   PURPOSE : Setup the 1D turbulent channel flow for RANS model
!
!                                                                 2016.11.17 K.Noh
!
!-----------------------------------------------------------------------------------!

        SUBROUTINE SETUP

            USE RANS_module,                                                    &
                ONLY : Ny, del, dy, Re_tau, nu, Cm, Ce1, Ce2, Sk, Se

            USE RANS_module,                                                    &
                ONLY : U, U_exac, Y, k, dis, nu_T

            IMPLICIT NONE
            INTEGER :: i,j

            Ny  = 100         ! the number of grid cells
            del = 1           ! the channel-half height
            dy  = (2*del)/NY  ! grid size

            Re_tau = 180      ! Reynolds number based on the friction velocity
            nu     = 1.004e-6 ! Kinematic viscosity of water at 20C

            Cm  = 0.09
            Ce1 = 1.44
            Ce2 = 1.92
            Sk  = 1.0
            Se  = 1.3

            ALLOCATE( U(0:NY), U_exac(0:Ny), Y(0:Ny) )
            ALLOCATE( k(0:Ny), dis(0:Ny), nu_T(0:Ny) )

            DO j = 0,Ny
              Y(j) = -del + j*dy
              U_exac(j) = -nu/2*(Re_tau/del)**2 * (Y(j)**2 - del**2)
            END DO

        END SUBROUTINE SETUP
