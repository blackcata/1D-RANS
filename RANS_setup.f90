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
                ONLY : Ny, del, dy, Re_tau, nu,                                 &
                       Cm, Ce1, Ce2, Sk, Se, alpha, beta

            USE RANS_module,                                                    &
                ONLY : U, U_exac, U_new, Y, k, dis, nu_T

            IMPLICIT NONE
            INTEGER :: i,j

            !--------------------------------!
            !     Constants for simulation
            !--------------------------------!
            Ny  = 100         ! the number of grid cells
            del = 2          ! the channel-half height
            dy  = (2*del)/NY  ! grid size

            Re_tau = 180      ! Reynolds number based on the friction velocity
            nu     = 3.5000e-4 ! Kinematic viscosity of reference data

            !------------------------------!
            !    Constants for k-e model
            !------------------------------!
            Cm  = 0.09
            Ce1 = 1.44
            Ce2 = 1.92
            Sk  = 1.0
            Se  = 1.3

            !------------------------------!
            !      Relaxation factors
            !------------------------------!
            alpha = 0.7
            beta  = 0.7

            ALLOCATE( U(0:NY), U_new(0:Ny), U_exac(0:Ny), Y(0:Ny) )
            ALLOCATE( k(0:Ny), dis(0:Ny), nu_T(0:Ny) )

            !--------------------------------!
            !       Initial Conditions
            !--------------------------------!
            DO j = 0,Ny
              Y(j)      = j*dy
              U(j)      = 0
              k(j)      = 0.0100
              dis(j)    = 0.0015
              nu_T(j)   = 0
              U_new(j)  = 0
              U_exac(j) = -(nu/(2*del))*(Re_tau/del)**2 * Y(j) * (Y(j)-2*del)
            END DO

        END SUBROUTINE SETUP
