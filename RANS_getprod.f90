!-----------------------------------------------------------------------------------!
!
!   PROGRAM : RANS_getprod.f90
!
!   PURPOSE : Get P(TKE Production) using k-e model relation
!
!                                                             2016.11.20 K.Noh
!
!-----------------------------------------------------------------------------------!

        SUBROUTINE GETPROD

            USE RANS_module,                                                    &
              ONLY : Ny, dy, u_tau, nu

            USE RANS_module,                                                    &
              ONLY : U_new, nu_T, prod

            IMPLICIT NONE
            INTEGER :: i,j

            DO j = 1,Ny-1
              prod(j) = nu_T(j) * ((U_new(j+1) - U_new(j-1))/(2*dy))**2.
            END DO

            !-----------------------------------------------------------!
            !                     Boundary conditions
            !-----------------------------------------------------------!
            prod(0)  = nu_T(0)  * (u_tau**2/nu)**2.
            prod(Ny) = nu_T(Ny) * (u_tau**2/nu)**2.

        END SUBROUTINE GETPROD
