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
              ONLY : Ny, dy

            USE RANS_module,                                                    &
              ONLY : U_new, nu_T, prod

            IMPLICIT NONE
            INTEGER :: i,j

            DO j = 1,Ny-1
              prod(j) = nu_T(j) * ((U_new(j+1) - U_new(j-1))/(2*dy))**2
            END DO

        END SUBROUTINE GETPROD
