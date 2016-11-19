!-----------------------------------------------------------------------------------!
!
!   PROGRAM : RANS_getnut.f90
!
!   PURPOSE : Get Nu_t(turbulent kinematic viscosity) using k-e model relation
!
!                                                             2016.11.19 K.Noh
!
!-----------------------------------------------------------------------------------!

        SUBROUTINE GETNUT

            USE RANS_module,                                                    &
              ONLY : Ny, dy, nu, del, Re_tau

            USE RANS_module,                                                    &
              ONLY : U, nu_T

            IMPLICIT NONE
            INTEGER :: i,j
            REAL(KIND=8) :: C0, u_tau

            u_tau = Re_tau*nu / del
            C0 = u_tau**2 / del

        END SUBROUTINE GETNUT
