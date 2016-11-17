!-----------------------------------------------------------------------------------!
!
!   PROGRAM : RANS_poiseuille.f90
!
!   PURPOSE : To make the initial poiseuille flow
!
!                                                                2016.11.18 K.Noh
!
!-----------------------------------------------------------------------------------!

      SUBROUTINE POISEUILLE

          USE RANS_module,                                                      &
              ONLY : Ny, dy, del, Re_tau, nu

          USE RANS_module,                                                      &
              ONLY : U, U_exac

          IMPLICIT NONE
          INTEGER :: i, j

      END SUBROUTINE POISEUILLE
