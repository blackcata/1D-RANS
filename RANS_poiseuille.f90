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
          REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: a,b,c,r,x

          ALLOCATE ( a(0:Ny), b(0:Ny), c(0:Ny), r(0:Ny), x(0:Ny))

          a(0:Ny) = 1
          b(0:Ny) = -2
          c(0:Ny) = 1
          x(0:Ny) = 0
          r(1:Ny-1) = -nu*(dy*Re_tau/del)**2

          b(0) = 1
          b(Ny) = 1
          r(0) = 0
          r(Ny) = 0

          CALL TDMA_Solver(a,b,c,r,x,Ny)

          U(0:Ny) = x(0:Ny)
          
          DO j = 0,Ny
            print*,j,x(j)
          END DO

      END SUBROUTINE POISEUILLE
