!-----------------------------------------------------------------------------------!
!
!  PROGRAM: TDMA_Solver.f90
!
!  PURPOSE: To solve the inverse problem at Tridiagonal matrix.
!           This is called the TDMA Solver.
!
!  p.s if you run TDMA in the loop, you have to initialize a,b,c arrays
!
!                                                                2016.03.03 K.Noh
!
!   log 
!   2016.03.03 First make a TDMA Solver and add some precautions
!              This subroutine has a point to be careful that  after subroutine 
!              a,b,c,d values are changed dramatically so if you use more than once
!              this subroutine you have to reset a,b,c,d value every time
!                                     
!-----------------------------------------------------------------------------------!

        SUBROUTINE TDMA_Solver(a,b,c,r,x,n)

            IMPLICIT NONE
            INTEGER :: it, n
            REAL(8) :: a(0:n), b(0:n), c(0:n), r(0:n), x(0:n)

            a(0) = 0.0
            c(n) = 0.0
    
            !---Forward Sweeping-----

            DO it = 1,n
                b(it) = b(it-1) * b(it) - a(it) * c(it-1)
                c(it) = b(it-1) * c(it)
                r(it) = b(it-1) * r(it) - a(it) * r(it-1)
            END DO

            !---Backward Sweeping----

            x(n) = r(n) / b(n)

            DO it = n-1,0,-1
                x(it) = (r(it) - c(it) * x(it+1))/b(it)
            END DO

        END SUBROUTINE TDMA_Solver

