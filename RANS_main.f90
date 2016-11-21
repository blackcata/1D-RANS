!-----------------------------------------------------------------------------------!
!
!   PROGRAM : RANS_main.f90
!
!   PURPOSE : To solve the 1D fully developed channel flow
!           with RANS model and k-epsilon model.
!
!                                                             2016.11.17 K.Noh
!
!-----------------------------------------------------------------------------------!


        PROGRAM RANS_main

          USE RANS_module,                                                      &
              ONLY : path_name, itmax, resi, tol

          IMPLICIT NONE
          INTEGER :: it

          !---------------------------------------------!
          !             Make Result folder
          !---------------------------------------------!
          path_name = 'RESULT'
          CALL SYSTEM('mkdir '//TRIM(path_name))
          CALL SYSTEM('rm -rf ./'//TRIM(path_name)//'/*')

          !---------------------------------------------!
          !               Initial Setting
          !---------------------------------------------!
          CALL SETUP
          CALL POISEUILLE

          !---------------------------------------------!
          !                  Main loop
          !---------------------------------------------!
          DO it = 0, itmax
            CALL GETNUT
            CALL GETU
            CALL GETPROD
            CALL GETK
            CALL GETDIS
            IF (resi < tol) EXIT
          END DO

          !---------------------------------------------!
          !              Write final result
          !---------------------------------------------!
          CALL OUTPUT

        END PROGRAM RANS_main
