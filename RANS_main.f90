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
              ONLY : path_name

          IMPLICIT NONE

          path_name = 'RESULT'
          CALL SYSTEM('mkdir '//TRIM(path_name))
          CALL SYSTEM('rm -rf ./'//TRIM(path_name)//'/*')

          CALL SETUP
          CALL POISEUILLE
          CALL OUTPUT

        END PROGRAM RANS_main
