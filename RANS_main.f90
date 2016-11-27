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
              ONLY : path_name, itmax, resi, tol, dis_new, nu, u_tau, U, Ny

          IMPLICIT NONE
          INTEGER :: it
          REAL(8) :: time_sta, time_end

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
          CALL CPU_TIME(time_sta)
          DO it = 0, itmax
            CALL GETFM
            CALL GETNUT
            CALL GETU
            CALL GETPROD
            CALL GETK
            CALL GETDIS
            WRITE(*,"(A,I8.7,A,E14.7,2X,F10.7,2X,E14.7)")'Iteration(',it,') : '            &
                   ,resi,U(Ny)/u_tau,dis_new(0)/(nu*(u_tau**2/nu)**2)
            IF (resi < tol) EXIT
          END DO
          CALL CPU_TIME(time_end)
          WRITE(*,"(A,F12.7,A)")'Total calculation time : ',time_end-time_sta,'s'

          !---------------------------------------------!
          !              Write final result
          !---------------------------------------------!
          CALL OUTPUT

        END PROGRAM RANS_main
