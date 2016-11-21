!-----------------------------------------------------------------------------------!
!
!   PROGRAM : RANS_output.f90
!
!   PURPOSE : Make each variables 1D plot in the RESULT folder
!
!                                                             2016.11.18 K.Noh
!
!-----------------------------------------------------------------------------------!

    SUBROUTINE OUTPUT

        USE RANS_module,                                                        &
            ONLY : file_name, path_name

        USE RANS_module,                                                        &
            ONLY : Ny,del,Re_tau,u_tau,nu,band

        USE RANS_module,                                                        &
            ONLY : Y, U, U_new, U_exac, k, k_new, dis, dis_new, nu_T

        IMPLICIT NONE
        INTEGER :: i, j

        !----------------------------------------------------------!
        !              Outputs for U(Mean velocity)
        !----------------------------------------------------------!
        file_name = TRIM(path_name)//'/U.plt'
        OPEN(100,FILE=file_name,FORM='FORMATTED',POSITION='APPEND')
        WRITE(100,*) 'VARIABLES="Y","Y+","U","U+","U_exac"'
        DO j = 0+band,Ny-band
            WRITE(100,"(5F15.9)") Y(j),Y(j)/(del/Re_tau),                       &
                                  U_new(j), U_new(j)/u_tau, U_exac(j)
        END DO
        CLOSE(100)

        !----------------------------------------------------------!
        !         Outputs for k(Turbulent kinetic energy)
        !----------------------------------------------------------!
        file_name = TRIM(path_name)//'/k.plt'
        OPEN(100,FILE=file_name,FORM='FORMATTED',POSITION='APPEND')
        WRITE(100,*) 'VARIABLES="Y","Y+","k","k+"'
        DO j = 0+band,Ny-band
            WRITE(100,"(4F15.9)") Y(j), Y(j)/(del/Re_tau),                      &
                                  k_new(j), k_new(j)/(u_tau**2)
        END DO
        CLOSE(100)

        !----------------------------------------------------------!
        !            Outputs for epsilon(dissipation)
        !----------------------------------------------------------!
        file_name = TRIM(path_name)//'/dissipation.plt'
        OPEN(100,FILE=file_name,FORM='FORMATTED',POSITION='APPEND')
        WRITE(100,*) 'VARIABLES="Y","Y+","dis","dis+"'
        DO j = 0+band,Ny-band
            WRITE(100,"(4F15.9)") Y(j), Y(j)/(del/Re_tau),                      &
                                  dis_new(j), dis_new(j)/(nu*(u_tau**2/nu)**2)
        END DO
        CLOSE(100)

        DEALLOCATE(U,U_new,U_exac,Y,k,k_new,dis,dis_new,nu_T)

    END SUBROUTINE OUTPUT
