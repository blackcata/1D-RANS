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

        USE RANS_module,                                                         &
            ONLY : Ny,del,Re_tau

        USE RANS_module,                                                         &
            ONLY : Y, U, U_new, U_exac, k, k_new

        IMPLICIT NONE
        INTEGER :: i, j

        file_name = TRIM(path_name)//'/U.plt'
        OPEN(100,FILE=file_name,FORM='FORMATTED',POSITION='APPEND')
        WRITE(100,*) 'VARIABLES="Y","Y+",U","U_new","U_exac"'
        DO j = 0,Ny
            WRITE(100,"(5F15.9)") Y(j),Y(j)/(del/Re_tau), U(j), U_new(j), U_exac(j)
        END DO
        CLOSE(100)

        file_name = TRIM(path_name)//'/k.plt'
        OPEN(100,FILE=file_name,FORM='FORMATTED',POSITION='APPEND')
        WRITE(100,*) 'VARIABLES="Y","Y+","k","k_new"'
        DO j = 0,Ny
            WRITE(100,*) Y(j), Y(j)/(del/Re_tau), k(j), k_new(j)
        END DO
        CLOSE(100)

    END SUBROUTINE OUTPUT
