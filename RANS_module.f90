!-----------------------------------------------------------------------------------!
!
!   PROGRAM : RANS_module.f90
!
!   PURPOSE : Moudles for RANS with k-e model
!
!                                                                2016.11.17 K.Noh
!
!-----------------------------------------------------------------------------------!

        MODULE RANS_module

            INTEGER :: Ny
            REAL(KIND=8) :: del, dy, Re_tau, nu
            CHARACTER(LEN=65) :: file_name, path_name

            REAL(KIND=8), DIMENSION(:), ALLOCATABLE   :: U, U_exac, Y


        END MODULE
