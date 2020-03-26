!  SchroedingerGleichung.f90 
!
!  FUNCTIONS:
!  SchroedingerGleichung - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: SchroedingerGleichung
!
!  PURPOSE:  Entry point for the console application.
!  COMPILED WITH INTEL FORTRAN COMPILER AND THE MKL TTF Libray   
!
!
!****************************************************************************
include "mkl_dfti.f90"

    
    program SchroedingerGleichung
    use MKL_DFTI
    implicit none

    
    type(DFTI_DESCRIPTOR), POINTER :: My_Desc1_Handle, My_Desc2_Handle
    Integer :: Status, i
    integer, parameter :: n = 128, t = 20
    Complex, dimension(n) :: psi, psi_t, psi_lsg
    Complex, parameter :: j = cmplx(0.0_8, 1.0_8)
    double precision, parameter :: delta_x = 5d-1, x_0 = -n/2 * delta_x, delta_k = 2.0_8  * acos(-1.0_8)/ (delta_x * n), k_a = 0d0, k_0=1.0_8 !k_a =1.0_8 - n/2* delta_k
    double precision :: x,k, a, b, c
    
    open(3, file="psi_lsg.dat", action="write")
    open(4, file="psi_t_x.dat", action="write")
    !Bestimmen der Werte
    do i = 0, n-1
        x = i*delta_x + x_0
        psi(i+1) = exp(-x**2 / 2.0_8) * exp(j * k_0 * x) 
        
        a = 1.0_8/ sqrt(sqrt(cmplx(1.0_8, t*1.0_8)))
        b = exp(-1.0_8*((x - k_0 * t) **2 /(2.0_8* cmplx(1.0_8, t))))
        c = exp(j * (k_0 * x - k_0**2 * t /2.0_8))
        psi_lsg(i+1) =  a* b * c
        write(3, *) x, real(psi_lsg(i+1))
    end do
    
    
    !Bestimme die Fouriertranformierte
    Status = DftiCreateDescriptor( My_Desc1_Handle, DFTI_SINGLE, DFTI_COMPLEX, 1, n )
    Status = DftiCommitDescriptor( My_Desc1_Handle )
    Status = DftiComputeForward( My_Desc1_Handle, psi )
    Status = DftiFreeDescriptor(My_Desc1_Handle)
    
    if (Status == 0) then
        !Zeitentwicklung durchführen-> 
        do i = 0, n-1
            k = i * delta_k  + k_a
            psi_t(i+1)= exp(-j *k**2 * t/ 2.0_8) * psi(i+1)
        end do
        
    end if
    
     !zurück transformieren
    Status = DftiCreateDescriptor( My_Desc2_Handle, DFTI_SINGLE, DFTI_COMPLEX, 1, n )
    Status = DftiCommitDescriptor( My_Desc2_Handle )
    Status = DftiComputeBackward( My_Desc2_Handle, psi_t )
    Status = DftiFreeDescriptor(My_Desc2_Handle)
    
    if (Status == 0) then    
        do i = 0, n-1
            x = i*delta_x + x_0
            write(4,*) x, real(psi_t(i+1)) /n
        end do
    end if
    
    close(3)
    close(4)
    
    end program SchroedingerGleichung

