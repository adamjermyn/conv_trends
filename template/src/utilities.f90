module utilities

use star_lib
use star_def
use const_def
use math_lib

implicit none

contains

    subroutine unity(s,k,val)
        type (star_info), pointer :: s
        integer, intent(in) :: k
        real(dp), intent(out) :: val
        val = 1d0
    end subroutine unity

    subroutine max_val(s,nzlo,nzhi,do1,maxi)
        type (star_info), pointer :: s
        integer, intent(in) :: nzlo, nzhi
        interface
        subroutine do1(s,k,val)
           use star_def
           type (star_info), pointer :: s
           integer, intent(in) :: k
           real(dp), intent(out) :: val
        end subroutine do1
        end interface
        real(dp), intent(out) :: maxi

        integer :: k
        real(dp) :: val

        maxi = -1d99
        do k=nzlo,nzhi
            call do1(s,k,val)
            if (val > maxi) maxi = val
        end do

    end subroutine max_val   

    subroutine integrate_dr(s,nzlo,nzhi,do1,integral)
        type (star_info), pointer :: s
        integer, intent(in) :: nzlo, nzhi
        interface
        subroutine do1(s,k,val)
           use star_def
           type (star_info), pointer :: s
           integer, intent(in) :: k
           real(dp), intent(out) :: val
        end subroutine do1
        end interface
        real(dp), intent(out) :: avg

        integer :: k
        real(dp) :: dr, val

        integral = 0d0
        do k=nzlo,nzhi
            dr = s%dm(k) / (4d0 * pi * pow2(s%rmid(k)) * s%rho(k))
            call do1(s,k,val)
            integral = integral + dr * val
        end do

    end subroutine integrate_dr   

    subroutine integrate_dm(s,nzlo,nzhi,do1,integral)
        type (star_info), pointer :: s
        integer, intent(in) :: nzlo, nzhi
        interface
        subroutine do1(s,k,val)
           use star_def
           type (star_info), pointer :: s
           integer, intent(in) :: k
           real(dp), intent(out) :: val
        end subroutine do1
        end interface
        real(dp), intent(out) :: avg

        integer :: k
        real(dp) :: val

        integral = 0d0
        do k=nzlo,nzhi
            call do1(s,k,val)
            integral = integral + s%dm(k) * val
        end do

    end subroutine integrate_dm

    subroutine r_average(s,nzlo,nzhi,do1,avg)
        type (star_info), pointer :: s
        integer, intent(in) :: nzlo, nzhi
        interface
        subroutine do1(s,k,val)
           use star_def
           type (star_info), pointer :: s
           integer, intent(in) :: k
           real(dp), intent(out) :: val
        end subroutine do1
        end interface
        real(dp), intent(out) :: avg

        integer :: k
        real(dp) :: dr, val, integral_dr

        integral_dr = integrate_dr(s,nzlo,nzhi,do1,val)
        dr = integrate_dr(s,nzlo,nzhi,unity,val)
        avg = integral_dr / dz
    end subroutine r_average

        

end module utilities
