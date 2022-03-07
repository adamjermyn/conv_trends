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
        real(dp), intent(out) :: integral

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
        real(dp), intent(out) :: integral

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

        call integrate_dr(s,nzlo,nzhi,do1,integral_dr)
        call integrate_dr(s,nzlo,nzhi,unity,dr)
        avg = integral_dr / dr
    end subroutine r_average

    subroutine r_average_hp_frac(s,k_top,k_bottom,h_frac,in_CZ,at_top,do1,avg)
        type (star_info), pointer :: s
        integer, intent(in) :: k_top,k_bottom
        real(dp), intent(in) :: h_frac
        logical, intent(in) :: in_CZ ! true means average on the CZ side of the boundary, false means RZ side
        logical, intent(in) :: at_top ! true means average at the outer boundary, false means inner boundary
        interface
        subroutine do1(s,k,val)
           use star_def
           type (star_info), pointer :: s
           integer, intent(in) :: k
           real(dp), intent(out) :: val
        end subroutine do1
        end interface
        real(dp), intent(out) :: avg
        real(dp) :: dr, dz, val
        integer :: k

        avg = 0d0        
        dz = 0d0
        if (at_top) then
            if (in_CZ) then
                do k=k_top,k_bottom
                    dr = s%dm(k) / (4d0 * pi * pow2(s%r(k)) * s%rho(k))

                    call do1(s,k,val)
                    avg = avg + dr * val
                    dz = dz + dr

                    if (dz > h_frac * s%scale_height(k_top) .or. k == k_bottom .or. s%brunt_N2(k) >= 0d0) then
                        avg = avg / dz
                        exit
                    end if
                end do
            else
                do k=k_top-1,1,-1
                    dr = s%dm(k) / (4d0 * pi * pow2(s%r(k)) * s%rho(k))

                    call do1(s,k,val)
                    avg = avg + dr * val
                    dz = dz + dr

                    if (dz > h_frac * s%scale_height(k_top) .or. k == 1 .or. s%brunt_N2(k) < 0d0) then
                        avg = avg / dz
                        exit
                    end if
                end do
            end if
        else
            if (in_CZ) then
                do k=k_bottom,k_top,-1
                    dr = s%dm(k) / (4d0 * pi * pow2(s%r(k)) * s%rho(k))

                    call do1(s,k,val)
                    avg = avg + dr * val
                    dz = dz + dr

                    if (dz > h_frac * s%scale_height(k_bottom) .or. k == k_top .or. s%brunt_N2(k) >= 0d0) then
                        avg = avg / dz
                        exit
                    end if
                end do
            else
                do k=k_bottom+1,s%nz
                    dr = s%dm(k) / (4d0 * pi * pow2(s%r(k)) * s%rho(k))

                    call do1(s,k,val)
                    avg = avg + dr * val
                    dz = dz + dr

                    if (dz > h_frac * s%scale_height(k_bottom) .or. k == s%nz .or. s%brunt_N2(k) < 0d0) then
                        avg = avg / dz
                        exit
                    end if
                end do
            end if
        end if
    end subroutine r_average_hp_frac

end module utilities
