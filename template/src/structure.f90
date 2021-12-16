module structure

use star_lib
use star_def
use const_def
use math_lib

implicit none
contains

    subroutine scale_height(s,k,h)
        type (star_info), pointer :: s
        integer, intent(in) :: k
        real(dp), intent(out) :: h
        h = s%scale_height(k)       
    end subroutine scale_height

    subroutine density(s,k,rho)
        type (star_info), pointer :: s
        integer, intent(in) :: k
        real(dp), intent(out) :: rho
        rho = s%rho(k)       
    end subroutine density

    subroutine dtau_dr(s,k,dtdr)
        type (star_info), pointer :: s
        integer, intent(in) :: k
        real(dp), intent(out) :: dtdr

        dtdr = s%opacity(k) * s%rho(k)
    end subroutine dtau_dr

    subroutine conv_vel(s,k,v)
        type (star_info), pointer :: s
        integer, intent(in) :: k
        real(dp), intent(out) :: v

        v = s%conv_vel(k)
    end subroutine conv_vel

    subroutine sound_speed(s,k,v)
        type (star_info), pointer :: s
        integer, intent(in) :: k
        real(dp), intent(out) :: v

        v = s%csound(k)
    end subroutine sound_speed

    subroutine gravity(s,k,g)
        type (star_info), pointer :: s
        integer, intent(in) :: k
        real(dp), intent(out) :: g

        g = s%grav(k)
    end subroutine gravity

    subroutine gradr_sub_grada(s,k,g)
        type (star_info), pointer :: s
        integer, intent(in) :: k
        real(dp), intent(out) :: g

        g = max(0d0,s%gradr(k) - s%grada(k))
    end subroutine gradr_sub_grada
    

end module structure