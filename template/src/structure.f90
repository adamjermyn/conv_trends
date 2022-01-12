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

    subroutine Prad_div_P(s,k,beta)
        type (star_info), pointer :: s
        integer, intent(in) :: k
        real(dp), intent(out) :: beta
        beta = s%Prad(k)/s%P(k)       
    end subroutine Prad_div_P

    subroutine dtau_dr(s,k,dtdr)
        type (star_info), pointer :: s
        integer, intent(in) :: k
        real(dp), intent(out) :: dtdr

        dtdr = s%opacity(k) * s%rho(k)
    end subroutine dtau_dr

    subroutine bottom_r(s,k,r)
        type (star_info), pointer :: s
        integer, intent(in) :: k
        real(dp), intent(out) :: r

        if (k == s%nz) then
            r = 0d0
        else
            r = s%r(k+1)
        end if
    end subroutine bottom_r

    subroutine bottom_tau(s,k,tau)
        type (star_info), pointer :: s
        integer, intent(in) :: k
        real(dp), intent(out) :: tau

        if (k == s%nz) then
            tau = s%tau(k)+s%opacity(k)*s%rho(k)*s%r(k)
        else
            tau = s%tau(k+1)
        end if
    end subroutine bottom_tau

    subroutine Nusselt(s,k,Nu)
        type (star_info), pointer :: s
        integer, intent(in) :: k
        real(dp), intent(out) :: Nu

        ! Nu = 1 + F_c / (F_r - F_{r,ad})
        ! = (F_r + F_c - F_{r,ad}) / (F_r - F_{r,ad})
        ! = (F - F_{r,ad}) / (F_r - F_{r,ad})
        !
        ! Now we can write
        !
        ! F_r = F (grad/gradR)
        ! F_{r,ad} = F (gradA/gradR)
        !
        ! So
        !
        ! Nu = (1 - gradA/gradR)/(grad/gradR - gradA/gradR)
        ! = (gradR - gradA) / (grad - gradA) 

        Nu = (s%gradR(k) - s%gradA(k)) / s%gradT_sub_grada(k)
    end subroutine Nusselt

    subroutine conv_vel(s,k,v)
        type (star_info), pointer :: s
        integer, intent(in) :: k
        real(dp), intent(out) :: v

        v = s%conv_vel(k)
    end subroutine conv_vel

    subroutine conv_vel_div_r(s,k,v)
        type (star_info), pointer :: s
        integer, intent(in) :: k
        real(dp), intent(out) :: v

        v = s%conv_vel(k)/s%r(k)
    end subroutine conv_vel_div_r

    subroutine L_conv_div_L(s,k,val)
        type (star_info), pointer :: s
        integer, intent(in) :: k
        real(dp), intent(out) :: val

        val = s%L_conv(k) / s%L(k)
    end subroutine L_conv_div_L

    subroutine L_div_Ledd(s,k,val)
        type (star_info), pointer :: s
        integer, intent(in) :: k
        real(dp), intent(out) :: val

        val = s%L(k) / (4d0 * pi * standard_cgrav * s%m(k) * clight / s%opacity(k))
    end subroutine L_div_Ledd


    subroutine stiffness_top(s,k_top,k_bottom,val)
        type (star_info), pointer :: s
        integer, intent(in) :: k_top, k_bottom
        real(dp), intent(out) :: val
        real(dp) :: brunt2_CZ, brunt2_RZ, dz, dr, vc
        integer :: k

        if (k_top == 1) then
            val = 0d0
            return
        end if

        dz = 0d0
        vc = 0d0
        do k=k_top,k_bottom
            dr = s%dm(k) / (4d0 * pi * pow2(s%r(k)) * s%rho(k))
            vc = vc + s%conv_vel(k) * dr
            dz = dz + dr
            if (dz > 0.2d0*s%scale_height(k_top) .or. k == k_bottom) then
                vc = vc / dz
                brunt2_CZ = pow2(vc/s%scale_height(k_top))
                exit
            end if
        end do

        dz = 0d0
        do k=k_top,1,-1
            dr = s%dm(k) / (4d0 * pi * pow2(s%r(k)) * s%rho(k))
            dz = dz + dr
            if (dz > 0.2d0*s%scale_height(k_top) .or. k == 1) then
                brunt2_RZ = s%brunt_N2(k)
                exit
            end if
        end do

        val = brunt2_RZ / brunt2_CZ

    end subroutine stiffness_top

    subroutine stiffness_bottom(s,k_top,k_bottom,val)
        type (star_info), pointer :: s
        integer, intent(in) :: k_top, k_bottom
        real(dp), intent(out) :: val
        real(dp) :: brunt2_CZ, brunt2_RZ, dz, dr, vc
        integer :: k

        if (k_bottom == s%nz) then
            val = 0d0
            return
        end if

        dz = 0d0
        vc = 0d0
        do k=k_bottom,k_top,-1
            dr = s%dm(k) / (4d0 * pi * pow2(s%r(k)) * s%rho(k))
            vc = vc + s%conv_vel(k) * dr
            dz = dz + dr
            if (dz > 0.2d0*s%scale_height(k_top) .or. k == k_top) then
                vc = vc / dz
                brunt2_CZ = pow2(vc/s%scale_height(k_top))
                exit
            end if
        end do

        dz = 0d0
        do k=k_bottom,s%nz-1
            dr = s%dm(k) / (4d0 * pi * pow2(s%r(k)) * s%rho(k))
            dz = dz + dr
            if (dz > 0.2d0*s%scale_height(k_top) .or. k == s%nz-1) then
                brunt2_RZ = s%brunt_N2(k)
                exit
            end if
        end do

        val = brunt2_RZ / brunt2_CZ

    end subroutine stiffness_bottom

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