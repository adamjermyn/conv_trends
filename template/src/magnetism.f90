module magnetism

use star_lib
use star_def
use const_def
use math_lib
use diffusivities

implicit none
contains

    subroutine compute_B_diffusion_time(s, k_hi, tau)
       type (star_info), pointer :: s
       integer, intent(in) :: k_hi 
       real(dp), intent(out) :: tau
       
       integer :: k
       real(dp) :: dr, D

       tau = 0d0
       do k=1,k_hi

         ! Spitzer resistivity per https://en.wikipedia.org/wiki/Spitzer_resistivity
         call magnetic_diffusivity(s, k, D)

         ! Integrate d(r^2)/D through the RZ.
         dr = s%dm(k) / (4d0 * pi * pow2(s%r(k)) * s%rho(k))

         tau = tau + 2d0 * dr * s%r(k) / D
       end do

    end subroutine compute_B_diffusion_time

    subroutine B_equipartition_conv(s, k, B)
        type (star_info), pointer :: s
        integer, intent(in) :: k
        real(dp), intent(out) :: B

        B = s%conv_vel(k)/sqrt(4d0 * pi * s%rho(k))

    end subroutine B_equipartition_conv

    subroutine B_equipartition_pressure(s, k, B)
        type (star_info), pointer :: s
        integer, intent(in) :: k
        real(dp), intent(out) :: B

        B = sqrt(8d0*pi*s%Peos(k))

    end subroutine B_equipartition_pressure

    subroutine B_shutoff_conv(s, k, B)
        type (star_info), pointer :: s
        integer, intent(in) :: k
        real(dp), intent(out) :: B

        real(dp) :: delta_grad
        
        delta_grad = s%gradr(k) - s%grada(k)
        B = 4d0 * pi * s%rho(k) * pow2(s%csound(k)) * delta_grad
        B = max(0d0, B)
        B = sqrt(B)
    end subroutine B_shutoff_conv

    subroutine B_shutoff_down_to_tau(s,tau,B)
        type (star_info), pointer :: s
        real(dp), intent(in) :: tau
        real(dp), intent(out) :: B

        integer :: k
        real(dp) :: B_shut

        B = 0d0
        do k=1,s%nz
            call B_shutoff_conv(s,k,B_shut)
            B = max(B, B_shut)
            if (s%tau(k) > tau) return
        end do

    end subroutine B_shutoff_down_to_tau

end module magnetism