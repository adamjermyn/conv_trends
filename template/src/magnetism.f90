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
        use const_def
        use eos_def
        use eos_lib

        type (star_info), pointer :: s
        integer, intent(in) :: k
        real(dp), intent(out) :: B

        integer :: ierr, id
        real(dp) :: delta_grad
        real(dp) :: Q
        real(dp) :: Gamma_00, Gamma_p1, Gamma_m1
        real(dp) :: dGamma
        
        real(dp) :: d_dxa_eos(num_eos_basic_results,s%species)
        real(dp) :: eos_results(num_eos_basic_results), d_dlnRho_const_T(num_eos_basic_results),d_dlnT_const_Rho(num_eos_basic_results)

        id = s%id

        call star_get_eos( &
           id, k, s%xa(:,k), & 
           s%Rho(k), arg_not_provided, s%T(k), arg_not_provided, & 
           eos_results, d_dlnRho_const_T, d_dlnT_const_Rho, d_dxa_eos, ierr)

        Gamma_00 = eos_results(i_gamma1)

        if (k < s%nz) then
            call star_get_eos( &
               id, k+1, s%xa(:,k+1), &
               s%Rho(k+1), arg_not_provided, s%T(k+1), arg_not_provided, & 
               eos_results, d_dlnRho_const_T, d_dlnT_const_Rho, d_dxa_eos, ierr)

            Gamma_p1 = eos_results(i_gamma1)
            dGamma = 2d0 * (Gamma_p1 - Gamma_00) / (Gamma_p1 + Gamma_00)
            dGamma = dGamma * 2d0 * (s%Peos(k+1) - s%Peos(k)) / (s%Peos(k+1) + s%Peos(k)) 
        else
            call star_get_eos( &
               id, k-1, s%xa(:,k-1), &
               s%Rho(k-1), arg_not_provided, s%T(k-1), arg_not_provided, & 
               eos_results, d_dlnRho_const_T, d_dlnT_const_Rho, d_dxa_eos, ierr)

            Gamma_m1 = eos_results(i_gamma1)
            dGamma = 2d0 * (Gamma_00 - Gamma_m1) / (Gamma_p1 + Gamma_00) 
            dGamma = dGamma * 2d0 * (s%Peos(k) - s%Peos(k-1)) / (s%Peos(k) + s%Peos(k-1))               
        end if

        delta_grad = s%gradr(k) - s%grada(k)
        Q = 1 + 4d0 * s%Prad(k) / (s%Peos(k) - s%Prad(k))
        B = 4 * pi * s%rho(k) * pow2(s%csound(k))
        B = B * Q * delta_grad / (1 - Q * delta_grad + dGamma)

        if (B < 0d0) then ! Approximate formula for handling extreme cases
            B = P * delta_grad
        end if
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