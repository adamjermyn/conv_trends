module magnetism

use star_lib
use star_def
use const_def
use math_lib

implicit none
contains


    subroutine get_B_shutoff_conv(s, k, B)
        use const_def
        use eos_def
        use eos_lib

        type (star_info), pointer :: s
        integer, intent(in) :: id
        integer, intent(out) :: ierr
        real(dp), intent(out) :: B

        integer :: n, k, sc_top, sc_bottom
        real(dp) :: B
        real(dp) :: delta_grad
        real(dp) :: Q
        real(dp) :: Gamma_00, Gamma_p1
        real(dp) :: dGamma

        real(dp) :: eos_results(num_eos_basic_results), d_dlnRho_const_T(num_eos_basic_results),d_dlnT_const_Rho(num_eos_basic_results),d_dabar_const_TRho(num_eos_basic_results), d_dzbar_const_TRho(num_eos_basic_results)

        ierr = 0
        call star_ptr(id, s, ierr)
        if (ierr /= 0) return

        call eosDT_get( &
            s%eos_handle, s%Z(k), s%X(k), s%abar(k), s%zbar(k),  &
            s%species, s%chem_id, s%net_iso, s%xa(:,k), &
            s%Rho(k), arg_not_provided, s%T(k), arg_not_provided,  &
            eos_results, d_dlnRho_const_T, d_dlnT_const_Rho, &
            d_dabar_const_TRho, d_dzbar_const_TRho, ierr)

        Gamma = eos_results(i_gamma1)

        if (k < s%nz) then
            call eosDT_get( &
                s%eos_handle, s%Z(k+1), s%X(k+1), s%abar(k+1), s%zbar(k+1),  &
                s%species, s%chem_id, s%net_iso, s%xa(:,k+1), &
                s%Rho(k+1), arg_not_provided, s%T(k+1), arg_not_provided,  &
                eos_results, d_dlnRho_const_T, d_dlnT_const_Rho, &
                d_dabar_const_TRho, d_dzbar_const_TRho, ierr)

            Gamma_p1 = eos_results(i_gamma1)
            dGamma = 2d0 * (Gamma_p1 - Gamma_00) / (Gamma_p1 + Gamma_00)
            dGamma = dGamma * 2d0 * (s%p(k+1) - s%p(k)) / (s%p(k+1) + s%p(k)) 
        else
            call eosDT_get( &
                s%eos_handle, s%Z(k-1), s%X(k-1), s%abar(k-1), s%zbar(k-1),  &
                s%species, s%chem_id, s%net_iso, s%xa(:,k-1), &
                s%Rho(k-1), arg_not_provided, s%T(k-1), arg_not_provided,  &
                eos_results, d_dlnRho_const_T, d_dlnT_const_Rho, &
                d_dabar_const_TRho, d_dzbar_const_TRho, ierr)

            Gamma_m1 = eos_results(i_gamma1)
            dGamma = 2d0 * (Gamma_00 - Gamma_m1) / (Gamma_p1 + Gamma_00) 
            dGamma = dGamma * 2d0 * (s%p(k) - s%p(k-1)) / (s%p(k) + s%p(k-1))               
        end if

        delta_grad = s%gradr(k) - s%grada(k)
        Q = 1 + 4d0 * s%Prad(k) / (s%P(k) - s%Prad(k))
        B = 4 * pi * s%rho(k) * pow2(s%csound(k))
        B = B * Q * delta_grad / (1 - Q * delta_grad + dGamma)

        if (B < 0d0) B = 0d0
        B = sqrt(B)

    end subroutine get_B_shutoff_conv_region_above_T


end module magnetism