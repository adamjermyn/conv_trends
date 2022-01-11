module diffusivities

use star_lib
use star_def
use const_def
use math_lib

implicit none

contains

    ! Using the Spitzer 1962 formula, with the Coulomb logarithm from Wendell+1987
    subroutine viscosity(s, k, nu)
        type (star_info), pointer :: s
        integer, intent(in) :: k
        real(dp), intent(out) :: nu
        real(dp) :: lnLambda

        if (s%T(k) > 4.2d5) then
           lnLambda = -17.9d0 + 1.5d0 * log(s%T(k)) - 0.5d0 * log(s%rho(k))
        else
           lnLambda = -11.5d0 + log(s%T(k)) - 0.5d0 * log(s%rho(k))
        end if
        nu = 2.21d-15 * pow(s%T(k),2.5d0) / s%rho(k) / lnLambda
        nu = nu + 4d0 * crad * pow4(s%T(k)) / (15d0 * clight * s%opacity(k) * pow2(s%rho(k)))

    end subroutine viscosity

    subroutine thermal_diffusivity(s, k, alpha)
        type (star_info), pointer :: s
        integer, intent(in) :: k
        real(dp), intent(out) :: alpha

        alpha = 16d0 * boltz_sigma * pow3(s%T(k)) / (3d0 * s%opacity(k) * s%Cp(k) * pow2(s%rho(k)))

    end subroutine thermal_diffusivity

    subroutine prandtl(s, k, Pr)
        type (star_info), pointer :: s
        integer, intent(in) :: k
        real(dp), intent(out) :: Pr
        real(dp) :: nu, alpha

        call viscosity(s,k,nu)
        call thermal_diffusivity(s,k,alpha)
        Pr = nu / alpha
    end subroutine prandtl

    !> Computes the electrical conductivity following
    !! S.-C. YOON Oct. 10, 2003.
    !!
    !! @param abar The mean atomic mass number.
    !! @param zbar The mean atomic charge.
    !! @param rho The density (g/cm^3).
    !! @param T The temperature (K).
    !! @param eta The magnetic diffusivity (cm^2/s)
    subroutine magnetic_diffusivity(s, k, eta)
        type (star_info), pointer :: s
        integer, intent(in) :: k
        real(dp), intent(out) :: eta

        real(dp) :: abar, zbar, rho, T
        real(dp) :: gamma, xlambda, f

        abar = s%abar(k)
        zbar = s%zbar(k)
        rho = s%rho(k)
        T = s%T(k)

        gamma = 0.2275d0*pow2(zbar) * pow(rho * 1.d-6 / abar, one_third)*1.d8/T

        if (T >= 4.2d5) then
            f = sqrt(4.2d5/T)
        else
            f = 1.d0
        end if
        xlambda = sqrt(3d0*pow3(zbar))*pow(gamma,-1.5d0)*f + 1d0
        eta = 3.d11*zbar*log(xlambda)*pow(t,-1.5d0)             ! magnetic diffusivity
        eta = eta/(1.d0-1.20487d0*exp(-1.0576d0*pow(zbar,0.347044d0))) ! correction: gammae
    end subroutine magnetic_diffusivity



end module diffusivities
