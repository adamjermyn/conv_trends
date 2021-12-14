module mdot

use star_lib
use star_def
use const_def
use math_lib

implicit none
contains

    subroutine eval_Vink_wind(w, T1, M1, L1, Z)
        real(dp), intent(in) :: T1, M1, L1, Z
        real(dp), intent(inout) :: w
        real(dp) :: alfa, w1, w2, Teff_jump, logMdot, dT, vinf_div_vesc
        real(dp), parameter :: Zsolar = 0.019d0 ! for Vink et al formula

        ! alfa = 1 for hot side, = 0 for cool side
        if (T1 > 27500d0) then
           alfa = 1
        else if (T1 < 22500d0) then
           alfa = 0
        else ! use Vink et al 2001, eqns 14 and 15 to set "jump" temperature
           Teff_jump = 1d3*(61.2d0 + 2.59d0*(-13.636d0 + 0.889d0*log10(Z/Zsolar)))
           dT = 100d0
           if (T1 > Teff_jump + dT) then
              alfa = 1
           else if (T1 < Teff_jump - dT) then
              alfa = 0
           else
              alfa = (T1 - (Teff_jump - dT)) / (2*dT)
           end if
        end if

        if (alfa > 0) then ! eval hot side wind (eqn 24)
           vinf_div_vesc = 2.6d0 ! this is the hot side galactic value
           vinf_div_vesc = vinf_div_vesc*pow(Z/Zsolar,0.13d0) ! corrected for Z
           logMdot = &
              - 6.697d0 &
              + 2.194d0*log10(L1/Lsun/1d5) &
              - 1.313d0*log10(M1/Msun/30) &
              - 1.226d0*log10(vinf_div_vesc/2d0) &
              + 0.933d0*log10(T1/4d4) &
              - 10.92d0*pow2(log10(T1/4d4)) &
              + 0.85d0*log10(Z/Zsolar)
           w1 = exp10(logMdot)
        else
           w1 = 0
        end if

        if (alfa < 1) then ! eval cool side wind (eqn 25)
           vinf_div_vesc = 1.3d0 ! this is the cool side galactic value
           vinf_div_vesc = vinf_div_vesc*pow(Z/Zsolar,0.13d0) ! corrected for Z
           logMdot = &
              - 6.688d0 &
              + 2.210d0*log10(L1/Lsun/1d5) &
              - 1.339d0*log10(M1/Msun/30) &
              - 1.601d0*log10(vinf_div_vesc/2d0) &
              + 1.07d0*log10(T1/2d4) &
              + 0.85d0*log10(Z/Zsolar)
           w2 = exp10(logMdot)
        else
           w2 = 0
        end if

        w = alfa*w1 + (1 - alfa)*w2

    end subroutine eval_Vink_wind



end module mdot