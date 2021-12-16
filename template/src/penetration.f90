module penetration

use star_lib
use star_def
use const_def
use math_lib

implicit none
contains



  subroutine compute_dm_core(s, m_core, dm_core, dr_core, dr_core_div_h)
     type (star_info), pointer :: s
     real(dp), parameter :: f = 0.86d0
     real(dp), parameter :: xi = 0.6d0
     integer :: k, j
     real(dp) :: Lint, delta_r, Lavg, RHS, dr, h
     real(dp), intent(out) :: m_core, dm_core, dr_core, dr_core_div_h

     delta_r = 0d0
     Lint = 0d0

     ! Integrate over CZ
     do j=s%nz,1,-1
        if (s%brunt_N2(j) > 0d0) then
           ! Means we've hit a radiative zone
           m_core = s%m(j)
           h = s%scale_height(j)
           k = j
           exit
        end if

        dr = s%dm(j) / (4d0 * pi * pow2(s%r(j)) * s%rho(j))
        Lint = Lint + s%L_conv(j) * dr
        delta_r = delta_r + dr

     end do 

     ! Calculate target RHS
     RHS = (1d0 - f) * Lint
     Lavg = Lint / delta_r
     Lint = 0d0

     ! Integrate over RZ until we find the edge of the PZ
     delta_r = 0d0
     do j=min(s%nz,k+1),1,-1
        dr = s%dm(j) / (4d0 * pi * pow2(s%r(j)) * s%rho(j))
        delta_r = delta_r + dr
        Lint = Lint + (xi * f * Lavg + s%L(j) * (s%grada(j) / s%gradr(j) - 1d0)) * dr

        if (Lint > RHS) then
           dm_core = s%m(j) - m_core
           dr_core = delta_r
           dr_core_div_h = delta_r / h
           exit
        end if

     end do             
  
  end subroutine compute_dm_core



end module penetration