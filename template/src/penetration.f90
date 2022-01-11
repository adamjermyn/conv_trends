module penetration

use star_lib
use star_def
use const_def
use math_lib

implicit none
contains




      subroutine compute_dm_core(s, id, m_core, dm_core, dr_core, dr_core_div_h, r_core, rho_core_top)
         use eos_def
         use star_lib
         type (star_info), pointer :: s
         integer, intent(in) :: id
         real(dp), parameter :: f = 0.86d0
         real(dp), parameter :: xi = 0.6d0
         integer :: k, j, nz, ierr
         real(dp) :: Lint, delta_r, V_CZ, Favg, RHS, dr, h
         real(dp) :: Rho, T, logRho, logT, Pr
         real(dp), intent(out) :: dm_core, dr_core, dr_core_div_h, r_core, m_core, rho_core_top
         real(dp), dimension(num_eos_basic_results) :: res, dres_dlnRho, dres_dlnT
         real(dp) :: dres_dxa(num_eos_d_dxa_results,s% species)
         real(dp) :: kap, dlnkap_dlnRho, dlnkap_dlnT, frac_Type2
         real(dp) :: gradr(s%nz), grada(s%nz)

         nz = s%nz

         ! Recalculate gradR and gradA assuming the composition of the core.
         do j=1,nz
            ! Call the EOS with the composition of the convective core
            Rho = s%rho(j)
            T = s%T(j)
            logRho = log10(Rho)
            logT = log10(T)
            ierr = 0
            call star_get_eos( &
               id, 0, s%xa(:,nz), & ! k = 0 means not being called for a particular cell
               Rho, logRho, T, logT, & 
               res, dres_dlnRho, dres_dlnT, &
               dres_dxa, ierr)
            grada(j) = res(i_grad_ad)

            ! Call the opacity with the composition of the convective core.
            ierr = 0
            call star_get_kap( &
               id, 0, s%zbar(nz), s%xa(:,nz), logRho, logT, &
               res(i_lnfree_e), dres_dlnRho(i_lnfree_e), dres_dlnT(i_lnfree_e), &
               kap, dlnkap_dlnRho, dlnkap_dlnT, frac_Type2, ierr)

            Pr = one_third*crad*T*T*T*T
            gradr(j) = s%P(j)*kap*s%L(j) / (16*pi*clight*s%m(j)*s%cgrav(j)*Pr)

         end do

         delta_r = 0d0
         V_CZ = 0d0
         Lint = 0d0

         ! Integrate over CZ
         do j=nz,1,-1
            if (gradr(j) < grada(j)) then
               ! Means we've hit a radiative zone
               m_core = s%m(j)
               r_core = s%r(j)
               rho_core_top = s%rho(j)
               h = s%scale_height(j)
               k = j
               exit
            end if

            dr = s%dm(j) / (4d0 * pi * pow2(s%r(j)) * s%rho(j))
            Lint = Lint + s%L_conv(j) * dr
            delta_r = delta_r + dr
            V_CZ = V_CZ + s%dm(j)/s%rho(j)

         end do 

         ! Calculate target RHS
         Favg = Lint / V_CZ
         RHS = (1d0 - f) * V_CZ * Favg
         Lint = 0d0

         ! Integrate over RZ until we find the edge of the PZ
         delta_r = 0d0
         do j=min(nz,k+1),1,-1
            dr = s%dm(j) / (4d0 * pi * pow2(s%r(j)) * s%rho(j))
            delta_r = delta_r + dr
            Lint = Lint + (xi * f * 4d0 * pi * pow2(s%r(j)) * Favg + s%L(j) * (grada(j) / gradr(j) - 1d0)) * dr

            if (Lint > RHS) then
               dm_core = s%m(j) - m_core
               dr_core = delta_r
               dr_core_div_h = delta_r / h
               exit
            end if

         end do             
      
      end subroutine compute_dm_core




end module penetration