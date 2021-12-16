module identification

use star_lib
use star_def
use const_def
use math_lib

implicit none
contains


        subroutine get_conv_regions(s, max_num_conv_regions, num_conv_regions)
           ! use mlt_def, only: convective_mixing
           type (star_info), pointer :: s
           integer, intent(in) :: max_num_conv_regions
           integer, intent(out) :: num_conv_regions
           integer :: prev_type, cur_type, cur_top, n, k
           include 'formats'

           cur_type = s% mixing_type(1)
           cur_top = 1
           n = 0

           ! Find all convective regions in the outer layers down to T_limit
           do k = 2, s%nz
              prev_type = cur_type
              cur_type = s% mixing_type(k)
              if (cur_type == prev_type .and. k < s%nz) cycle
              ! change of type from k-1 to k
              if (prev_type == convective_mixing) then
                 n = n + 1
                 s% mixing_region_type(n) = prev_type
                 s% mixing_region_top(n) = cur_top
                 if (k == s%nz) then
                    s% mixing_region_bottom(n) = k
                 else
                    s% mixing_region_bottom(n) = k-1
                 end if
                 if (n == max_num_conv_regions) exit
              end if
              cur_top = k
           end do

           num_conv_regions = n

        end subroutine get_conv_regions


  	  subroutine classify_conv_region(s, sc_top, sc_bottom, sc_type)
           type (star_info), pointer :: s
           integer, intent(in) :: sc_top, sc_bottom
           character (len=100), intent(out) ::  sc_type

           real(dp), DIMENSION(2) :: T_HI, T_HeI, T_HeII, T_FeCZ
           integer :: n, k
           include 'formats'

           ! Pass upper and lower gridpoints of convective regions, check temperature and classify

           T_HI    = (/ 3000,11000 /)     ! Rough T range for H Conv Region
           T_HeI   = (/ 11000,35000 /)    ! Rough T range for HeI Conv Region
           T_HeII  = (/ 35000,100000 /)   ! Rough T range for HeII Conv Region
           T_FeCZ  = (/ 100000,500000 /)  ! Rough T range for FeCZ Conv Region

           !write(*,*)   T_HI(1), T_HI(2), MAXVAL(s% T(sc_top:sc_bottom))

           ! Find gridpoint corresponding to max temperature (select only outer layers)

           sc_type = 'UNKNOWN'
           if ( sc_top > 0 ) then
             if (s% T(sc_top) < T_HI(2)) then
               sc_type = 'HI'
             else if (s%m(sc_bottom)/s%m(1) < 3d-2 .or. sc_bottom > s%nz - 3) then
               sc_type = 'Core'
             else if (s%T(sc_bottom) > T_FeCZ(2)) then ! It's not a core, but it's also not one of the four subsurface zones
               sc_type = 'Envelope'
             else
               do k = sc_top, sc_bottom
                  if  (s% T(k) > T_HeI(1) .AND. s% T(k) < T_HeI(2)) then
                    sc_type = 'HeI'
              	  else if (s% T(k) > T_HeII(1) .AND. s% T(k) < T_HeII(2)) then
                    sc_type = 'HeII'
              	  else if (s% T(k) > T_FeCZ(1) .AND. s% T(k) < T_FeCZ(2)) then
                    sc_type = 'FeCZ'
              	  end if
           	   end do
             end if
           end if
        end subroutine classify_conv_region


end module identification