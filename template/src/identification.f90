module identification

use star_lib
use star_def
use const_def
use math_lib

implicit none
contains


        subroutine get_conv_regions(s, max_num_conv_regions, num_conv_regions)
            type (star_info), pointer :: s
            integer, intent(in) :: max_num_conv_regions
            integer, intent(out) :: num_conv_regions
            logical :: in_CZ
            integer :: n, k
            include 'formats'

            s% mixing_region_top = -1
            s% mixing_region_bottom = -1

            in_CZ = (s%brunt_N2(1) < 0d0)
            if (in_CZ) then
                n = 1
                s% mixing_region_top(n) = 1
            else
                n = 0
            end if 

            do k=2,s%nz
!                write(*,*) k, s%m(k)/Msun, s%brunt_N2(k), in_CZ
               if (in_CZ .and. s%brunt_N2(k) >= 0d0) then
                  ! change of type from k-1 to k, no longer convective
                  in_CZ = .false.
                  s% mixing_region_bottom(n) = k-1
               else if (.not. in_CZ .and. s%brunt_N2(k) < 0d0) then
                  n = n+1
                  s% mixing_region_top(n) = k
                  in_CZ = .true.
               end if
               if (n == max_num_conv_regions+1) exit
            end do
            if (in_CZ) s%mixing_region_bottom(n) = s%nz


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
             if (s%m(sc_bottom)/s%m(1) < 3d-2 .or. sc_bottom > s%nz - 3 .and. s%m(1) > Msun) then ! We call it an envelope CZ at low masses.
               sc_type = 'Core'
             else if (s%T(sc_bottom) > T_FeCZ(2) .and. s%m(1) < 2d0*Msun) then ! It's not a core, but it's also not one of the four subsurface zones
               sc_type = 'Envelope'
             else if (s% T(sc_top) < T_HI(2)) then
               sc_type = 'HI'
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

           write(*,*) sc_top, sc_bottom, sc_type
        end subroutine classify_conv_region


end module identification