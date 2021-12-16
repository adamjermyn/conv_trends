module identification

use star_lib
use star_def
use const_def
use math_lib

implicit none
contains


        subroutine get_conv_regions_above_T(id, T_limit, ierr, num_conv_regions)
           ! use mlt_def, only: convective_mixing
           type (star_info), pointer :: s
           integer, intent(in) :: id
           integer, intent(out) :: ierr
           real(dp) :: T_limit
           integer :: prev_type, cur_type, cur_top, n, k, num_conv_regions, max_num_conv_regions, n_limit
           include 'formats'

           ierr = 0
           call star_ptr(id, s, ierr)
           if (ierr /= 0) return

           ierr = 0
           cur_type = s% mixing_type(1)
           cur_top = 1
           n = 0
           n_limit = 0
           max_num_conv_regions = 4 ! Max number of convective regions allowed


           ! Find gridpoint corresponding to max temperature (select only outer layers)
           do k = 1, s% nz
              if (s% T(k) < T_limit) then
                    n_limit = k
              end if
           end do

           ! Find all convective regions in the outer layers down to T_limit
           do k = 2, n_limit
              prev_type = cur_type
              cur_type = s% mixing_type(k)
              if (cur_type == prev_type .and. k < n_limit) cycle
              ! change of type from k-1 to k
              if (prev_type == convective_mixing) then
                 n = n + 1
                 s% mixing_region_type(n) = prev_type
                 s% mixing_region_top(n) = cur_top
                 if (k == n_limit) then
                    s% mixing_region_bottom(n) = k
                 else
                    s% mixing_region_bottom(n) = k-1
                 end if
                 if (n == max_num_conv_regions) exit
              end if
              cur_top = k
           end do

           num_conv_regions = n

        end subroutine get_conv_regions_above_T



      subroutine get_convective_core(id, sc_convective_core,ierr)
           type (star_info), pointer :: s
           integer, intent(in) :: id
           integer, intent(out) :: ierr, sc_convective_core
           integer :: cur_type, k
           include 'formats'

           ierr = 0
           call star_ptr(id, s, ierr)
           if (ierr /= 0) return

           k = s% nz
           cur_type = s% mixing_type(k)
           write(*,*) 'Convective type', cur_type
           do while (cur_type == 1 .and. k > 2)
             cur_type = s% mixing_type(k)
             k = k - 1
           end do
           sc_convective_core = k
      end subroutine get_convective_core


  	  subroutine classify_conv_region_above_T(id, ierr, sc_top, sc_bottom, sc_type)

           type (star_info), pointer :: s
           integer, intent(in) :: id
           integer, intent(out) :: ierr

           character (len=7) ::  sc_type
           real(dp), DIMENSION(2) :: T_HI, T_HeI, T_HeII, T_FeCZ
           integer :: n, k, sc_top, sc_bottom
           include 'formats'

           ! Pass upper and lower gridpoints of convective regions, check temperature and classify

           ierr = 0
           call star_ptr(id, s, ierr)
           if (ierr /= 0) return


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
             else
               do k = sc_top, sc_bottom
                  if  (s% T(k) > T_HeI(1) .AND. s% T(k) < T_HeI(2)) then
                    sc_type = 'HeI'
              	else if (s% T(k) > T_HeII(1) .AND. s% T(k) < T_HeII(2)) then
                    sc_type = 'HeII'
              	else if (s% T(k) > T_FeCZ(1) .AND. s% T(k) < T_FeCZ(2)) then
                    sc_type = 'FeCZ'
              	else
                    sc_type = 'UNKNOWN'
              	end if
           	 end do
             end if
           end if
           write(*,*) 'Type: ', s% T(sc_top), s% T(sc_bottom), sc_type
        end subroutine classify_conv_region_above_T


end module identification