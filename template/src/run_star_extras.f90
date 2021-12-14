! ***********************************************************************
!
!   Copyright (C) 2010  Bill Paxton
!
!   this file is part of mesa.
!
!   mesa is free software; you can redistribute it and/or modify
!   it under the terms of the gnu general library public license as published
!   by the free software foundation; either version 2 of the license, or
!   (at your option) any later version.
!
!   mesa is distributed in the hope that it will be useful,
!   but without any warranty; without even the implied warranty of
!   merchantability or fitness for a particular purpose.  see the
!   gnu library general public license for more details.
!
!   you should have received a copy of the gnu library general public license
!   along with this software; if not, write to the free software
!   foundation, inc., 59 temple place, suite 330, boston, ma 02111-1307 usa
!
! ***********************************************************************

      module run_star_extras

      use star_lib
      use star_def
      use const_def
      use math_lib

      implicit none

      ! these routines are called by the standard run_star check_model
      contains

        subroutine extras_controls(id, ierr)
           integer, intent(in) :: id
           integer, intent(out) :: ierr
           type (star_info), pointer :: s
           ierr = 0
           call star_ptr(id, s, ierr)
           if (ierr /= 0) return

           ! this is the place to set any procedure pointers you want to change
           ! e.g., other_wind, other_mixing, other_energy  (see star_data.inc)

           ! Uncomment these lines if you wish to use the functions in this file,
           ! otherwise we use a null_ version which does nothing.
           s% extras_startup => extras_startup
           s% extras_start_step => extras_start_step
           s% extras_check_model => extras_check_model
           s% extras_finish_step => extras_finish_step
           s% extras_after_evolve => extras_after_evolve
           s% how_many_extra_history_columns => how_many_extra_history_columns
           s% data_for_extra_history_columns => data_for_extra_history_columns
           s% how_many_extra_profile_columns => how_many_extra_profile_columns
           s% data_for_extra_profile_columns => data_for_extra_profile_columns

           s% how_many_extra_history_header_items => how_many_extra_history_header_items
           s% data_for_extra_history_header_items => data_for_extra_history_header_items
           s% how_many_extra_profile_header_items => how_many_extra_profile_header_items
           s% data_for_extra_profile_header_items => data_for_extra_profile_header_items

           s% how_many_other_mesh_fcns => how_many_my_other_mesh_fcns
           s% other_mesh_fcn_data => other_mesh_fcn_data

           ! Once you have set the function pointers you want,
           ! then uncomment this (or set it in your star_job inlist)
           ! to disable the printed warning message,
            s% job% warn_run_star_extras =.false.

        end subroutine extras_controls

      subroutine how_many_my_other_mesh_fcns(id, n)
         integer, intent(in) :: id
         integer, intent(out) :: n
         n = 1
      end subroutine how_many_my_other_mesh_fcns

      subroutine other_mesh_fcn_data( &
            id, nfcns, names, gval_is_xa_function, vals1, ierr)
         integer, intent(in) :: id
         integer, intent(in) :: nfcns
         character (len=*) :: names(:)
         logical, intent(out) :: gval_is_xa_function(:) ! (nfcns)
         real(dp), pointer :: vals1(:) ! =(nz, nfcns)
         integer, intent(out) :: ierr
         integer :: nz, k
         real(dp), pointer :: vals(:,:)
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         names(1) = 'N2_function'
         gval_is_xa_function(1) = .false.
         nz = s% nz
         vals(1:nz,1:nfcns) => vals1(1:nz*nfcns)
         do k=1,nz
            vals(k,1) = 15d0 * log10(s%brunt_N2(k) + 1d-14)
         end do
      end subroutine other_mesh_fcn_data  


      subroutine extras_startup(id, restart, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
      end subroutine extras_startup


      integer function extras_start_step(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_start_step = 0
      end function extras_start_step

      ! returns either keep_going, retry, or terminate.
      integer function extras_check_model(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_check_model = keep_going         
         if (.false. .and. s% star_mass_h1 < 0.35d0) then
            ! stop when star hydrogen mass drops to specified level
            extras_check_model = terminate
            write(*, *) 'have reached desired hydrogen mass'
            return
         end if


         ! if you want to check multiple conditions, it can be useful
         ! to set a different termination code depending on which
         ! condition was triggered.  MESA provides 9 customizeable
         ! termination codes, named t_xtra1 .. t_xtra9.  You can
         ! customize the messages that will be printed upon exit by
         ! setting the corresponding termination_code_str value.
         ! termination_code_str(t_xtra1) = 'my termination condition'

         ! by default, indicate where (in the code) MESA terminated
         if (extras_check_model == terminate) s% termination_code = t_extras_check_model
      end function extras_check_model


        integer function how_many_extra_history_columns(id)
           integer, intent(in) :: id
           integer :: ierr
           type (star_info), pointer :: s
           ierr = 0
           call star_ptr(id, s, ierr)
           if (ierr /= 0) return
           how_many_extra_history_columns = 135
        end function how_many_extra_history_columns


        subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
           use diffusivities
           use utilities
           use mdot
           integer, intent(in) :: id, n

           ! Outputs
           integer, intent(out) :: ierr
           character (len=maxlen_history_column_name) :: names(n)
           real(dp) :: vals(n)

           ! Intermediates
           type (star_info), pointer :: s
           integer :: i,j
           integer ::  k, num_conv_regions
           character(len=100) :: name

           ! Quantities that are one per CZ
           integer, parameter :: nQs = 30
           integer, parameter :: nZs = 5 ! Max # of CZs
           logical :: sc_exists(nZs)
           integer, dimension(nZs) :: sc_top, sc_bottom
           character(len=100) :: Q_names(nQs)
           character(len=100) :: cz_names(nZs)
           real(dp) :: outputs(nQs, nZs)
           real(dp) :: Ra(nZs), Pm(nZs), Pr(nZs), nu(nZs), alpha(nZs), eta(nZs), dr(nZs)

           ! Quantities that are one per star
           integer, parameter :: nOffs = 10
           character(len=100) :: one_off_names(nOffs)
           real(dp) :: one_off_outputs(nOffs)
           real(dp) :: wind

           ! Set up star pointer
           ierr = 0
           call star_ptr(id, s, ierr)
           if (ierr /= 0) return

           ! Name the CZs
           i = 1
           cz_names(i) = 'HI'; i = i+1
           cz_names(i) = 'HeI'; i = i+1
           cz_names(i) = 'HeII'; i = i+1
           cz_names(i) = 'FeCZ'; i = i+1
           cz_names(i) = 'Core'; i = i+1

           ! Identify which CZs, if any, exist

           ! Quantities in the CZs
           ! Quantities with profiles are radially averaged, weighted by dr unless otherwise specified.  
           do k=1,nZs
            if (sc_exists(k)) then
               i = 0

               i = i+1
               Q_names(i) = 'viscosity'
               call r_average(s, sc_top(k), sc_bottom(k), viscosity, outputs(i,k))

               i = i+1
               Q_names(i) = 'thermal_diffusivity'
               call r_average(s, sc_top(k), sc_bottom(k), thermal_diffusivity, outputs(i,k))

               i = i+1
               Q_names(i) = 'magnetic_diffusivity'
               call r_average(s, sc_top(k), sc_bottom(k), magnetic_diffusivity, outputs(i,k))

               i = i+1
               Q_names(i) = 'cz_dm'
               call integrate_dm(s, sc_top(k), sc_bottom(k), unity, outputs(i,k))

               i = i+1
               Q_names(i) = 'cz_dr'
               call integrate_dr(s, sc_top(k), sc_bottom(k), unity, outputs(i,k))

               i = i+1
               Q_names(i) = 'cz_top_r'
               outputs(i,k) = s%r(s%sc_top(k))

               i = i+1
               Q_names(i) = 'cz_bottom_r'
               outputs(i,k) = s%r(s%sc_bottom(k))

               i = i+1
               Q_names(i) = 'B_shutoff'
               outputs(i,k) = max_val(s, sc_top(k), sc_bottom(k), B_shutoff_conv, outputs(i,k))

            end if
           end do

           ! One-off quantities
           i = 1
           one_off_outputs(i) = s% mixing_length_alpha
           one_off_names(i) = 'mixing_length_alpha'; i=i+1

           one_off_outputs(i) = s% rho(1)
           one_off_names(i) = 'rho_surf'; i=i+1

           call eval_Vink_wind(wind, s%Teff, s%m(1), s%L(1), 1d0 - s%surface_h1 - s%surface_he4)
           one_off_outputs(i) = wind
           one_off_names(i) = 'wind_mdot'; i=i+1


           ! Pack output 
           k = 1
           do i=1,nQs
            do j=1,nZs
               name = cz_names(j) // '_' // Q_names(i)
               names(k) = name
               vals(k) = outputs(i,j)
               k = k + 1
            end do
           end do

           do i=1,nOffs
            names(k) = one_off_names(k)
            vals(k) = one_off_outputs(k)
            k = k + 1
           end do

        end subroutine data_for_extra_history_columns



        integer function how_many_extra_profile_columns(id)
           use star_def, only: star_info
           integer, intent(in) :: id
           integer :: ierr
           type (star_info), pointer :: s
           ierr = 0
           call star_ptr(id, s, ierr)
           if (ierr /= 0) return
           how_many_extra_profile_columns = 3
        end function how_many_extra_profile_columns

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


      !> Compute the magnetic diffusivity from the electric conductivity.
      !! @param sig The electrical conductivity (1/s).
      !! @param eta The magnetic diffusivity (output, cm^2/s).
      real(dp) function calc_eta(sig) result(eta)
         real(dp), intent(in) :: sig

      end function calc_eta


      subroutine compute_diffusivities(id, ierr, sc_top, sc_bottom, Ra, Pm, Pr, nu_avg, alpha_avg, eta_avg, dz)
         use const_def

         type (star_info), pointer :: s
         integer, intent(in) :: id, sc_top, sc_bottom
         integer, intent(out) :: ierr
         real(dp), intent(out) :: Ra, Pm, Pr, alpha_avg, nu_avg, eta_avg, dz

         integer :: k
         real(dp) :: alpha, lnLambda, nu
         real(dp) :: dr
         real(dp) :: gradr_sub_grada_avg, g_avg, eta, sig

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         dz = s%r(sc_top) - s%r(sc_bottom)

         eta_avg = 0d0
         g_avg = 0d0
         gradr_sub_grada_avg = 0d0
         do k=sc_top,sc_bottom

            sig = calc_sige(s%abar(k),s%zbar(k),s%rho(k),s%T(k))
            eta = calc_eta(sig)

            dr = s%dm(k) / (4d0 * pi * pow2(s%rmid(k)) * s%rho(k))

            gradr_sub_grada_avg = gradr_sub_grada_avg + dr * max(0d0,s%gradr(k) - s%grada(k))
            g_avg = g_avg + dr * s%grav(k)
            eta_avg = eta_avg + dz * eta
         end do
         gradr_sub_grada_avg = gradr_sub_grada_avg / dz
         g_avg = g_avg / dz
         eta_avg = eta_avg / dz

         Pr = nu_avg / alpha_avg
         Pm = nu_avg / eta_avg
         Ra = pow3(dz) * g_avg * gradr_sub_grada_avg / (nu_avg * alpha_avg)

      end subroutine compute_diffusivities



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

  	  subroutine get_conv_radii(id, ierr, sc_top, sc_bottom, r_top, r_bottom)
           type (star_info), pointer :: s
           integer, intent(in) :: id
           integer, intent(out) :: ierr
           real(dp) :: r_top, r_bottom
           integer :: n, k, sc_top, sc_bottom
           include 'formats'

           ierr = 0
           call star_ptr(id, s, ierr)
           if (ierr /= 0) return
           r_top = 0d0
           r_bottom = 0d0

           if ( sc_top > 0 ) then
           	r_top = s% r(sc_top)
           	r_bottom = s% r(sc_bottom)
           end if
        end subroutine get_conv_radii


        subroutine get_hp_radii(id, ierr, hp, k_hp)
           type (star_info), pointer :: s
           integer, intent(in) :: id
           integer, intent(out) :: ierr
           integer :: k_hp
           real(dp) :: hp
           include 'formats'

           ierr = 0
           call star_ptr(id, s, ierr)
           if (ierr /= 0) return
           k_hp = 1
  !        write(*,*) s% lnP
  !         do while (((s% lnP(k_hp) - s% lnP(1)) < LOG10(hp)) .and. (k_hp < s% nz))
           do while ((LOG(EXP(s% lnP(k_hp))/EXP(s% lnP(1))) < hp).and.(k_hp < s% nz))
           !           write(*,*)  s% lnP(k_hp),  s% lnP(1), LOG10(hp)
           	k_hp = k_hp + 1
           end do
        end subroutine get_hp_radii

        subroutine get_max_fc(id, ierr, fcmax, sc_top, sc_bottom)
           type (star_info), pointer :: s
           integer, intent(in) :: id
           integer, intent(out) :: ierr
           real(dp) :: fcmax
           integer :: n, k, sc_top, sc_bottom
           include 'formats'

           ierr = 0
           call star_ptr(id, s, ierr)
           if (ierr /= 0) return
           fcmax = 0
           ! fcmax = MAXVAL((s% conv_vel(sc_top:sc_bottom)**3.0) * s% rho(sc_top:sc_bottom))
           fcmax = MAXVAL(s% L_conv(sc_top:sc_bottom)/s% L(sc_top:sc_bottom))
        end subroutine get_max_fc

        subroutine get_conv_mass(id, ierr, sc_top, sc_bottom, total_mass)
        type (star_info), pointer :: s
           integer, intent(in) :: id
           integer, intent(out) :: ierr
           real(dp) :: total_mass
           integer :: n, k, sc_top, sc_bottom
           include 'formats'

           ierr = 0
           call star_ptr(id, s, ierr)
           if (ierr /= 0) return

           total_mass = 0d0
           if ( sc_top > 0 ) then
           	total_mass = SUM(s% dm(sc_top:sc_bottom))
           end if

        end subroutine  get_conv_mass

        subroutine get_conv_velocities(id, ierr, v_max, v_aver, sc_top, sc_bottom,b_eq,b_max,rho_aver)
           type (star_info), pointer :: s
           integer, intent(in) :: id
           integer, intent(out) :: ierr
           real(dp) :: v_max, v_aver, b_eq, b_max, rho_aver, rho_v_max
           integer :: n, k, sc_top, sc_bottom, i_v_max
           include 'formats'

           ierr = 0
           call star_ptr(id, s, ierr)
           if (ierr /= 0) return

           v_max = 0d0
           v_aver = 0d0
           rho_aver = 0d0
           b_eq = 0d0
           b_max = 0d0
           rho_aver = 0d0
           i_v_max = 0
           rho_v_max = 0d0
           ! Calculate Max and average V_conv in Conv region (average in dr, see Eq. 6 in Cantiello et al. 2009)


           if ( sc_top > 0 ) then
           	do k = sc_top, sc_bottom
              	v_aver = v_aver + (s% r(k) - s% r(k+1)) * s% conv_vel(k)
              	rho_aver = rho_aver + (s% r(k) - s% r(k+1)) * s% rho(k)
              !  write(*,*) 'DRHO: ', (s% r(k) - s% r(k+1)) * s% rho(k), s% rho(k)
            !    write(*,*) 'DV: ',(s% r(k) - s% r(k+1))*s% conv_vel(k),s% conv_vel(k)
           	end do
           	v_max = MAXVAL(s% conv_vel(sc_top:sc_bottom))
            i_v_max = MAXLOC (s% conv_vel(sc_top:sc_bottom), DIM=1)
            rho_v_max = s% rho(i_v_max)
           	v_aver = v_aver /( s% r(sc_top) - s% r(sc_bottom) )
           	rho_aver = rho_aver /( s% r(sc_top) - s% r(sc_bottom) )


           end if
           ! Calculate B_equipartition and B_max

           b_eq = (v_aver)*(4.0*pi*rho_aver)**(0.5)
           b_max = (v_max)*(4.0*pi*rho_v_max)**(0.5)
           !b_max = (v_max)*(4.0*pi*rho_aver)**(0.5) ! For convective core this would work better

           !write(*,*) v_aver, v_max, rho_aver, b_eq, b_max

        end subroutine get_conv_velocities

        subroutine get_pressure_eq_field(id, ierr, sc_top, sc_bottom,b_p_eq,b_p_max)
           type (star_info), pointer :: s
           integer, intent(in) :: id
           integer, intent(out) :: ierr
           real(dp) :: p_max, p_aver, b_p_eq, b_p_max
           integer :: n, k, sc_top, sc_bottom
           include 'formats'

           ierr = 0
           call star_ptr(id, s, ierr)
           if (ierr /= 0) return


           b_p_eq = 0d0
           b_p_max = 0d0
           p_aver = 0d0

           if ( sc_top > 0 ) then
           	do k = sc_top, sc_bottom
              	p_aver = p_aver + (s% r(k) - s% r(k+1)) * EXP(s% lnP(k))
           	end do
           	p_max = EXP(MAXVAL(s% lnP(sc_top:sc_bottom)))
           	p_aver = p_aver /( s% r(sc_top) - s% r(sc_bottom) )
           end if
           ! Calculate B_Pressure_equipartition and B_max

           b_p_eq = (p_aver*8.0*pi)**(0.5)
           b_p_max = (p_max*8*pi)**(0.5)

  !         write(*,*) b_p_eq, b_p_max

        end subroutine get_pressure_eq_field

        subroutine get_average_hp(id, ierr, sc_top, sc_bottom, hp_aver)
           type (star_info), pointer :: s
           integer, intent(in) :: id
           integer, intent(out) :: ierr
           real(dp) :: hp_aver
           integer :: n, k, sc_top, sc_bottom
           include 'formats'

           ierr = 0
           call star_ptr(id, s, ierr)
           if (ierr /= 0) return

           hp_aver = 0d0
           ! Calculate average hp in convective region (dr weighted)

           if ( sc_top > 0 ) then
           	do k = sc_top, sc_bottom
              	hp_aver = hp_aver + (s% r(k) - s% r(k+1)) *s% scale_height(k)
           	end do
              hp_aver = hp_aver /( s% r(sc_top) - s% r(sc_bottom) )
           end if

        end subroutine get_average_hp

        subroutine compute_optical_depths(s, sc_top, sc_bottom, tau_cz, tau_surf)
           type (star_info), pointer :: s
           integer, intent(in) :: sc_top, sc_bottom
           real(dp), intent(out) :: tau_cz, tau_surf
           integer :: k
           real(dp) :: dr

           tau_cz = 0d0
           tau_surf = 0d0

           do k=sc_top,sc_bottom
             dr = s%dm(k) / (4d0 * pi * pow2(s%r(k)) * s%rho(k))
             tau_cz = tau_cz + dr * s%rho(k) * s%opacity(k)
           end do

           do k=1,sc_top
             dr = s%dm(k) / (4d0 * pi * pow2(s%r(k)) * s%rho(k))
             tau_surf = tau_surf + dr * s%rho(k) * s%opacity(k)
           end do

        end subroutine compute_optical_depths


        real(dp) function compute_B_diffusion_time(s, k_hi) result(tau)
           type (star_info), pointer :: s
           integer, intent(in) :: k_hi 
           
           integer :: k
           real(dp) :: eta, ne, lambdaD, dr, D

           tau = 0d0
           do k=1,k_hi
             ne = exp(s%lnfree_e(k)) * avo * (s%rho(k) / s%abar(k))
             lambdaD = sqrt(kerg * s%T(k) / (4d0 * pi * pow2(qe) * (s%Z2bar(k))))

             ! Spitzer resistivity per https://en.wikipedia.org/wiki/Spitzer_resistivity
             eta = 4d0 * sqrt(2d0 * pi) / 3d0
             eta = eta * s%zbar(k) * pow2(qe) * sqrt(me) / pow(kerg * s%T(k), 1.5d0) ! Transverse Spitzer resistivity
             eta = eta * log(4d0 * pi * pow3(lambdaD) * ne / 3d0)

             ! eta is currently a resistivity, with CGS units of erg*cm / (erg^(3/2) / g^(1/2)) = erg*cm/(erg*cm/s) = s
             ! But we want a diffusivity. The MHD diffusivity is
             D = eta * pow2(clight) / (4d0 * pi) ! This has units cm^2/s.

             ! Integrate d(r^2)/D through the RZ.
             dr = s%dm(k) / (4d0 * pi * pow2(s%r(k)) * s%rho(k))

             tau = tau + 2d0 * dr * s%r(k) / D
           end do

        end function compute_B_diffusion_time

        real(dp) function compute_buoyant_time(s, k_hi, B_dynamo) result(tau)
           type (star_info), pointer :: s
           integer, intent(in) :: k_hi 
           real(dp), intent(in) :: B_dynamo
           
           integer :: k
           real(dp) :: dr, Pmag, alpha, beta, v_therm, v_drag, v

           Pmag = (pow2(B_dynamo) / (8d0 * pi))
           tau = 0d0
           do k=1,k_hi
            alpha = 16d0 * boltz_sigma * pow3(s%T(k)) / (3d0 * s%opacity(k) * pow2(s%rho(k)) * s%cp(k))
            beta = s%P(k) / Pmag
            v_therm = 2d0 * alpha / (s%scale_height(k) * beta * abs(s%grada(k) - s%gradr(k)) * abs(4d0 - 3d0 / s%chiRho(k)))
            v_drag = sqrt(2d0 * Pmag / s%rho(k))

            v = min(v_drag, v_therm)

            ! Integrate d(r^2)/D through the RZ.
            dr = s%dm(k) / (4d0 * pi * pow2(s%r(k)) * s%rho(k))
            tau = tau + dr / v
           end do

        end function compute_buoyant_time


        subroutine get_microturb(mach1_aver_ahp, rho1_aver, rho_surf,v1_aver_ahp, v1_surf_aver)
           real(dp) :: v1_surf_aver, v1_aver_ahp, mach1_aver_ahp, rho1_aver, rho_surf

           v1_surf_aver = 0d0
           if (rho_surf /= 0) then
  	         v1_surf_aver = v1_aver_ahp * (mach1_aver_ahp * rho1_aver/rho_surf)**(1./2.)
           end if

        end subroutine get_microturb

  	  subroutine get_turnover(mixing_length_alpha,v_HI_aver, HI_hp_aver, turnover_HI)
           real(dp) :: v_HI_aver, HI_hp_aver, turnover_HI, mixing_length_alpha


            turnover_HI = 0d0
            if (v_HI_aver /= 0) then
              turnover_HI = mixing_length_alpha*HI_hp_aver/v_HI_aver
           endif


        end subroutine get_turnover

        subroutine get_bsurf(rho_surf, rho_HI_aver, b_HI_aver, b_HI_max, b_HI_surf, b_HI_surf_max)
           real(dp) :: rho_surf, rho_HI_aver, b_HI_aver, b_HI_max, b_HI_surf, b_HI_surf_max
           b_HI_surf = 0d0
           b_HI_surf_max = 0d0
            if (rho_HI_aver /= 0) then
              b_HI_surf = b_HI_aver * (rho_surf/rho_HI_aver)**(2./3.)
              b_HI_surf_max = b_HI_max * (rho_surf/rho_HI_aver)**(2./3.)
           endif


        end subroutine get_bsurf



  	  subroutine get_conv_ahp(id, ierr, sc_top, sc_bottom, v_aver_ahp, mach_top_cz, mach_aver_ahp, rho_aver_ahp)
           type (star_info), pointer :: s
           integer, intent(in) :: id
           integer, intent(out) :: ierr
           real(dp) ::  v_aver_ahp, mach_aver_ahp, rho_aver_ahp, cs_aver_ahp, delta_r, mach_top_cz
           integer :: n, k, sc_top, sc_bottom, kk
           include 'formats'

           ierr = 0
           call star_ptr(id, s, ierr)
           if (ierr /= 0) return

           v_aver_ahp = 0d0
           rho_aver_ahp = 0d0
           cs_aver_ahp = 0d0
           mach_aver_ahp = 0d0
           mach_top_cz = 0d0
           kk = 0

           ! Calculate rho_aver and v_aver in top alpha*Hp of convective zone ! Follows Eq.6 in Cantiello et al. 2009

           if ( sc_top > 0 ) then
           	do k = sc_top, sc_bottom
              	if (s% r(k) > s% r(sc_top) - s% mixing_length_alpha * s% scale_height(sc_top)) then
                 		rho_aver_ahp = rho_aver_ahp + (s% r(k) - s% r(k+1)) * s% rho(k)
                 		v_aver_ahp =  v_aver_ahp + (s% r(k) - s% r(k+1)) * s%  conv_vel(k)
                 		cs_aver_ahp = cs_aver_ahp + (s% r(k) - s% r(k+1)) * s% csound(k)
                 		kk = k
              	end if
           	end do
           end if
           rho_aver_ahp = rho_aver_ahp / ( s% r(sc_top) - s% r(kk) )
           v_aver_ahp = v_aver_ahp / ( s% r(sc_top) - s% r(kk) )
           cs_aver_ahp = cs_aver_ahp / ( s% r(sc_top) - s% r(kk) )

           if (cs_aver_ahp /=0) then
           	mach_aver_ahp = v_aver_ahp/cs_aver_ahp
           end if

           if (s% csound(sc_top) /=0) then
           	mach_top_cz = s%  conv_vel(sc_top) / s% csound(sc_top)
           end if


        end subroutine get_conv_ahp


        subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
           use star_def, only: star_info, maxlen_profile_column_name
           use const_def, only: dp
           use math_lib, only: safe_log10
           integer, intent(in) :: id, n, nz
           character (len=maxlen_profile_column_name) :: names(n)
           real(dp) :: vals(nz,n), rsun , pi, clight
           integer, intent(out) :: ierr
           type (star_info), pointer :: s
           integer :: k
           ierr = 0
           call star_ptr(id, s, ierr)
           if (ierr /= 0) return

           pi = 3.1415
           rsun = 6.96e10
           clight = 2.99792458e10

           !!!

           if (n /= 3) stop 'data_for_extra_profile_columns'
           names(1) = 'log_f_div_rhocs3'

           do k = 1, nz
              vals(k,1) = (s% L(k)/(4*pi*(s% r(k)**2)))/(s% Rho(k) * s% csound(k)**3)
              vals(k,1) = safe_log10(vals(k,1))
           end do

           names(2) = 'log_fconv_div_rhocs3'

           do k = 1, nz
              vals(k,2) = (s% conv_vel(k) / s% csound(k))**3
              vals(k,2) = safe_log10(vals(k,2))
           end do

           names(3) = 'log_Fedd_div_rhocs3'

           do k = 1, nz
             vals(k,3) = 4*pi*clight*s% cgrav(k)*s% m_grav(k)/s% opacity(k) ! L_edd
             vals(k,3) = vals(k,3) / (4*pi*(s% r(k)**2))                    ! F_edd
             vals(k,3) = vals(k,3) /(s% Rho(k) * s% csound(k)**3)           ! F_edd/F_max
             vals(k,3) = safe_log10(vals(k,3))
           end do

           !!!!


        end subroutine data_for_extra_profile_columns



      integer function how_many_extra_history_header_items(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_history_header_items = 0
      end function how_many_extra_history_header_items


      subroutine data_for_extra_history_header_items(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         type(star_info), pointer :: s
         integer, intent(out) :: ierr
         ierr = 0
         call star_ptr(id,s,ierr)
         if(ierr/=0) return

         ! here is an example for adding an extra history header item
         ! also set how_many_extra_history_header_items
         ! names(1) = 'mixing_length_alpha'
         ! vals(1) = s% mixing_length_alpha

      end subroutine data_for_extra_history_header_items


      integer function how_many_extra_profile_header_items(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_header_items = 0
      end function how_many_extra_profile_header_items


      subroutine data_for_extra_profile_header_items(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(n)
         type(star_info), pointer :: s
         integer, intent(out) :: ierr
         ierr = 0
         call star_ptr(id,s,ierr)
         if(ierr/=0) return

         ! here is an example for adding an extra profile header item
         ! also set how_many_extra_profile_header_items
         ! names(1) = 'mixing_length_alpha'
         ! vals(1) = s% mixing_length_alpha

      end subroutine data_for_extra_profile_header_items


      ! returns either keep_going or terminate.
      ! note: cannot request retry; extras_check_model can do that.
      integer function extras_finish_step(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_finish_step = keep_going

         ! to save a profile, 
            ! s% need_to_save_profiles_now = .true.
         ! to update the star log,
            ! s% need_to_update_history_now = .true.

         ! see extras_check_model for information about custom termination codes
         ! by default, indicate where (in the code) MESA terminated
         if (extras_finish_step == terminate) s% termination_code = t_extras_finish_step
      end function extras_finish_step
      
      
      subroutine extras_after_evolve(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
      end subroutine extras_after_evolve


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

      end module run_star_extras
