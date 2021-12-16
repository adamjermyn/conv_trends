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
           use identification
           use diffusivities
           use utilities
           use mdot
           use structure
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
               Q_names(i) = 'pressure_scale_height'
               call r_average(s, sc_top(k), sc_bottom(k), scale_height, outputs(i,k))

               i = i+1
               Q_names(i) = 'pressure_scale_height_top'
               outputs(i,k) = s%scale_height(s%sc_top(k))

               i = i+1
               Q_names(i) = 'pressure_scale_height_bottom'
               outputs(i,k) = s%scale_height(s%sc_bottom(k))

               i = i+1
               Q_names(i) = 'density'
               call r_average(s, sc_top(k), sc_bottom(k), density, outputs(i,k))

               i = i+1
               Q_names(i) = 'sound_speed'
               call r_average(s, sc_top(k), sc_bottom(k), sound_speed, outputs(i,k))

               i = i+1
               Q_names(i) = 'gravity'
               call r_average(s, sc_top(k), sc_bottom(k), gravity, outputs(i,k))

               i = i+1
               Q_names(i) = 'gradr_sub_grada'
               call r_average(s, sc_top(k), sc_bottom(k), gradr_sub_grada, outputs(i,k))

               i = i+1
               Q_names(i) = 'conv_vel'
               call r_average(s, sc_top(k), sc_bottom(k), conv_vel, outputs(i,k))

               i = i+1
               Q_names(i) = 'conv_vel_max'
               call max_val(s, sc_top(k), sc_bottom(k), conv_vel, outputs(i,k))

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
               Q_names(i) = 'cz_dtau'
               outputs(i,k) = integrate_dr(s, sc_top(k), sc_bottom(k), dtau_dr, outputs(i,k))

               i = i+1
               Q_names(i) = 'cz_top_xm'
               outputs(i,k) = integrate_dm(s, 1, sc_top(k), unity, outputs(i,k))

               i = i+1
               Q_names(i) = 'cz_bottom_xm'
               outputs(i,k) = integrate_dm(s, 1, sc_bottom(k), unity, outputs(i,k))

               i = i+1
               Q_names(i) = 'cz_top_tau'
               outputs(i,k) = s%tau(sc_top(k))

               i = i+1
               Q_names(i) = 'cz_bottom_tau'
               outputs(i,k) = s%tau(sc_bottom(k))

               i = i+1
               Q_names(i) = 'B_shutoff'
               outputs(i,k) = max_val(s, sc_top(k), sc_bottom(k), B_shutoff_conv, outputs(i,k))

               i = i+1
               Q_names(i) = 'B_equipartition_conv_max'
               outputs(i,k) = max_val(s, sc_top(k), sc_bottom(k), B_equipartition_conv, outputs(i,k))

               i = i+1
               Q_names(i) = 'B_equipartition_conv'
               call r_average(s, sc_top(k), sc_bottom(k), B_equipartition_conv, outputs(i,k))

               i = i+1
               Q_names(i) = 'B_equipartition_pressure'
               call r_average(s, sc_top(k), sc_bottom(k), B_equipartition_pressure, outputs(i,k))

               i = i+1
               Q_names(i) = 'B_diffusion_to_CZ_top'
               outputs(i,k) = compute_B_diffusion_time(s, sc_top(k))

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
           how_many_extra_profile_columns = 5
        end function how_many_extra_profile_columns


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

           names(1) = 'log_F_div_rhocs3'

           do k = 1, nz
              vals(k,1) = (s% L(k)/(4*pi*(s% r(k)**2)))/(s% Rho(k) * s% csound(k)**3)
              vals(k,1) = safe_log10(vals(k,1))
           end do

           names(2) = 'log_Fconv_div_rhocs3'

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

           names(4) = 'nu'
           names(5) = 'alpha'
           do k=1,nz
            call viscosity(s, k, vals(k,4))
            call thermal_diffusivity(s, k, vals(k,5))
         end do

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
