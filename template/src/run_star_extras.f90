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
            vals(k,1) = 15d0 * log10(s%brunt_N2(k) + 1d-18)
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
           how_many_extra_history_columns = 368
        end function how_many_extra_history_columns


        subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
           use identification
           use diffusivities
           use utilities
           use magnetism
           use penetration
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
           integer, parameter :: nQs = 60
           integer, parameter :: nZs = 6 ! Max # of CZs
           integer :: nFound
           logical :: sc_exists(nZs)
           integer, dimension(nZs) :: sc_top, sc_bottom
           character(len=100) :: Q_names(nQs)
           character(len=100) :: cz_names(nZs)
           character(len=100) :: cz_name
           real(dp) :: outputs(nQs, nZs)
           real(dp) :: Ra(nZs), Pm(nZs), Pr(nZs), nu(nZs), alpha(nZs), eta(nZs), dr(nZs)

           ! Quantities that are one per star
           integer, parameter :: nOffs = 8
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
           cz_names(i) = 'Envelope'; i = i+1

           ! Identify which CZs, if any, exist
           sc_exists = .false.
           sc_top = -1
           sc_bottom = -1
           call get_conv_regions(s, nZs, nFound)
           do i=1,nFound
            call classify_conv_region(s, s% mixing_region_top(i), s% mixing_region_bottom(i), cz_name)
            do k=1,nZs
               if (cz_names(k) == cz_name) then
                  ! These are the indices of the inner-most and outer-most cells of each cz.
                  ! That means the cz's lie between the faces (top, bottom-1)
                  sc_top(k) = s% mixing_region_top(i)
                  sc_bottom(k) = s% mixing_region_bottom(i)
                  sc_exists(k) = .true.
               end if
            end do
           end do

           ! Quantities in the CZs
           ! Quantities with profiles are radially averaged, weighted by dr unless otherwise specified.  
           do k=1,nZs
            if (sc_exists(k)) then
               i = 0

               i = i+1
               Q_names(i) = 'bottom_m_div_mstar'
               outputs(i,k) = s%m(sc_bottom(k))/s%m(1)

               i = i+1
               Q_names(i) = 'top_m_div_mstar'
               outputs(i,k) = s%m(sc_top(k))/s%m(1)

               i = i+1
               Q_names(i) = 'bottom_index'
               outputs(i,k) = sc_bottom(k)

               i = i+1
               Q_names(i) = 'top_index'
               outputs(i,k) = sc_top(k)

               i = i+1
               Q_names(i) = 'viscosity'
               call r_average(s, sc_top(k), sc_bottom(k), viscosity, outputs(i,k))

               i = i+1
               Q_names(i) = 'thermal_diffusivity'
               call r_average(s, sc_top(k), sc_bottom(k), thermal_diffusivity, outputs(i,k))

               i = i+1
               Q_names(i) = 'prandtl'
               call r_average(s, sc_top(k), sc_bottom(k), prandtl, outputs(i,k))

               i = i+1
               Q_names(i) = 'magnetic_diffusivity'
               call r_average(s, sc_top(k), sc_bottom(k), magnetic_diffusivity, outputs(i,k))

               i = i+1
               Q_names(i) = 'pressure_scale_height'
               call r_average(s, sc_top(k), sc_bottom(k), scale_height, outputs(i,k))

               i = i+1
               Q_names(i) = 'pressure_scale_height_top'
               outputs(i,k) = s%scale_height(sc_top(k))

               i = i+1
               Q_names(i) = 'pressure_scale_height_bottom'
               outputs(i,k) = s%scale_height(sc_bottom(k))

               i = i+1
               Q_names(i) = 'Prad_div_P'
               call r_average(s, sc_top(k), sc_bottom(k), Prad_div_P, outputs(i,k))

               i = i+1
               Q_names(i) = 'density'
               call r_average(s, sc_top(k), sc_bottom(k), density, outputs(i,k))

               i = i+1
               Q_names(i) = 'density_ratio'
               outputs(i,k) = s%rho(sc_bottom(k))/s%rho(sc_top(k))

               i = i+1
               Q_names(i) = 'sound_speed_adiabatic'
               call r_average(s, sc_top(k), sc_bottom(k), sound_speed_adiabatic, outputs(i,k))

               i = i+1
               Q_names(i) = 'sound_speed_isothermal'
               call r_average(s, sc_top(k), sc_bottom(k), sound_speed_isothermal, outputs(i,k))

               i = i+1
               Q_names(i) = 'gravity'
               call r_average(s, sc_top(k), sc_bottom(k), gravity, outputs(i,k))

               i = i+1
               Q_names(i) = 'gradr_sub_grada'
               call r_average(s, sc_top(k), sc_bottom(k), gradr_sub_grada, outputs(i,k))

               i = i+1
               Q_names(i) = 'conv_vel' ! This is a face quantity defined on the inside of the CZ, so want to go from (top+1...bottom-1)
               call r_average(s, sc_top(k)+1, sc_bottom(k)-1, conv_vel, outputs(i,k))

               i = i+1
               Q_names(i) = 'conv_vel_div_r' ! This is a face quantity defined on the inside of the CZ, so want to go from (top+1...bottom-1)
               call r_average(s, sc_top(k)+1, sc_bottom(k)-1, conv_vel_div_r, outputs(i,k))

               i = i+1
               Q_names(i) = 'L_conv_div_L' ! This is a face quantity defined on the inside of the CZ, so want to go from (top+1...bottom-1)
               call r_average(s, sc_top(k)+1, sc_bottom(k)-1, L_conv_div_L, outputs(i,k))

               i = i+1
               Q_names(i) = 'conv_vel_max' ! This is a face quantity defined on the inside of the CZ, so want to go from (top+1...bottom-1)
               call max_val(s, sc_top(k)+1, sc_bottom(k)-1, conv_vel, outputs(i,k))

               i = i+1
               Q_names(i) = 'stiffness_bottom'
               call stiffness_bottom(s, sc_top(k), sc_bottom(k), outputs(i,k))

               i = i+1
               Q_names(i) = 'stiffness_top'
               call stiffness_top(s, sc_top(k), sc_bottom(k), outputs(i,k))

               i = i+1
               Q_names(i) = 'cz_dm'
               call integrate_dm(s, sc_top(k), sc_bottom(k), unity, outputs(i,k))

               i = i+1
               Q_names(i) = 'cz_dr'
               call integrate_dr(s, sc_top(k), sc_bottom(k), unity, outputs(i,k))

               i = i+1
               Q_names(i) = 'turnover_time'
               call turnover(s, sc_top(k), sc_bottom(k), outputs(i,k))

               i = i+1
               Q_names(i) = 'cz_top_r'
               outputs(i,k) = s%r(sc_top(k))

               i = i+1
               Q_names(i) = 'cz_bottom_r'
               call bottom_r(s, sc_bottom(k), outputs(i,k))

               i = i+1
               Q_names(i) = 'cz_dtau'
               call integrate_dr(s, sc_top(k), sc_bottom(k), dtau_dr, outputs(i,k))

               i = i+1
               Q_names(i) = 'cz_top_xm' ! Integrate from cell 1 down to cell top-1 because we want xm of the top face. 
               call integrate_dm(s, 1, sc_top(k)-1, unity, outputs(i,k))

               i = i+1
               Q_names(i) = 'cz_bottom_xm' ! Integrate from cell 1 down to cell bottom because we want xm of the bottom face (which is face bottom+1). 
               call integrate_dm(s, 1, sc_bottom(k), unity, outputs(i,k))

               i = i+1
               Q_names(i) = 'cz_top_tau'
               outputs(i,k) = s%tau(sc_top(k))

               i = i+1
               Q_names(i) = 'cz_bottom_tau'
               call bottom_tau(s,k,outputs(i,k))

               i = i+1
               Q_names(i) = 'B_shutoff'! This is a face quantity defined on the inside of the CZ, so want to go from (top+1...bottom-1)
               call max_val(s, sc_top(k)+1, sc_bottom(k)-1, B_shutoff_conv, outputs(i,k))

               i = i+1
               Q_names(i) = 'B_equipartition_conv_max'! This is a face quantity defined on the inside of the CZ, so want to go from (top+1...bottom-1)
               call max_val(s, sc_top(k)+1, sc_bottom(k)-1, B_equipartition_conv, outputs(i,k))

               i = i+1
               Q_names(i) = 'B_equipartition_conv'! This is a face quantity defined on the inside of the CZ, so want to go from (top+1...bottom-1)
               call r_average(s, sc_top(k)+1, sc_bottom(k)-1, B_equipartition_conv, outputs(i,k))

               i = i+1
               Q_names(i) = 'B_equipartition_pressure'
               call r_average(s, sc_top(k), sc_bottom(k), B_equipartition_pressure, outputs(i,k))

               i = i+1
               Q_names(i) = 'B_diffusion_to_CZ_top'
               call compute_B_diffusion_time(s, sc_top(k), outputs(i,k)) ! This one handles going to the outer face inside the call.

               i = i+1
               Q_names(i) = 'L_div_Ledd'! This is a face quantity but it's defined perfectly well on the edges of the CZ, so we go from (top,bottom+1)
                                         ! L=0 at the center so we can just exclude that if the CZ runs all the way to the center.
               call r_average(s, sc_top(k), min(s%nz,sc_bottom(k)+1), L_div_Ledd, outputs(i,k))

               i = i+1
               Q_names(i) = 'Nusselt'! This is a face quantity defined on the inside of the CZ, so want to go from (top+1...bottom-1)
               call r_average(s, sc_top(k)+1, sc_bottom(k)-1, Nusselt, outputs(i,k))

               i = i+1
               Q_names(i) = 'L_div_Ledd_max'! This is a face quantity but it's defined perfectly well on the edges of the CZ, so we go from (top,bottom+1).
                                            ! L=0 at the center so we can just exclude that if the CZ runs all the way to the center.
               call max_val(s, sc_top(k), min(s%nz,sc_bottom(k)+1), L_div_Ledd, outputs(i,k))

               i = i+1
               Q_names(i) = 'viscosity_bottom'
               call r_average_hp_frac(s, sc_top(k), sc_bottom(k), 0.3d0, .true., .false., viscosity, outputs(i,k))

               i = i+1
               Q_names(i) = 'viscosity_top'
               call r_average_hp_frac(s, sc_top(k), sc_bottom(k), 0.3d0, .true., .true., viscosity, outputs(i,k))

               i = i+1
               Q_names(i) = 'thermal_diffusivity_bottom'
               call r_average_hp_frac(s, sc_top(k), sc_bottom(k), 0.3d0, .true., .false., thermal_diffusivity, outputs(i,k))

               i = i+1
               Q_names(i) = 'thermal_diffusivity_top'
               call r_average_hp_frac(s, sc_top(k), sc_bottom(k), 0.3d0, .true., .true., thermal_diffusivity, outputs(i,k))

               i = i+1
               Q_names(i) = 'magnetic_diffusivity_bottom'
               call r_average_hp_frac(s, sc_top(k), sc_bottom(k), 0.3d0, .true., .false., magnetic_diffusivity, outputs(i,k))

               i = i+1
               Q_names(i) = 'magnetic_diffusivity_top'
               call r_average_hp_frac(s, sc_top(k), sc_bottom(k), 0.3d0, .true., .true., magnetic_diffusivity, outputs(i,k))

               i = i+1
               Q_names(i) = 'Prad_div_P_bottom'
               call r_average_hp_frac(s, sc_top(k), sc_bottom(k), 0.3d0, .true., .false., Prad_div_P, outputs(i,k))

               i = i+1
               Q_names(i) = 'Prad_div_P_top'
               call r_average_hp_frac(s, sc_top(k), sc_bottom(k), 0.3d0, .true., .true., Prad_div_P, outputs(i,k))

               i = i+1
               Q_names(i) = 'sound_speed_adiabatic_bottom'
               call r_average_hp_frac(s, sc_top(k), sc_bottom(k), 0.3d0, .true., .false., sound_speed_adiabatic, outputs(i,k))

               i = i+1
               Q_names(i) = 'sound_speed_adiabatic_top'
               call r_average_hp_frac(s, sc_top(k), sc_bottom(k), 0.3d0, .true., .true., sound_speed_adiabatic, outputs(i,k))

               i = i+1
               Q_names(i) = 'conv_vel_bottom'
               call r_average_hp_frac(s, sc_top(k), sc_bottom(k), 0.3d0, .true., .false., conv_vel, outputs(i,k))

               i = i+1
               Q_names(i) = 'conv_vel_top'
               call r_average_hp_frac(s, sc_top(k), sc_bottom(k), 0.3d0, .true., .true., conv_vel, outputs(i,k))

               i = i+1
               Q_names(i) = 'L_conv_div_L_bottom'
               call r_average_hp_frac(s, sc_top(k), sc_bottom(k), 0.3d0, .true., .false., L_conv_div_L, outputs(i,k))

               i = i+1
               Q_names(i) = 'L_conv_div_L_top'
               call r_average_hp_frac(s, sc_top(k), sc_bottom(k), 0.3d0, .true., .true., L_conv_div_L, outputs(i,k))

               i = i+1
               Q_names(i) = 'L_div_Ledd_bottom'
               call r_average_hp_frac(s, sc_top(k), sc_bottom(k), 0.3d0, .true., .false., L_div_Ledd, outputs(i,k))

               i = i+1
               Q_names(i) = 'L_div_Ledd_top'
               call r_average_hp_frac(s, sc_top(k), sc_bottom(k), 0.3d0, .true., .true., L_div_Ledd, outputs(i,k))

               i = i+1
               Q_names(i) = 'gravity_bottom'
               call r_average_hp_frac(s, sc_top(k), sc_bottom(k), 0.3d0, .true., .false., gravity, outputs(i,k))

               i = i+1
               Q_names(i) = 'gravity_top'
               call r_average_hp_frac(s, sc_top(k), sc_bottom(k), 0.3d0, .true., .true., gravity, outputs(i,k))

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

           call B_shutoff_down_to_tau(s,3d0,one_off_outputs(i))
           one_off_names(i) = 'B_shutoff_down_to_tau_3'; i=i+1

           call B_shutoff_down_to_tau(s,1d1,one_off_outputs(i))
           one_off_names(i) = 'B_shutoff_down_to_tau_10'; i=i+1

           call B_shutoff_down_to_tau(s,1d2,one_off_outputs(i))
           one_off_names(i) = 'B_shutoff_down_to_tau_100'; i=i+1

           call B_shutoff_down_to_tau(s,1d3,one_off_outputs(i))
           one_off_names(i) = 'B_shutoff_down_to_tau_1000'; i=i+1

           call B_shutoff_down_to_tau(s,1d4,one_off_outputs(i))
           one_off_names(i) = 'B_shutoff_down_to_tau_1d4'; i=i+1


           ! Pack output 
           k = 1
           do i=1,nQs
            do j=1,nZs
               name = trim(cz_names(j)) // '_' // trim(Q_names(i))
               names(k) = name
               vals(k) = outputs(i,j)
               k = k + 1
            end do
           end do

           do i=1,nOffs
            names(k) = one_off_names(i)
            vals(k) = one_off_outputs(i)
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
           how_many_extra_profile_columns = 6
        end function how_many_extra_profile_columns


        subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
           use diffusivities
           use magnetism
           use star_def, only: star_info, maxlen_profile_column_name
           use const_def, only: dp
           use math_lib, only: safe_log10
           integer, intent(in) :: id, n, nz
           character (len=maxlen_profile_column_name) :: names(n)
           real(dp) :: vals(nz,n) 
           integer, intent(out) :: ierr
           type (star_info), pointer :: s
           integer :: k
           ierr = 0
           call star_ptr(id, s, ierr)
           if (ierr /= 0) return

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
           names(6) = 'B_shutoff'
           do k=1,nz
            call viscosity(s, k, vals(k,4))
            call thermal_diffusivity(s, k, vals(k,5))
            call B_shutoff_conv(s, k, vals(k,6))
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

      end module run_star_extras
