
           Q_names(i) = 'dm_penetration'; i = i+1


           ! Identify top of convective core (center has singular values of e.g. density. Use s% nz -1 )
           call get_convective_core(id, sc_convective_core, ierr)
           if (sc_convective_core < s% nz) then
              call get_conv_velocities(id, ierr, v_max_core, v_aver_core, sc_convective_core, &
                   s% nz - 1, b_eq_core,b_max_core,rho_aver_core)
              call get_average_hp(id, ierr, sc_convective_core, s% nz - 1, hp_aver_core)
              call get_turnover(mixing_length_alpha, v_aver_core, hp_aver_core, turnover_core)
           end if

           call compute_dm_core(s, m_conv_core, dm_core, dr_core, dr_core_div_h)

           ! Identify number of convective regions above a certain temperature  (Max 4, HI, HeI, HeII, FeCZ)

           call get_conv_regions_above_T(id,1d6,ierr,num_conv_regions)
           names(1) = 'subsurface_convective_regions'
           vals(1)  = num_conv_regions

           ! Calculate relevant column values
           do k = 1, num_conv_regions ! num_conv_regions should always be >= 1
             sc_top(k) = s% mixing_region_top(k)
             sc_bottom(k) = s% mixing_region_bottom(k)
             if (sc_top(k) .NE. 0) then
               call classify_conv_region_above_T(id, ierr, sc_top(k), sc_bottom(k), sc_type(k))
               if ( sc_type(k) == 'HI' ) then
                  HI_tau_eta = compute_B_diffusion_time(s,  sc_top(k))
                  HI_buoyant_time = compute_buoyant_time(s,  sc_top(k), b_HI_aver)
                  HI_xm = sum(s%dm(1:sc_top(k))) / Msun
                  !write(*,*) sc_top(k), sc_bottom(k), sc_type(k), v_HI_aver, rho_HI_aver, b_HI_max, b_HI_surf, v_HI_surf
              end if
            end if
           end do