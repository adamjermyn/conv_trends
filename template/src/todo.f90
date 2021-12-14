
           do k=1,nZs
            outputs(i,k) = Pr(k)
           end do
           Q_names(i) = 'Prandtl'; i = i+1
           do k=1,nZs
            outputs(i,k) = Pm(k)
           end do
           Q_names(i) = 'Magnetic_Prandtl'; i = i+1
           do k=1,nZs
            outputs(i,k) = Ra(k)
           end do
           Q_names(i) = 'Rayleigh'; i = i+1

           Q_names(i) = 'Reynolds'; i = i+1


           Q_names(i) = 'Mach'; i = i+1
           Q_names(i) = 'conv_vel'; i = i+1
           Q_names(i) = 'conv_vel_max'; i = i+1
           Q_names(i) = 'f_conv'; i = i+1
           Q_names(i) = 'f_conv_max'; i = i+1
           Q_names(i) = 'buoyant_time'; i = i+1
           Q_names(i) = 'B_equipartition_conv'; i = i+1
           Q_names(i) = 'B_equipartition_conv_max'; i = i+1
           Q_names(i) = 'B_equipartition_conv_surf'; i = i+1
           Q_names(i) = 'B_equipartition_pressure'; i = i+1
           Q_names(i) = 'magnetic_diffusive_time'; i = i+1
           Q_names(i) = 'turnover_time'; i = i+1
           Q_names(i) = 'pressure_scale_height'; i = i+1
           Q_names(i) = 'pressure_scale_height_top'; i = i+1
           Q_names(i) = 'pressure_scale_height_bottom'; i = i+1
           Q_names(i) = 'rho'; i = i+1
           Q_names(i) = 'cz_dtau'; i = i+1
           Q_names(i) = 'cz_top_xm'; i = i+1
           Q_names(i) = 'cz_bottom_xm'; i = i+1
           Q_names(i) = 'cz_top_tau'; i = i+1
           Q_names(i) = 'cz_bottom_tau'; i = i+1
           Q_names(i) = 'dm_penetration'; i = i+1






           ! Identify top of convective core (center has singular values of e.g. density. Use s% nz -1 )
           call get_convective_core(id, sc_convective_core, ierr)
           if (sc_convective_core < s% nz) then
              call get_conv_velocities(id, ierr, v_max_core, v_aver_core, sc_convective_core, &
                   s% nz - 1, b_eq_core,b_max_core,rho_aver_core)
              call get_average_hp(id, ierr, sc_convective_core, s% nz - 1, hp_aver_core)
              call get_turnover(mixing_length_alpha, v_aver_core, hp_aver_core, turnover_core)
              m_core = s% m(sc_convective_core)
              r_core = s% r(sc_convective_core)
              hp_core_top = s% scale_height(sc_convective_core)
              call get_conv_ahp(id, ierr, sc_convective_core, s% nz - 1, v_aver_ahp_core, &
                                mach_top_cz_core, mach_aver_ahp_core, rho_aver_ahp_core)

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
                  call get_conv_velocities(id, ierr, v_HI_max, v_HI_aver, sc_top(k), sc_bottom(k), b_HI_aver, b_HI_max, rho_HI_aver)
                  call get_average_hp(id, ierr, sc_top(k),  sc_bottom(k), HI_hp_aver)
                  call get_conv_ahp(id, ierr, sc_top(k),  sc_bottom(k), v_HI_aver_ahp, mach_HI_top, mach_HI_aver_ahp, rho_HI_aver)
                  call get_microturb(mach_HI_aver_ahp, rho_HI_aver, rho_surf,v_HI_aver_ahp, v_HI_surf)
                  call get_turnover(mixing_length_alpha, v_HI_aver, HI_hp_aver, turnover_HI)
                  call get_bsurf(rho_surf, rho_HI_aver, b_HI_aver, b_HI_max, b_HI_surf, b_HI_surf_max)
                  call get_conv_radii(id, ierr, sc_top(k), sc_bottom(k), HI_r_top, HI_r_bottom)
                  call get_conv_mass(id, ierr, sc_top(k), sc_bottom(k), HI_mass)
                  call get_max_fc(id, ierr, HI_fcmax, sc_top(k), sc_bottom(k))
                  call get_pressure_eq_field(id, ierr, sc_top(k), sc_bottom (k),HI_b_p_eq,HI_b_p_max)
                  call compute_Ra_Re(id, ierr, sc_top(k), sc_bottom(k), v_HI_aver, HI_Ra, HI_Re, HI_Pr, HI_nu, HI_alpha, HI_eta, HI_dz)
                  HI_tau_eta = compute_B_diffusion_time(s,  sc_top(k))
                  HI_buoyant_time = compute_buoyant_time(s,  sc_top(k), b_HI_aver)
                  HI_xm = sum(s%dm(1:sc_top(k))) / Msun
                  !write(*,*) sc_top(k), sc_bottom(k), sc_type(k), v_HI_aver, rho_HI_aver, b_HI_max, b_HI_surf, v_HI_surf
               else if ( sc_type(k) == 'HeI' ) then
                  call get_conv_velocities(id, ierr, v_HeI_max, v_HeI_aver, sc_top(k), sc_bottom(k), &
                  b_HeI_aver, b_HeI_max, rho_HeI_aver)
                  call get_average_hp(id, ierr, sc_top(k),  sc_bottom(k), HeI_hp_aver)
                  call get_conv_ahp(id, ierr, sc_top(k),  sc_bottom(k), v_HeI_aver_ahp, mach_HeI_top, &
                  mach_HeI_aver_ahp, rho_HeI_aver)
                  call get_microturb(mach_HeI_aver_ahp, rho_HeI_aver, rho_surf,v_HeI_aver_ahp, v_HeI_surf)
                  call get_turnover(mixing_length_alpha, v_HeI_aver, HeI_hp_aver, turnover_HeI)
                  call get_bsurf(rho_surf, rho_HeI_aver, b_HeI_aver, b_HeI_max, b_HeI_surf, b_HeI_surf_max)
                  call get_conv_radii(id, ierr, sc_top(k), sc_bottom(k), HeI_r_top, HeI_r_bottom)
                  call get_conv_mass(id, ierr, sc_top(k), sc_bottom(k), HeI_mass)
                  call get_max_fc(id, ierr, HeI_fcmax, sc_top(k), sc_bottom(k))
                  call get_pressure_eq_field(id, ierr, sc_top(k), sc_bottom (k),HeI_b_p_eq,HeI_b_p_max)
                  call compute_Ra_Re(id, ierr, sc_top(k), sc_bottom(k), v_HeI_aver, HeI_Ra, HeI_Re, HeI_Pr, HeI_nu, HeI_alpha, HeI_eta, HeI_dz)
                  HeI_tau_eta = compute_B_diffusion_time(s,  sc_top(k))
                  HeI_buoyant_time = compute_buoyant_time(s,  sc_top(k), b_HeI_aver)
                  HeI_xm = sum(s%dm(1:sc_top(k))) / Msun
                  !write(*,*) sc_top(k), sc_bottom(k), sc_type(k), v_HeI_aver, rho_HeI_aver, b_HeI_max, b_HeI_surf, v_HeI_surf
               else if ( sc_type(k) == 'HeII' ) then
                  call get_conv_velocities(id, ierr, v_HeII_max, v_HeII_aver, sc_top(k), sc_bottom(k), &
                  b_HeII_aver, b_HeII_max, rho_HeII_aver)
                  call get_average_hp(id, ierr, sc_top(k),  sc_bottom(k), HeII_hp_aver)
                  call get_conv_ahp(id, ierr, sc_top(k),  sc_bottom(k), v_HeII_aver_ahp, mach_HeII_top, &
                  mach_HeII_aver_ahp, rho_HeII_aver)
                  call get_microturb(mach_HeII_aver_ahp, rho_HeII_aver, rho_surf,v_HeII_aver_ahp, v_HeII_surf)
                  call get_turnover(mixing_length_alpha, v_HeII_aver, HeII_hp_aver, turnover_HeII)
                  call get_bsurf(rho_surf, rho_HeII_aver, b_HeII_aver, b_HeII_max, b_HeII_surf, b_HeII_surf_max)
                  call get_conv_radii(id, ierr, sc_top(k), sc_bottom(k), HeII_r_top, HeII_r_bottom)
                  call get_conv_mass(id, ierr, sc_top(k), sc_bottom(k), HeII_mass)
                  call get_max_fc(id, ierr, HeII_fcmax, sc_top(k), sc_bottom(k))
                  call get_pressure_eq_field(id, ierr, sc_top(k), sc_bottom (k),HeII_b_p_eq,HeII_b_p_max)
                  call compute_Ra_Re(id, ierr, sc_top(k), sc_bottom(k), v_HeII_aver, HeII_Ra, HeII_Re, HeII_Pr, HeII_nu, HeII_alpha, HeII_eta, HeII_dz)
                  HeII_tau_eta = compute_B_diffusion_time(s,  sc_top(k))
                  HeII_buoyant_time = compute_buoyant_time(s,  sc_top(k), b_HeII_aver)
                  HeII_xm = sum(s%dm(1:sc_top(k))) / Msun
                  !write(*,*) sc_top(k), sc_bottom(k), sc_type(k), v_HeII_aver, rho_HeII_aver, b_HeII_max, b_HeII_surf, v_HeII_surf
               else if ( sc_type(k) == 'FeCZ' ) then
                 call get_conv_velocities(id, ierr, v_FeCZ_max, v_FeCZ_aver, sc_top(k), sc_bottom(k), &
                  b_FeCZ_aver, b_FeCZ_max, rho_FeCZ_aver)
                  call get_average_hp(id, ierr, sc_top(k),  sc_bottom(k), FeCZ_hp_aver)
                  call get_conv_ahp(id, ierr, sc_top(k),  sc_bottom(k), v_FeCZ_aver_ahp, mach_FeCZ_top, &
                  mach_FeCZ_aver_ahp, rho_FeCZ_aver)
                  call get_microturb(mach_FeCZ_aver_ahp, rho_FeCZ_aver, rho_surf,v_FeCZ_aver_ahp, v_FeCZ_surf)
                  call get_turnover(mixing_length_alpha, v_FeCZ_aver, FeCZ_hp_aver, turnover_FeCZ)
                  call get_bsurf(rho_surf, rho_FeCZ_aver, b_FeCZ_aver, b_FeCZ_max, b_FeCZ_surf, b_FeCZ_surf_max)
                  call get_conv_radii(id, ierr, sc_top(k), sc_bottom(k), FeCZ_r_top, FeCZ_r_bottom)
                  call get_conv_mass(id, ierr, sc_top(k), sc_bottom(k), FeCZ_mass)
                  call get_max_fc(id, ierr, FeCZ_fcmax, sc_top(k), sc_bottom(k))
                  call get_pressure_eq_field(id, ierr, sc_top(k), sc_bottom (k),FeCZ_b_p_eq,FeCZ_b_p_max)
                  call compute_Ra_Re(id, ierr, sc_top(k), sc_bottom(k), v_FeCZ_aver, FeCZ_Ra, FeCZ_Re, FeCZ_Pr, FeCZ_nu, FeCZ_alpha, FeCZ_eta, FeCZ_dz)
                  FeCZ_tau_eta = compute_B_diffusion_time(s,  sc_top(k))
                  FeCZ_buoyant_time = compute_buoyant_time(s,  sc_top(k), b_FeCZ_aver)
                  FeCZ_xm = sum(s%dm(1:sc_top(k))) / Msun
                  !write(*,*) sc_top(k), sc_bottom(k), sc_type(k), v_FeCZ_aver, rho_FeCZ_aver, b_FeCZ_max, b_FeCZ_surf, v_FeCZ_surf
               end if
            end if
           end do