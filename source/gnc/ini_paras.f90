subroutine apply_paras()
	use model_config
	use model_basic
	use md_star_pot
	use md_mbh_evl_acc
	implicit none
	integer i,ier
	character*(50) str
	character*(50) pname
	character*(200) pvalue
	character*(400) root_path
	do i=1, pa_now_used%n
		pname=trim(adjustl(pa_now_used%tp(i)%name))
		pvalue=trim(adjustl(pa_now_used%tp(i)%str))
		associate(tp=>pa_now_used%tp(i))
			select case(trim(adjustl(pname)))
			case("version number")
				read(unit=pvalue,fmt=*) ctl%current_version_number
			case("all exact in sample para")
				read(unit=pvalue,fmt=*) ctl%all_exact
			case("energy bin type")
				read(unit=pvalue,fmt=*) ctl%str_barge_grid_type
				select case(trim(adjustl(ctl%str_barge_grid_type)))
				case("IREG_P")
					ctl%barge_grid_type=barge_grid_type_iregular_phi
				case("IREG_J")
					ctl%barge_grid_type=barge_grid_type_iregular_jc
					call tp%get_sub_para("--energy max min steps", str,ier) 
					read(unit=str,fmt=*) gx_func_max_step,gx_func_min_step
				case default
					print*, "error ! define barge_grid type:",&
						trim(adjustl(ctl%str_barge_grid_type))
					stop
				end select
			case("barg evaluate method")
				read(unit=pvalue,fmt=*) ctl%str_barge_evl_method
				select case(trim(adjustl(ctl%str_barge_evl_method))) 
				case("GRID_2D")
					ctl%barge_evl_method=barge_evl_method_grid_2d
				case default
					print*, "error! define barge evel method: ", ctl%str_barge_evl_method
					stop
				end select
				
			case("diffusion coefficient bin type")
				read(unit=pvalue,fmt=*) ctl%str_dc_grid_type
				select case(trim(adjustl(ctl%str_dc_grid_type))) 
				case("IREG")
					ctl%dc_grid_type=dc_grid_irregular
				case default
					print*, "error ! define dc_grid type:", trim(adjustl(ctl%str_dc_grid_type))
					stop
				end select
			case("diffusion coefficient integral method")
				read(unit=pvalue,fmt=*) ctl%str_dc_method
				!ctl%str_dc_method=trim(adjustl(str_dc_method))
				select case(trim(adjustl( ctl%str_dc_method)))
				case("FAST")
					ctl%get_dc_method=get_dc_method_fast
				case("GENERAL")
					ctl%get_dc_method=get_dc_method_general
				case default
					print*, "error ! define dc_grid type:", trim(adjustl( ctl%str_dc_method))
					stop
				end select
			case("diffusion coefficient integral acc")
				read(unit=pvalue,fmt=*) star_dc_int_acc_a,star_dc_int_acc_r			
			
			case("aux bin sizes")
				read(unit=pvalue,fmt=*) aux_function_bin_size, f12_function_bin_size_1, f12_function_bin_size_2
			case("func r and pd int acc")
				read(unit=pvalue,fmt=*) fr_funcs_int_acc_a, fr_funcs_int_acc_r, &
						pd_int_acc_a, pd_int_acc_r
			case("adiabatic invariant method")
				read(unit=pvalue,fmt=*) ctl%str_adb_est_method
				select case(trim(adjustl(ctl%str_adb_est_method)))
				case("FAST")
					ctl%adb_est_method=adb_est_method_fast
				case("ACC")
					ctl%adb_est_method=adb_est_method_acc
				case default
					print*, "error! define adb est method", trim(adjustl(ctl%str_adb_est_method))
					stop
				end select
			case("mbh mass growth")
				read(unit=pvalue, fmt=*) ctl%enable_evl_mbh
				if(ctl%enable_evl_mbh.ge.1)then
					call tp%get_sub_para("--initial gas reservoir", str,ier)
					!print*, "get sub str=", str
					read(unit=str,fmt=*) mbh_mmg%gas_reservior_left
				end if
			case("init adiabatically response of mbh")
				read(unit=pvalue, fmt=*) ctl%init_adb_mbh_inc
				if(ctl%init_adb_mbh_inc.ge.1)then
					call tp%get_sub_para("---init mbh log mass",str,ier)
					read(unit=str,fmt=*) ctl%init_adb_mbh_log_factor
					call tp%get_sub_para("---n steps",str,ier)
					read(unit=str,fmt=*) ctl%ini_adb_increase_nt
					call tp%get_sub_para("---target mbh log mass",str,ier)
					read(unit=str,fmt=*) ctl%init_adb_mbh_log_factor_target
				end if
			case("energy and angular momentum in bins of")
				read(unit=pvalue,fmt=*) ctl%str_dejmodel
				select case (trim(adjustl(ctl%str_dejmodel)))
				case ("EJ")
					ctl%dejmodel=dejmodel_EJ 
				case default
					print*, "ERROR, define dejmodel:", trim(adjustl(ctl%str_dejmodel))
					stop
				end select
			case("interpolate method of energy and J ")
				read(unit=pvalue,fmt=*) ctl%str_method_interpolate
				select case(trim(adjustl(ctl%str_method_interpolate)))
				case("near")
					ctl%method_interpolate=method_int_nearst
					print*, "warnning, you are using a near interpolot method, which could be inaccurate, use li2d instead"
				case("li2d")
					ctl%method_interpolate=method_int_linear
				case default
					print*, "ERROR, define interpolation method:", trim(adjustl(ctl%str_method_interpolate))
				end select
			case("minimum and maximum Jdmless")
				read(unit=pvalue,fmt=*) jmin_value, jmax_value
				log10jmin_value=log10(jmin_value)
				log10jmax_value=log10(jmax_value)
			case("relativistic loss cone table dir")
				read(unit=pvalue,fmt="(A200)") ctl%rv_nw_conv_dir
				call GETENV("GWPATH",root_path)
				ctl%rv_nw_conv_dir=trim(adjustl(root_path))//trim(adjustl(ctl%rv_nw_conv_dir))
			case("mass unit")
				read(unit=pvalue,fmt=*) m0_cl
				call get_rh_in_pc(m0_cl,r0_cl)
				r0_cl=r0_cl*pc	
				!print*, "apply_paras:",m0_cl,r0_cl
			case("mbh mass")
				read(unit=pvalue,fmt=*) spp_new%mbh_dmless
				if(spp_new%mbh_dmless.eq.0)then
					call tp%get_sub_para("--replace sample exceed energy max",str,ier)
					read(unit=str,fmt=*) ctl%replace_sample_eceed_emax
				end if
			case("mbh spin")
				read(unit=pvalue,fmt=*) mbh_spin 
			case("rmin rmax")
				read(unit=pvalue,fmt=*) ctl%rmin_factor, ctl%rmax_factor
				ctl%log10rmin_factor=log10(ctl%rmin_factor)
				ctl%log10rmax_factor=log10(ctl%rmax_factor)
				!print*, ctl%rmin_factor,ctl%rmax_factor, pvalue
			case("set emax")
				read(unit=pvalue,fmt=*) ctl%emax_boundary
			case("bin type of angular momentum")
				read(unit=pvalue,fmt=*) ctl%str_jbin_bd
				select case(trim(adjustl(ctl%str_jbin_bd))) 
				case("LOG")
					ctl%jbin_type=Jbin_type_log 
				case default
					print*, "error! define jbin type", ctl%str_jbin_bd
					stop
				end select
			case("init type of f(j)")
				read(unit=pvalue,fmt=*) ctl%str_fj_bd
				select case(trim(adjustl(ctl%str_fj_bd)))
				case("ISO")
					ctl%boundary_fj=boundary_fj_iso 
				case default
					print*, "error! define fj type", ctl%str_fj_bd
					stop
				end select
			case("bin type of energy")
				read(unit=pvalue,fmt=*) ctl%str_ebin_type
				select case(trim(adjustl(ctl%str_ebin_type))) 
				case("LOG")
					ctl%ebin_type=ebin_type_log
				case default
					print*, "error! define ebin type", ctl%str_ebin_type
					stop
				end select
			case("init random seed")
				read(unit=pvalue,fmt=*) ctl%same_rseed_ini
			case("running random seed")
				read(unit=pvalue,fmt=*) ctl%seed_value
			case("bin number of f(r)")
				read(unit=pvalue,fmt=*) ctl%dstr_bins_r
			case("bin number of f(e)") 
				read(unit=pvalue,fmt=*)  ctl%dstr_bins_e
			case("bin number of f(j)") 
				read(unit=pvalue,fmt=*)  ctl%dstr_bins_j
			case("bin number of dc(e,j)") 
				read(unit=pvalue,fmt=*)  ctl%diff_coeff_bins
				if(mod(ctl%diff_coeff_bins,ctl%ntasks).ne.0)then
					print*, "diff_coeff_bins, ntasks=", ctl%diff_coeff_bins, ctl%ntasks
					print*, "warnning! grid_bin number should be integer times of ntasks"
				end if
			case("use same running random seed")
				read(unit=pvalue,fmt=*)  ctl%same_rseed_evl
			case("minimum number of samples in mass bin")
				read(unit=pvalue,fmt=*) ctl%min_sample_in_mass_bin
			case("task mode")
				read(unit=pvalue,fmt=*) ctl%str_task_mode
				select case (trim(adjustl(ctl%str_task_mode)))
				case ("NEW")
					ctl%intaskmode=task_mode_new
					ctl%ntask_bg=0
				case ("APP")
					ctl%intaskmode=task_mode_append
					call tp%get_sub_para("--task begin index",str,ier)
					read(unit=str,fmt=*) ctl%ntask_bg

				case default
					print*, "define taskmode", trim(adjustl(ctl%str_task_mode))
					stop
				end select
				ctl%ntask_total=ctl%ntask_bg+ctl%ntasks
			case("snapshot mode")
				read(unit=pvalue,fmt=*) ctl%str_snap_mode

				select case (trim(adjustl(ctl%str_snap_mode)))
				case ("NEW")
					ctl%insnapmode=snap_mode_new
				case ("APP")
					ctl%insnapmode=snap_mode_append
					call tp%get_sub_para("--snapshot begin index",str,ier)
					read(unit=str,fmt=*) ctl%n_spshot_bg
				case ("ONE")
					ctl%insnapmode=snap_mode_one
				case default
					print*, "define snapmode ", trim(adjustl(ctl%str_snap_mode))
					stop
				end select  
			case("num of particle update per snap")
				read(unit=pvalue,fmt=*) ctl%num_step_per_update
			case("timestep unit")
				read(unit=pvalue,fmt=*) ctl%time_unit
				select case(trim(adjustl(ctl%time_unit)))
				case("TNR")
					ctl%timestep_mode=timestep_mode_tnr
				case("TRH")
					ctl%timestep_mode=timestep_mode_trh
				end select
				call tp%get_sub_para("--update 2body relaxation timestep",str,ier)
				read(unit=str,fmt=*) ctl%update_2bdrlx_per_shapshot
			case("fraction of timestep of snapshot")
				read(unit=pvalue,fmt=*) ctl%ts_snap_dt_per_snap
			case("max timestep of snapshot")
				read(unit=pvalue,fmt=*) ctl%tmax_timestep
			case ("running mode of time")
				read(unit=pvalue,fmt=*) ctl%str_end_time_mode
				select case(trim(adjustl(ctl%str_end_time_mode)))
				case("SNAP")
					call tp%get_sub_para("--num of total snapshots",str,ier)
					read(unit=str, fmt=*) ctl%n_spshot
					ctl%n_spshot_total=ctl%n_spshot+ctl%n_spshot_bg
					ctl%time_run_mode=time_run_mode_snap

				case("TTOT")
					call tp%get_sub_para("--total time (Myr)",str,ier)
					read(unit=str, fmt=*) ctl%total_time
					ctl%time_run_mode=time_run_mode_ttot
				end select
			case("freq of data output")
				read(unit=pvalue,fmt=*) ctl%output_data_freq
			case("freq of dms output")
				read(unit=pvalue,fmt=*) ctl%output_dms_freq
			case("two body relaxation")
				read(unit=pvalue,fmt=*) ctl%two_body_relaxation_on
			case("convergence of potential critical value")
				read(unit=pvalue,fmt=*)  ctl%gx_conv_cri 
			case("max iteration of potential convergence")
				read(unit=pvalue,fmt=*)  max_self_con_iter
			case("alpha cri")
				read(unit=pvalue,fmt=*) ctl%md_alpha_cri
			case("fine rtd range")
				read(unit=pvalue,fmt=*) ctl%md_fine_rtd_range
			case("GW outer orbit")
				read(unit=pvalue, fmt=*) ctl%gw_radiation_otby
				if(ctl%gw_radiation_otby.ge.1)then
					call tp%get_sub_para("--EMRI criteria",str,ier)
					read(unit=str,fmt=*) ctl%emri_cri
					call tp%get_sub_para("--output EMRI data bin",str,ier)
					read(unit=str,fmt=*) ctl%output_gw_emri_data_bin
				end if
			case("density estimation method")
				read(unit=pvalue, fmt=*) ctl%str_source_den
				select case(trim(adjustl(ctl%str_source_den)))
				case("G(X)")
					ctl%source_fden=source_ana
				case("INDVD")
					ctl%source_fden=source_simu
					call tp%get_sub_para("--density bin criteria",str,ier)
					read(unit=str, fmt=*) ctl%den_bin_cri
				case default
					print*, "error! define density estimation method:", trim(adjustl(ctl%str_source_den))
					stop
				end select
			case("density integration method")
				read(unit=pvalue, fmt=*) ctl%str_density_est_method
				select case(trim(adjustl(ctl%str_density_est_method)))
				case("1D")
					ctl%fden_ana_est_method=fden_ana_est_method_1d_iso
				case("2D")
					ctl%fden_ana_est_method=fden_ana_est_method_2d
				case default
					print*, "error! define the fden ana est method=", &
						trim(adjustl(ctl%str_density_est_method))
					stop
				end select
			case("loss cone")
				read(unit=pvalue, fmt=*) ctl%include_loss_cone
				if(ctl%include_loss_cone.ge.1)then
					call tp%get_sub_para("--rtd  eta",str,ier)
					read(unit=str, fmt=*) ctl%eta_td_star
				end if
			case("clone scheme")
				read(unit=pvalue, fmt=*) ctl%clone_scheme
				if(ctl%clone_scheme.ge.1)then
					call tp%get_sub_para("--e0",str,ier)
					read(unit=str, fmt=*) clone_e0_factor_input
					call tp%get_sub_para("--clone separation in energy",str,ier)
					read(unit=str, fmt=*) clone_bd_sep
					log_clone_bd_sep=log(clone_bd_sep)
				end if
				ctl%del_cross_clone=ctl%clone_scheme
			case("del exit emin")
				read(unit=pvalue, fmt=*) ctl%del_exit_min 
			
			case("stellar evolution")
				read(unit=pvalue, fmt=*) ctl%include_stellar_evolution
				if(ctl%include_stellar_evolution.ge.1)then
					call tp%get_sub_para("--fraction of mass loss to gas reservoir",str,ier)
					read(unit=str,fmt=*) ctl%stellar_evolution_fraction_mass_loss_to_reservior
				end if
			case("chattery")
				read(unit=pvalue, fmt=*) ctl%chattery
			case("trace all samples")
				read(unit=pvalue, fmt=*) ctl%trace_all_sample

				if(ctl%trace_all_sample.ge.1)then
					call tp%get_sub_para("--track tidal disruption",str,ier)
					read(unit=str, fmt=*) ctl%output_track_td
					call tp%get_sub_para("--track plunge",str,ier)
					read(unit=str, fmt=*) ctl%output_track_plunge
					call tp%get_sub_para("--track emri",str,ier)
					read(unit=str, fmt=*) ctl%output_track_emri
				end if 
			case("radiative efficiency")
				read(unit=pvalue,fmt=*) ctl%radiative_efficiency
			case default
				print*, "error! para not defined:", trim(adjustl(tp%name))
				stop
			end select
		end associate
	end do
	spp_new%mbh=spp_new%mbh_dmless*m0_cl
	mbh_radius=spp_new%mbh/(my_unit_vel_c**2)

end subroutine