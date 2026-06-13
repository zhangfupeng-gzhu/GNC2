module model_basic
    use md_chain_pointer
    use md_bk_species
	use md_by_particle
	! use md_events
	use md_dms
	use constant
    use md_mobse_stellar_single
	use md_fc_type
	use md_iso
	use md_chain

	real(8):: log10rh, mbh_iso_radius, mbh_radius
	real(8):: m0_cl, r0_cl
	integer mbh_spin
    real(8) rhmin,rhmax, nsc_radius_eff
	real(8) sample_logemin, sample_logemax
	real(8) sample_logrmin, sample_logrmax
	real(8) sample_emin,sample_emax
	real(8) sample_nxgx_logemin,sample_nxgx_logemax
	real(8) energy_transition ! when phi_star(r_t)=Mbh/r_t=x_t
	real(8) radius_transition ! when phi_star(r_t)=Mbh/r_t
	real(8)::clone_bd_sep,log_clone_bd_sep
	integer::sample_threads_number
	type type_mdehnen
		real(8) mtot
		real(8) ra_crit
		real(8) gamma
	end type
	type type_mplummer
		real(8) mtot
		real(8) ra_crit
	end type
	! type(core_comp_type),pointer::cct_share
	type(s1d_type),pointer::fc_share
	type(s1d_type)::common_aux
	!type(s1d_ird_type)::common_jc
	type(s2d_type)::common_dee_log, common_djj_log, common_pd_log
	!type(s2d_type)::common_rp, common_ra
	!type(s1d_type)::common_gx
	!type(s1d_type)::common_barp
	!type(s1d_ird_type)::common_gx_ir
	type(s1d_ird_type),pointer::fc_ir_share
	type(s2d_type),pointer::gxj_share
	type(s2d_type)::gxjcr
    real(8) fgx_g0, delta_phi_correction 
    

	type control_type 
		type(s1d_type)::ini_nper_bin(20)
		type(s1d_type)::ini_nx_log(20)   ! in log bin of x
		type(s1d_type)::ini_nx_tot   ! in log bin of x
		type(s1d_type)::ini_ge(20)
		type(s1d_type)::ini_ge_tot
		type(s1d_type)::ini_frho(20)  !density of stellar components
		type(s1d_type)::ini_fna(20)  !cumulative number of stellar components
		type(s1d_type)::ini_fna_tot
		type(s1d_type)::ini_frho_tot  !density of stellar components
		type(s1d_type)::ini_fphi(20)
		type(s1d_type)::ini_fphi_tot 
		type(s1d_type)::ini_fma_tot
		type(type_mdehnen)::dehnen(20) 
		type(type_mplummer)::plummer(20) 
		real(8) rmin_factor, rmax_factor
		real(8) loge_max_factor, loge_min_factor! energy boundary to remove particle
		real(8) log10rmin_factor, log10rmax_factor
		real(8):: mass_ref
		real(8):: bin_mass_min, bin_mass_max
        real(8) n_basic
        !real(8) weight_asym
		real(8):: bin_mass(20)
		real(8):: bin_mass_m1(20), bin_mass_m2(20), asymptot_ini(n_tot_comp+1,20)
		real(8):: bin_total_number(20)
		real(8):: bin_mass_particle_number(20)  ! the number of particles at each bin
		real(8):: bin_fracmass(20)   ! the fraction of mass in this bin respect to the whole cluster
		real(8):: bin_fracmass_now(20)
		real(8):: asymptot_now(n_tot_comp+1,20)  ! current value of asymptot
		
		real(8):: num_particle_tot, mass_particle_tot
		real(8):: type_num_particle_tot(n_tot_comp)
		real(8):: type_mass_particle_tot(n_tot_comp)
		
		real(8):: ini_weight_n(20) 
		real(8):: Weight_n(20) 
        real(8):: bin_mass_N(20,4)=0
		real(8)  total_time
        real(8) burn_in_time, burn_in_time_input        
		real(8) sigma0_bh, energy0, energy_min, energy_max 
		real(8) run_snap_time_i, run_snap_time_f, run_snap_time_0
		real(8) run_snap_time_inner_steps_i, run_snap_time_inner_steps_f 
		real(8) total_energy,phi_energy,energy_lost_emin!, total_energy_new
		real(8) clone_e0, v0,  n0!, rh_min, rh_max
		real(8) rbd, gx_conv_cri, md_fine_rtd_range, md_alpha_cri
		real(8)::den_bin_cri=0.1
		real(8)::emri_cri=1d-2 ! should be less than 0.1
		real(8) x_boundary, energy_boundary!, !birth_position, x_birth
        real(8) ts_spshot_dt, ts_snap_dt_per_snap    !ts_spshot_tnr, the shapshottime in unit of 
                                                !two body relaxation timescale, in unit of Myr
		real(8) tmax_timestep 
		real(8)::mass_move_out_of_emin 
		real(8)::init_adb_mbh_log_factor=-4
		real(8)::trlx_rh0 
		real(8)::emax_boundary 
		real(8)::eta_td_star

		! real(8)::king_model_rho0, king_model_rt, king_model_r0, king_model_sigma
		! real(8)::king_model_w, king_model_rho1,king_model_phirt
		
		real(8)::init_adb_mbh_log_factor_target

		real(8)::metal_z
		real(8)::stellar_evolution_fraction_mass_loss_to_reservior 
		real(8)::radiative_efficiency

		character*(3) str_jbin_bd, str_fj_bd, str_ebin_type
		character*(3) str_task_mode, str_snap_mode,str_byctype
		character*(6) str_dejmodel,time_unit, str_model_intej
		character*(6) str_end_time_mode
        character*(4)   str_method_interpolate
        character*(200) default_para_file_dir 
		character*(200) rv_nw_conv_dir
		character*(8) str_ini_den_model(20),  str_dc_grid_type
		character*(8) str_barge_grid_type,str_barge_evl_method,str_adb_est_method
		character*(8) str_dc_method, str_source_den,str_density_est_method 
		
		integer:: timestep_mode
		integer:: update_2bdrlx_per_shapshot

		integer:: bin_mass_simulation_particle_number_tot(20)
		integer:: bin_mass_simulation_particle_number_comp_tot(20,n_tot_comp_sg)

		integer:: clone_factor(20)
		integer:: idx_stellar_type(20)  ! the list of types 
		integer:: idx_stellar_type_sg(20), idx_stellar_type_by(20)
		integer:: exist_stellar_type(20) ! whether the specific stellar type exists
		integer,allocatable:: ini_stellar_tot(:) ! tot number of stellar object of each type
		integer,allocatable:: ini_stellar_each_mass(:,:) ! tot number of stellar object of each type and mass
		integer::ini_adb_increase_nt=20 
		integer fden_ana_est_method
		integer random_seed_size
		integer,allocatable::random_seed(:)
		
		integer:: min_sample_in_mass_bin=40
		
		!integer num_mdehnen
		!integer num_mking
		!integer num_
		
		integer::current_version_number 
		integer run_seq_idx
		integer output_dms_freq
		integer output_data_freq
        integer intaskmode, insnapmode 
		integer ini_model_list(20) 
		integer ini_model_list_in(20)
		integer ini_den_model(20) ! the type number of the model of mass bin i
		!integer ini_den_model_idx(10) ! the index of the model in the same type of mass bin i
		!integer ini_den_model_idx_in(10)
		integer num_of_ini_den_model
		integer num_mdehnen 
		integer num_mspow 
		integer num_mplummer 
		integer time_run_mode
		integer enable_evl_mbh
		 
        integer num_clone_created, num_clone_elim, num_boundary_created
		integer num_boundary_elim  
        integer boundary_method, boundary_fj 
		integer num_step_per_update
		integer n_spshot, n_spshot_bg, n_spshot_total
		integer Dejmodel, include_loss_cone 
		integer num_get_sample_para_exact, num_get_sample_para_kpl, num_of_loops,num_of_lc
		integer model_intej
		integer i_spshot_bg_last 
		integer::binary_3body_encounter, include_bhb_mbh_3body 
		integer::same_rseed_evl, same_rseed_ini	
		integer::gw_radiation_otby 
		integer::output_gw_emri_data_bin 
		integer::clone_scheme, del_exit_min
		integer::trace_all_sample  
		integer::all_exact
		integer::include_stellar_evolution
		integer::init_adb_mbh_inc 
		integer::two_body_relaxation_on

		integer stellar_collision_pair_mode
		integer stellar_collision_method
		integer collision_consider_weight
		
		integer::del_cross_clone
		integer::consider_by_types(5)
        integer::num_cby_types
		integer::use_tidal_spin
		integer::collection_data_method
        integer burn_in_snap
		integer jbin_type,ebin_type
		integer idxstar, idxsbh, idxbbh, idxmsb, idxns, idxwd,idxbd,idxrg, idxdm, idxnHe 
		integer barge_grid_type
		integer dc_grid_type
		integer barge_evl_method
		integer adb_est_method
		integer ini_mass_bin_mode 
		integer method_interpolate
		!integer::npar_sam_tot
		integer::m_bins ! number of mass bin
		integer::diff_coeff_bins ! for diffusion coefficients
		integer::dstr_bins_r       ! for other 1D distribution functions of r
		integer::dstr_bins_e       ! for other 1D distribution functions of e 
		integer::dstr_bins_j ! for gxj_ir
        integer::idx_ref
        !integer::consider_merge_kicks
		integer:: debug=0
        integer:: ini_sample_sg_mode  
        integer:: bin_mass_emax_out(20,4)=0
		integer:: chattery ! <=2, normal; =3, add tdial details
        integer:: ntasks, ntask_total 
        integer:: seed_value, ntask_bg
		integer::nblock_mpi_bg, nblock_mpi_ed, nblock_size
		integer::output_track_td, output_track_emri
		integer::output_track_plunge, output_track_norm
		integer:: source_fden
		integer::get_dc_method
		integer::n_tot_samples		
		logical replace_sample_eceed_emax
	end type
	!type(ini_par_type)::ipt

	type(control_type),target::ctl 
    type(chain_type)::Allsams

	type(particle_samples_arr_type)::bksams_arr_ini
	type(chain_type):: bkstars,   bksbhs
	type(chain_type):: bksams 
	type(chain_type):: bksams_norm 

	type(particle_samples_arr_type):: bkstars_arr 
	type(particle_samples_arr_type):: bksbhs_arr
	type(particle_samples_arr_type):: bksbhs_arr_norm 
	type(particle_samples_arr_type),target:: bksams_arr_norm   ! only normal samples
	! type(particle_samples_arr_type),target:: bksams_arr_emin   ! only exit to emin samples
	type(particle_samples_arr_type),target:: bksams_arr_norm_sbh
	type(particle_samples_arr_type):: bksams_arr    ! include all samples
	type(particle_samples_arr_type):: bksams_arr_merge  ! due to gw capture   
	type(samples_type_pointer_arr),target::bksams_pointer_arr
	type(diffuse_mspec),target::dms 
	
    integer,parameter:: chattery_out_unit_0=1383829393
    integer chattery_out_unit  
    real(8),parameter::my_unit_vel_c5=my_unit_vel_c**5
	real(8)::dc_grid_xstep, dc_grid_ystep
	 
	integer::aux_function_bin_size=200
	integer::f12_function_bin_size_1=100
	integer::f12_function_bin_size_2=30
	real(8)::star_dc_int_acc_r=1d-9
	real(8)::star_dc_int_acc_a=1d-12
	real(8)::fr_funcs_int_acc_r=1d-25
	real(8)::fr_funcs_int_acc_a=1d-15
	real(8)::pd_int_acc_a=1d-12
	real(8)::pd_int_acc_r=1d-8

	integer::max_self_con_iter=15

	real(8)::gx_func_max_step=10
	real(8)::gx_func_min_step=0.05  
	real(8) clone_e0_factor     
	real(8) clone_e0_factor_input 
	real(8):: clone_emax
	real(8)::jmin_value, jmax_value
	real(8)::log10jmin_value,log10jmax_value
	
	real(8)::emin_value, emax_value 
	real(8) nx_exmax, nx_logemin,nx_logemax 

	integer,parameter:: fden_ana_est_method_1d_iso=1, fden_ana_est_method_2d=2
	integer,parameter:: timestep_mode_tnr=1, timestep_mode_trh=2 

	integer,parameter::ini_mass_bin_mode_given=1, ini_mass_bin_mode_kroupa=2, ini_mass_bin_mode_pow=3
	integer,parameter::ini_mass_bin_mode_topheavy=4

	integer,parameter:: get_dc_method_fast=1,get_dc_method_general=2 
  
	
	integer,parameter::source_simu=2, source_ana=1 
	integer,parameter::barge_grid_type_iregular_phi=3!, barge_grid_type_iregular_z=4
	integer,parameter::barge_grid_type_iregular_jc=5
	integer,parameter::ini_den_model_dehnen=1, ini_den_model_plummer=2
	integer,parameter:: barge_evl_method_direct=1
	integer,parameter:: barge_evl_method_grid_2d=3  
	integer,parameter:: dc_grid_irregular=2 
	integer,parameter::adb_est_method_fast=2, adb_est_method_acc=1 

	integer,parameter::method_int_nearst=1, method_int_linear=2
	integer,parameter::time_run_mode_snap=1,time_run_mode_ttot=2 
	integer,parameter::boundary_fj_iso=1 
	integer,parameter::ini_sample_mode_given=1, ini_sample_mode_mobse=2 
	integer,parameter::ini_sample_mode_bse=3
	integer,parameter::dejmodel_EJ=1 
	integer,parameter::task_mode_new=1, task_mode_append=2
	integer,parameter::snap_mode_new=1, snap_mode_append=2, snap_mode_one=3   
     
	integer,parameter::MAX_LENGTH=100000
!    integer,parameter::default_track_size=10000
    integer(8),parameter::MAX_RUN_LENGTH=int(1d7) 
	integer,parameter::n_orb_track_block=100000
    integer,parameter::n_orb_track_block_max=n_orb_track_block*16 
    integer nsize_chain_bk, nsize_chain_by, nsize_arr_bk, nsize_arr_by, &
    nsize_arr_bk_norm, nsize_arr_by_norm, nsize_arr_bk_pointer, &
    nsize_arr_by_pointer, nsize_tot_bk, nsize_tot_by, nsize_tot 
	integer::boundary_sts_emax_cros=0
	integer::update_correction_emax=0, running_correction_emax=0, self_correction_emax=0
	!logical::correction_emax_at_current_snap=.false.

	real(8),parameter::my_unit_vel_c3=my_unit_vel_c**3
	real(8),parameter::my_unit_vel_c2=my_unit_vel_c*my_unit_vel_c

	real(8) imf_para_nt(4), imf_para_nc(3), imf_para_nq
	real(8) imf_para_mt(4), imf_para_mc(3), imf_para_mq
	real(8),parameter::topheavy_alpha(3)=(/-0.3d0,-1.3d0,-1.6d0/)
	real(8),parameter::topheavy_xb(4)=(/0.01d0,0.08d0,0.5d0,150d0/)

! #if defined(MOBSE)
	! type(obj_massbin)::omsamples
! #endif	
	type(ISO_table_type)::iso_kerr!, iso_schw 
contains

	subroutine all_chain_to_arr_single(sps, sps_arr)
		implicit none
		type(chain_type),intent(in)::sps
		type(particle_samples_arr_type),intent(out)::sps_arr
	!	type(chain_type)::sps_tmp
		type(chain_pointer_type),pointer::sp
		integer i,n
		integer::nr
	!	print*, associated(sps_tmp%sp)
	!	call copy_sams(sps,sps_tmp)

		call sps%get_length(nr,type=1)

        !print*, "nr=", nr
        if(nr.eq.0)return
!		sp=>sps_tmp%sp(1)
!		print*, sps_tmp%sp(1)%ob%ac

		call sps_arr%init(nr)
		
        sp=>sps%head
		call sp%chain_to_arr_single(sps_arr%sp(1:nr),nr)

	end subroutine
	
	subroutine get_exit_flag_str(exit_flag, str_flag)
		implicit none
		integer exit_flag
		character*(100) str_flag 
		select case(exit_flag)
		
		case (exit_normal)
			str_flag="NORMAL"
		case (exit_tidal)
			str_flag="TIDAL"
		case (exit_max_reach)
			str_flag="NMAX"
		case (exit_boundary_min)
			str_flag="BOUNDAY_MIN"
		case(exit_boundary_max)
			str_flag="BOUNDARY MAX" 
		case(exit_invtransit)
			str_flag="INVERSE TRANSIT"
		case(exit_tidal_empty)
			str_flag="TD EMPTY"
		case(exit_tidal_full)
			str_flag="TD FULL"        
		case default
			str_flag="Null"
		end select
	end subroutine 
	
end module


module md_star_pot
	use com_sts_type
	implicit none
	type star_pot_para
		real(8) M_r_within_max, spt_rho_rmin, phi_r1r2_s, phi_r1r2_s2, m_r_within_min
		real(8) N_r_within_max, phi_star0
		!real(8) emin_factor, emax_factor
		!real(8) sample_logemax, sample_logemin
		logical::has_set_density=.false.
		logical::has_init_tables=.false.
		logical::has_set_rhomin=.false.
		real(8) mbh, mbh_dmless
		!logical::has_set_beta=.false.
		type(s1d_type)::fphi_star,fma_star,frho_star!, beta
	contains
		procedure::init=>init_spp_tables
	end type
	type,extends(star_pot_para)::star_pot_para_record
		type(s1d_type)::fgx_ir
	end type 
	logical::has_made_trade_off
	real(8),parameter::weight_correct_factor=3d0
	real(8)::current_weight_correct_factor=weight_correct_factor
	type(star_pot_para)::spp_new, spp_old
	type(star_pot_para_record)::spp_record(20)
	integer::n_record=0
	!real(8) star_energy_emax_limit
	private::init_spp_tables
contains 
	subroutine copy_spp_to_record(spp,spp_record)
		implicit none
		class(star_pot_para_record)::spp_record
		class(star_pot_para)::spp
		spp_record%fphi_star=spp%fphi_star
		spp_record%fma_star=spp%fma_star
		spp_record%frho_star=spp%frho_star
		spp_record%M_r_within_max=spp%M_r_within_max
		spp_record%mbh=spp%mbh
		spp_record%mbh_dmless=spp%mbh_dmless
	end subroutine
	subroutine init_spp_tables(spp,logrmin,logrmax,dstr_bins)
		implicit none
		class(star_pot_para)::spp
		real(8) logrmin,logrmax
		integer dstr_bins
		call spp%fphi_star%init(logrmin,logrmax,dstr_bins,sts_type_dstr)
		call spp%fphi_star%set_range()
		call spp%fma_star%init(logrmin,logrmax,dstr_bins,sts_type_dstr)
		call spp%fma_star%set_range()
		call spp%frho_star%init(logrmin,logrmax,dstr_bins,sts_type_dstr)
		call spp%frho_star%set_range()
		!call spp%beta%init(logrmin,logrmax,dstr_bins,sts_type_dstr)
		!call spp%beta%set_range()

		spp%has_init_tables=.true.
		spp%has_set_density=.false.
		!spp%has_set_rhomin=.false.
	end subroutine
	subroutine save_spp_tables_bin(spp,funit)
		implicit none
		type(star_pot_para)::spp
		integer funit
		integer mbh_evl
		write(funit) spp%fphi_star,spp%fma_star,spp%frho_star
		write(funit) spp%M_r_within_max, spp%spt_rho_rmin, spp%phi_r1r2_s, &
			spp%phi_r1r2_s2, spp%m_r_within_min,spp%has_set_density,spp%has_init_tables, &
			spp%has_set_rhomin, spp%N_r_within_max!, spp%spt_rho_rmax
		!write(funit) spp%emin_factor, spp%emax_factor, spp%sample_logemax, spp%sample_logemin
		!if(mbh_evl.ge.1)then
		write(funit) spp%mbh, spp%mbh_dmless
		!end if
	end subroutine
	subroutine read_spp_tables_bin(spp,funit)
		implicit none
		type(star_pot_para)::spp
		integer funit
		integer mbh_evl
		read(funit) spp%fphi_star,spp%fma_star,spp%frho_star 
		read(funit) spp%M_r_within_max, spp%spt_rho_rmin, spp%phi_r1r2_s, &
			spp%phi_r1r2_s2, spp%m_r_within_min,spp%has_set_density,spp%has_init_tables, &
			spp%has_set_rhomin, spp%N_r_within_max!, spp%spt_rho_rmax
		!read(funit) spp%emin_factor, spp%emax_factor, spp%sample_logemax, spp%sample_logemin
		!spp%spt_rho_rmax=0
		read(funit) spp%mbh, spp%mbh_dmless
		!end if
	end subroutine
	subroutine check_spp_data(spp,ier)
		implicit none
		type(star_pot_para)::spp
		integer ier
		ier=0
		if(.not.spp%has_init_tables)then
			print*, "error! spp has not init tables"
			ier=-1
		end if
		if(.not.spp%has_set_density)then
			print*, "error! spp%frho_star has not been set"
			ier=-2
		end if
		if(.not.spp%has_set_rhomin)then
			print*, "error! spprhomin has not been set"
			ier=-3
		end if
	end subroutine
end module


 
