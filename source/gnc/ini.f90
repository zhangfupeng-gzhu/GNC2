 
subroutine get_rh_in_pc(m0, r0)
	implicit none
	real(8) r0,m0
	r0=3.1*(m0/4e6)**0.55
end subroutine
   
subroutine end()
	use com_main_gw
	print*, "stop mpi"
	call stop_mpi()
	print*, "deallocate"
	call deallocate_chains_arrs()
end subroutine

subroutine head(fmodel,fmbin)
	use com_main_gw
	use model_config
	implicit none
	character*(*) fmodel, fmbin
	! print*, "start"
	call readin_model_par(trim(adjustl(fmodel)))
	call init_mpi()
	!print*, "apply paras"
	call apply_paras()
	call readin_mass_bins(trim(adjustl(fmbin)))
	!print*, "print models"
	call print_models()	
	if(rid.eq.0)then
	 	call pa_now_used%print("Used Parameters")
		if(ctl%intaskmode.ne.task_mode_append.and.ctl%insnapmode.ne.snap_mode_append)then
			call pa_now_used%save_to_txt("./output/used_parameters.txt")
		else
			call pa_now_used%save_to_txt("./output/used_parameters_append.txt")
		end if
	end if
	call init_model()
	!!stop

	if(rid.eq.0)then
	    !call print_model_par()
        call print_current_code_version()
    end if
    !rid=2
	if(rid.eq.0)then
		write(*,*) 'root process id = ', proc_id
	end if
	if(rid.eq.1)then
		write(*,*) 'cld process id = ', proc_id
	end if

	if(ctl%include_loss_cone.ge.1)then
		call input_rp_iso(iso_kerr, ctl%rv_nw_conv_dir)
	endif
        
	mpi_master_id=0

end subroutine
subroutine print_models()
	use com_main_gw
	implicit none
	integer i
	if(rid.eq.0)then
		do i=1, ctl%m_bins
			select case(trim(adjustl(ctl%str_ini_den_model(i))))
			case("Dehnen")
				print*, "i,model,mtot,ra,gamma=", i, "Dehnen", ctl%dehnen(i)%mtot, ctl%dehnen(i)%ra_crit, &
					ctl%dehnen(i)%gamma
			case("Plummer")
				print*, "i,model,mtot,ra=", i, "Plummer", ctl%dehnen(i)%mtot, ctl%dehnen(i)%ra_crit
			end select
		end do
	end if
end subroutine

subroutine get_rh_vh_nh(m0, r0_cl, nh, vh)
	use constant
	implicit none
	real(8),intent(in):: m0, r0_cl
	real(8),intent(out):: nh, vh
	real(8),parameter::rgc=8.32d3
	!real(8),parameter::r0=3.1
	real(8),parameter::n0=2d4

	!r0_cl=r0*(Mbh/4d6)**0.55*pc	
	vh=sqrt(m0/r0_cl)    
	nh=m0/r0_cl**3
	

end subroutine
subroutine init_model_ctl()
	use com_main_gw
	use md_mbh_evl_acc
	!use md_dms_saving_data
	implicit none
	integer i
	integer(8) num_int
    real(8) rh
 
    id_saver=0 
	emin_value=(1-jmax_value**2)**0.5
	emax_value=(1-jmin_value**2)**0.5
	
	call get_rh_vh_nh(m0_cl, r0_cl, ctl%n0, ctl%v0)
    log10rh=log10(r0_cl)
    rhmin=r0_cl/emax_factor/2d0; 
    rhmax=r0_cl/emin_factor/2d0
    !===================test============= 
	call get_rh_in_pc(spp_new%mbh,rh)
	ctl%sigma0_bh=sqrt(spp_new%mbh/(rh*pc))
	ctl%energy0=-m0_cl/r0_cl 
	ctl%rbd=r0_cl*0.5d0/ctl%x_boundary
	if(ctl%clone_scheme.ge.1)then 
		clone_e0_factor=clone_e0_factor_input
		ctl%clone_e0=clone_e0_factor*ctl%energy0 
	end if 
	
    ctl%energy_boundary=ctl%x_boundary*ctl%energy0
      
    ctl%n_basic=minval(ctl%ini_weight_n(1:ctl%m_bins))

	do i=1, ctl%m_bins
		num_int=ctl%ini_weight_n(i)/ctl%n_basic
		if(abs(num_int-(ctl%ini_weight_n(i)/ctl%n_basic))>1d-1)then
			print*, i, int(ctl%ini_weight_n(i)/ctl%n_basic),(ctl%ini_weight_n(i)/ctl%n_basic)
			print*, "error! set particle weighting an integer multiple of the minimium one"
			stop
		end if
	end do
    ctl%weight_n(1:ctl%m_bins)=int(ctl%ini_weight_n(1:ctl%m_bins)/ctl%n_basic)
 
    ctl%num_boundary_created=0
	ctl%num_boundary_elim=0
    ctl%num_clone_created=0  
	ctl%nblock_size=int(ctl%diff_coeff_bins/ctl%ntasks)
    ctl%nblock_mpi_bg=rid+1
    ctl%nblock_mpi_ed=ctl%diff_coeff_bins


	!call aux_for_de%init(0d0,1d0,ctl%diff_coeff_bins,sts_type_grid)
	!call aux_for_de%set_range()
	call common_aux%init(0d0, pi/2d0, aux_function_bin_size, sts_type_grid)
	do i=1, common_aux%nbin
		common_aux%xb(i)=pi/2d0*((i-1)/dble(aux_function_bin_size-1))**2
	end do
	if(allocated(ctl%ini_stellar_tot))then
		deallocate(ctl%ini_stellar_tot ,ctl%ini_stellar_each_mass)
	end if
	allocate(ctl%ini_stellar_tot(n_tot_comp),ctl%ini_stellar_each_mass(n_tot_comp,ctl%m_bins))

	call print_current_code_version()
	
	if(ctl%enable_evl_mbh.ge.1)then
		call mbh_mmg%init()
	end if

	ctl%idx_stellar_type_sg(1:n_tot_comp_sg)=(/star_type_ms,star_type_bh, &
	star_type_ns, star_type_wd, star_type_bd, star_type_rg,star_type_dark_matter,star_type_nakedHe/) 
	
	ctl%idx_stellar_type(1:n_tot_comp_sg)=ctl%idx_stellar_type_sg(1:n_tot_comp_sg) 
	
end subroutine

subroutine init_model()
	!use com_main_gw
	implicit none
	call init_model_ctl()
    call init_chattery()
	! call init_pro()
end subroutine
subroutine init_chattery()
    use com_main_gw
    implicit none
    character*(5) tmprid
    logical,save::first=.true.
    if(ctl%chattery.ge.2)then
        if(first)then
            first=.false.
        else
            close(unit=chattery_out_unit)
        end if
        chattery_out_unit=chattery_out_unit_0+rid
        write(unit=tmprid,fmt="(I5)") rid
        open(unit=chattery_out_unit,file="output/chattery_"//trim(adjustl(tmprid)), &
            status='replace')

    end if
    chattery_out_unit=6
end subroutine

subroutine set_seed(same_seed,seed_value)
	use com_main_gw
	implicit none
	logical same_seed
	integer seed_value
	if(same_seed)then
	    call same_random_seed(seed_value)
	else
		call random_seed()
	!   call same_random_seed(seed_value+rid)
	end if
end subroutine  


