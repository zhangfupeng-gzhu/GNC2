#ifdef HDF5
 
 
SUBROUTINE save_dms_aux_data_hdf5(dm,spp,group_id)
	use md_hdf5	
	use com_main_gw
	use md_event_datas
	implicit none
	type(diffuse_mspec)::dm
	type(star_pot_para)::spp
	integer i,error
	character*(4) tmpi
	integer(hid_t)::group_id,sub_group_id 
	call h5gcreate_f(group_id,"aux",sub_group_id,error)
	call dms%jc%save_hdf5(sub_group_id,  "jc_dmless")
	call spp%fphi_star%save_hdf5(sub_group_id,  "phi_star")
	call dms%alpha_r%save_hdf5(sub_group_id,  "alpha") 
	call dms%rc%save_hdf5(sub_group_id,  "rc")
	call spp%fma_star%save_hdf5(sub_group_id,"fma_star")
	call dms%fr_phi%save_hdf5(sub_group_id,  "rmax")  
	call dms%barp_ir%save_hdf5(sub_group_id,"barp_ir") 
	call dm%surface_den%save_hdf5(sub_group_id,"surface_den")
	call dm%cum_sur_den%save_hdf5(sub_group_id,"cum_sur_den") 
	call output_jlc(sub_group_id)
    call dm%pd%save_hdf5(sub_group_id,"pd")
    call dm%rp%save_hdf5(sub_group_id,"rp")
	call dm%ra%save_hdf5(sub_group_id,"ra")
 
	call h5gclose_f(sub_group_id,error)
	!call hdf5_file%close()
end subroutine
subroutine get_jlc(rp_dmless,jlc)
	use com_main_gw
	implicit none
	real(8) rp_dmless, ex, jc, phirp
	type(s1d_type)::jlc
	integer j

	jlc=dms%jc
	call get_phi_star_full_range(spp_new,log10(rp_dmless),phirp)
	do j=1, jlc%nbin
		ex=10**jlc%xb(j) 
		jc=dms%jc%fx(j)
		if(10**phirp+spp_new%mbh_dmless/rp_dmless-ex>0.and.rp_dmless<dms%rc%fx(j))then
			jlc%fx(j)=(2*(10**phirp+spp_new%mbh_dmless/rp_dmless-ex))**0.5*rp_dmless/jc
		else
			jlc%fx(j)=1d0
		end if 
	end do
end subroutine
subroutine output_jlc(group_id)
	use md_hdf5
	use com_main_gw
	implicit none
	integer(hid_t)::group_id
	real(8) jc, phirp,rp,ex
	integer i,j,stellar_type
	type(s1d_type)::jlc
	integer have_outputed(2)
	!jlc=dms%jc
	have_outputed=0
	do i=1, n_tot_comp_sg
		call get_type_from_ctl_obidx_sg(i,stellar_type)	
		if(ctl%exist_stellar_type(i)>0)then
			select case(stellar_type)
			case(star_type_ms,star_type_bd,star_type_rg)
				if(have_outputed(1).eq.0)then
					have_outputed(1)=1
				else
					cycle
				end if
				rp=rd_sun*(spp_new%mbh)**0.333333/r0_cl
				call get_jlc(rp,jlc)				
				call jlc%save_hdf5(group_id,"jlc_star")
			case(star_type_bh,star_type_wd,star_type_ns,star_type_nakedHe)
				if(have_outputed(2).eq.0)then
					have_outputed(2)=1
				else
					cycle
				end if
				rp=8*spp_new%mbh/1e8/r0_cl
				call get_jlc(rp,jlc)
				call jlc%save_hdf5(group_id,"jlc_compact")
			end select
		end if
	end do
end subroutine

subroutine output_ctl_ini_dstr_hdf5(fl)
	use md_hdf5	
	use com_main_gw
	implicit none
	type(diffuse_mspec)::dm
	character*(*) fl
	integer i,j, error
	character*(4) tmpi
	type(s1d_type)::fe,rc,jc, fpot,fma,fm_dstr
	type(hdf5_file_type)::hdf5_file
	type(hdf5_table_type)::hdf5_table
	INTEGER(HID_T) :: group_id,sub_group_id       ! File identifier	
 
	call hdf5_file%open(trim(adjustl(fl))//".hdf5")
	call ctl%ini_nx_tot%save_hdf5(hdf5_file%file_id,  "nx_tot")
	call ctl%ini_fma_tot%save_hdf5(hdf5_file%file_id,  "fma_tot")
	call ctl%ini_fna_tot%save_hdf5(hdf5_file%file_id,  "fna_tot")
	call h5gcreate_f(hdf5_file%file_id,"nx",sub_group_id,error)
	do i=1, ctl%m_bins
		write(unit=tmpi, fmt="(I4)") i
		call ctl%ini_nx_log(i)%save_hdf5(sub_group_id, trim(adjustl(tmpi))//"_ini_nx")
		call ctl%ini_nper_bin(i)%save_hdf5(sub_group_id, trim(adjustl(tmpi))//"_ini_nper_bin")
	end do

	call ctl%ini_frho_tot%save_hdf5(hdf5_file%file_id, "ini_frho")
	call h5gcreate_f(hdf5_file%file_id,"frho",sub_group_id,error)
	do i=1, ctl%m_bins
		write(unit=tmpi, fmt="(I4)") i
		call ctl%ini_frho(i)%save_hdf5(sub_group_id, trim(adjustl(tmpi))//"_ini_frho")
	end do
	call h5gclose_f(sub_group_id,error)
	call ctl%ini_fphi_tot%save_hdf5(hdf5_file%file_id,"ini_fphi")
	do i=1, ctl%num_mdehnen
		if(spp_new%mbh.eq.0d0.and.(ctl%dehnen(i)%gamma.eq.1d0.or.ctl%dehnen(i)%gamma.eq.0d0))then
			fe=dms%all%all%barge
			call get_dehnen_fe(fe, ctl%dehnen(i)%gamma, ctl%dehnen(i)%ra_crit, ctl%dehnen(i)%mtot)
			call fe%save_hdf5(hdf5_file%file_id,"fe_theory")
			rc=dms%rc
			call get_dehnen_rc(rc, ctl%dehnen(i)%gamma, ctl%dehnen(i)%ra_crit, ctl%dehnen(i)%mtot)
			call rc%save_hdf5(hdf5_file%file_id,"rc_theory")
			jc=dms%jc
			call get_dehnen_jc(jc,rc, ctl%dehnen(i)%gamma, ctl%dehnen(i)%ra_crit,&
				ctl%dehnen(i)%mtot)
			call jc%save_hdf5(hdf5_file%file_id,"jc_theory")
		end if
	end do
	do i=1, ctl%num_mplummer
		if(spp_new%mbh.eq.0d0.and.(ctl%dehnen(i)%gamma.eq.1d0.or.ctl%dehnen(i)%gamma.eq.0d0))then
			fe=dms%all%all%barge
			call get_plummer_fe(fe,ctl%plummer(i)%ra_crit, ctl%plummer(i)%mtot)
			call fe%save_hdf5(hdf5_file%file_id,"fe_theory")
			fpot=spp_new%fphi_star
			do j=1, fpot%nbin
				call get_plummer_pot(fpot%fx(j),10**fpot%xb(j),ctl%plummer(i)%mtot,ctl%plummer(i)%ra_crit)
			end do
			call fpot%save_hdf5(hdf5_file%file_id,"fpot_theory")
			fma=spp_new%fphi_star
			call get_plummer_fma_s1d(fma,ctl%plummer(i)%ra_crit, ctl%plummer(i)%mtot)
			call fma%save_hdf5(hdf5_file%file_id,"fma_theory") 
		end if
	end do
  
	do i=1, ctl%m_bins
		write(unit=tmpi,fmt="(I4)") i
		call ctl%ini_ge(i)%save_hdf5(hdf5_file%file_id,  "ge"//trim(adjustl(tmpi)))
		!call ctl%ini_ge(i)%print("ini ge")
	end do
	if(ctl%include_stellar_evolution.ge.1.and.ctl%ini_mass_bin_mode.eq.ini_mass_bin_mode_kroupa)then
		call fm_dstr%init(dms%all%all%fm_dstr%xmin, dms%all%all%fm_dstr%xmax, dms%all%all%fm_dstr%nbin,sts_type_dstr)
		fm_dstr%xb=dms%all%all%fm_dstr%xb
		call get_kroupa_imf(fm_dstr)
		call fm_dstr%save_hdf5(hdf5_file%file_id,"fm_kroupa")
	end if
	call hdf5_file%close()
end subroutine
subroutine get_kroupa_imf(fm)
	use com_main_gw
	implicit none
	type(s1d_type)::fm
	integer i
	real(8) m,fx,canonical_IMF_func
	do i=1, fm%nbin
		m=10**fm%xb(i)
		fx=canonical_IMF_func(imf_para_nt,imf_para_nc,imf_para_nq,m)
		fm%fx(i)=fx
	end do
end subroutine
subroutine write_down_attributions(dm,group_id)
	use md_hdf5
	use com_main_gw
	use md_dms_saving_data
	use md_mbh_evl_acc
	integer(HID_T)::group_id,attr_id  
	type(diffuse_mspec)::dm
	real(8) reff_now, rh_now
	type(hdf5_group_type)::hg

	call add_attr_dble(group_id,attr_id, "Time(Myr)", ctl%run_snap_time_f)
	call add_attr_dble(group_id,attr_id, "dT(Myr)", ctl%run_snap_time_f-ctl%run_snap_time_i)
	call add_attr_dble(group_id,attr_id, "Trlx_rh(Myr)", ctl%trlx_rh0)
	call add_attr_dble(group_id,attr_id, "MBH", spp_new%mbh)
	call add_attr_dble(group_id,attr_id, "MBH_dm", spp_new%mbh_dmless)
	call add_attr_dble(group_id,attr_id, "Mcluster(Msun)", spp_new%M_r_within_max*m0_cl)
	call add_attr_dble(group_id,attr_id, "Ncluster", spp_new%N_r_within_max*m0_cl)
	call add_attr_dble(group_id,attr_id, "r0(pc)", r0_cl/pc)
	call add_attr_dble(group_id,attr_id, "m0(msun)", m0_cl)
	call add_attr_dble(group_id,attr_id, "logemin", log10emin_factor)
	call add_attr_dble(group_id,attr_id, "logemax", log10emax_factor)
	call add_attr_dble(group_id,attr_id, "dmsemin", dms%emin)
	call add_attr_dble(group_id,attr_id, "dmsemax", dms%emax)
	call add_attr_dble(group_id,attr_id, "sampleemin", sample_logemin)
	call add_attr_dble(group_id,attr_id, "sampleemax", sample_logemax)
	if(spp_new%mbh_dmless.eq.0)then
		call add_attr_int(group_id,attr_id,"N_update_emax_cor", update_correction_emax)
		call add_attr_int(group_id,attr_id,"N_running_emax_cor", running_correction_emax)
		call add_attr_int(group_id,attr_id,"N_self_emax_cor", self_correction_emax)
	end if

	call add_attr_dble(group_id,attr_id, "reff(pc)", 10**nsc_radius_eff*r0_cl/pc)

	call get_rh_now(rh_now)
	call add_attr_dble(group_id,attr_id, "rh(pc)", 10**rh_now*r0_cl/pc)
	if(ctl%enable_evl_mbh.ge.1)then
		call hg%create(group_id,"dMbh")
		call write_hdf5_mass_mbh_growth(hg%group_id,mbh_mmg)
		call hg%close()
	end if
	
	call write_down_current_asymptot(group_id,"asymptot")
	call add_attr_dble(group_id,attr_id, "num_particle_tot", ctl%num_particle_tot)
	call add_attr_dble(group_id,attr_id, "mass_particle_tot(msun)", ctl%mass_particle_tot)
	call add_attr_dble_arr(group_id,attr_id, "num_type_particle", ctl%type_num_particle_tot,n_tot_comp_sg)
	call add_attr_dble_arr(group_id,attr_id, "mass_type_particle(msun)", ctl%type_mass_particle_tot,n_tot_comp_sg)
	call add_attr_dble(group_id,attr_id,"mass out of emin", ctl%mass_move_out_of_emin)
	call add_attr_int_arr(group_id, attr_id, "num_type_particle(simu)",ctl%bin_mass_simulation_particle_number_tot,n_tot_comp_sg) 

end subroutine
subroutine write_down_current_asymptot(group_id,tablename)
	use md_hdf5
	use model_basic
	implicit none
	integer(HID_T)::group_id
	type(hdf5_table_type)::h5table
	character*(*) tablename
	integer,parameter::nfields=n_tot_comp_sg+3+4
	real(8) tp(dms%n)
	integer i
	call h5table%init_table(nfields,dms%n,tablename)
	h5table%field_names=(/"      m1","      m2","    mbin","   nftot", "   mftot", &
			   "    ntot", "    mtot", &
	           "    star", "     sbh","      ns", "      wd", &
                  "      bd","      rg", "      DM", &
                  "     nHe"/)	
	h5table%field_types(1:nfields)=H5T_NATIVE_DOUBLE
	call h5table%prepare_write_table(group_id)

	call h5table%write_column_real(ctl%bin_mass_m1(1:dms%n))
	call h5table%write_column_real(ctl%bin_mass_m2(1:dms%n))
	call h5table%write_column_real(ctl%bin_mass(1:dms%n))
	tp(1:dms%n)=ctl%asymptot_now(1,1:dms%n)
	call h5table%write_column_real(tp)
	call h5table%write_column_real(ctl%bin_fracmass_now(1:dms%n))
	call h5table%write_column_real(tp*ctl%num_particle_tot)
	call h5table%write_column_real(ctl%bin_fracmass_now(1:dms%n)*ctl%mass_particle_tot)
	
	do i=2, n_tot_comp_sg+1
		tp(1:dms%n)=ctl%asymptot_now(i,1:dms%n)
		call h5table%write_column_real(tp)
	end do

end subroutine
subroutine get_reff_now(reff_now)
	use com_main_gw
	implicit none
	real(8) reff_now
	type(s1d_type)::s1d_tmp
	s1d_tmp=dms%cum_sur_den
	s1d_tmp%xb=dms%cum_sur_den%fx
	s1d_tmp%fx=dms%cum_sur_den%xb
	!s1d_tmp%xmin=minval(dm%cum_suf_den%fx)
	s1d_tmp%xmax=maxval(dms%cum_sur_den%fx)
	call s1d_tmp%get_value_l(0.5*s1d_tmp%xmax,reff_now)

end subroutine
subroutine get_rh_now(rh_now)
	use com_main_gw
	implicit none
	real(8) rh_now
	type(s1d_type)::s1d_tmp
	real(8) mw
	
	s1d_tmp=spp_new%fma_star
	s1d_tmp%xb=spp_new%fma_star%fx
	s1d_tmp%fx=spp_new%fma_star%xb
	s1d_tmp%xmin=minval(spp_new%fma_star%fx)
	s1d_tmp%xmax=maxval(spp_new%fma_star%fx)
	mw=2*spp_new%mbh/m0_cl
	!call s1d_tmp%print("s1d_tmp")
	if(mw<s1d_tmp%xmin)then
		rh_now=log10((mw*3/spp_new%spt_rho_rmin/4/pi)**(1d0/3d0))
	else
		call s1d_tmp%get_value_l(mw,rh_now)	
	end if 

end subroutine
 
subroutine write_down_oe_attributions(group_id,gname,star_type_number,oe)
	use md_hdf5
	use com_main_gw
	use md_dms_saving_data
	implicit none
	integer(HID_T)::group_id
	type(hdf5_group_type)::hg
	character*(*) gname
	integer star_type_number
	logical,external::test_if_component_exists
	type(obj_events)::oe

	if(test_if_component_exists(star_type_number))then
		!if(ctl%idxsbh.ne.-1)then
			call hg%create(group_id,trim(adjustl(gname)))
			call write_down_se_attributions(oe%se_emax,hg%group_id,"emax")
			if(ctl%include_loss_cone.ge.1)then
				call write_down_se_attributions(oe%se_lc,hg%group_id,"ls")
				call write_down_se_attributions(oe%se_td,hg%group_id,"td")
				if(oe%se_lc%nw>0)then
					call oe%se_lc%fdstr_x%save_hdf5(hg%group_id,"ls_x_dstr")
					call oe%se_lc%fdstr_m%save_hdf5(hg%group_id,"ls_m_dstr")
				end if
				if(oe%se_td%nw>0)then
					call oe%se_td%fdstr_x%save_hdf5(hg%group_id,"td_x_dstr")
					call oe%se_td%fdstr_m%save_hdf5(hg%group_id,"td_m_dstr")
				end if
			end if
			if(ctl%gw_radiation_otby.ge.1)then
				call write_down_se_attributions(oe%se_emris,hg%group_id,"emris")
				if(oe%se_emris%nw>0)then
					call oe%se_emris%fdstr_x%save_hdf5(hg%group_id,"emris_x_dstr")
					call oe%se_emris%fdstr_m%save_hdf5(hg%group_id,"emris_m_dstr")
					call oe%fd_emris_ecc%save_hdf5(hg%group_id,"emris_1-e_dstr")
					call oe%fd_emris_nxj_ir%save_hdf5(hg%group_id,"emris_nxj_ir")	
				end if
			end if
			call hg%close()
		!end if
	end if
end subroutine
subroutine collection_event_numbers(nr,nw)
	use MPI_comu
	implicit none
	real(8) nw, nr
	call collection_and_avg_real(nr)
	call collection_and_avg_real(nw)
end subroutine
subroutine write_down_se_attributions(se,group_id, aname)
	use md_hdf5
	use com_main_gw
	use md_dms_saving_data
	integer(HID_T)::group_id,attr_id  
	type(snap_event)::se
	character*(*) aname
	real(8) datas(3)
	datas=(/se%n, se%nw,se%rate/)
	call add_attr_dble_arr(group_id,attr_id, aname, datas,3)
end subroutine
 

subroutine output_dms_hdf5_pdf(dm, fl)
	use md_hdf5	
	use com_main_gw
	use md_dms_saving_data
	implicit none
	type(diffuse_mspec)::dm
	character*(*) fl
	integer i, error
	character*(4) tmpi
	type(hdf5_file_type)::hdf5_file
	type(hdf5_table_type)::hdf5_table
	INTEGER(HID_T) :: group_id,sub_group_id     ! File identifier	
    type(s1d_type)::s1d_trlx
	real(8) trlx_rh

	!print*, "0000",trim(adjustl(fl))//".hdf5"
	call hdf5_file%open(trim(adjustl(fl))//".hdf5") 
	call write_down_attributions(dm,hdf5_file%file_id)  
	do i=1, dm%n 
		write(unit=tmpi, fmt="(I4)") i
		CALL h5gcreate_f(hdf5_file%file_id, trim(adjustl(tmpi)), group_id, error)
		call save_dms_hdf5_pdf(dm%mb(i)%all, group_id, "all") 
		call save_dms_hdf5_pdf(dm%mb(i)%star, group_id, "star")
		call save_dms_hdf5_pdf(dm%mb(i)%bd, group_id, "bd")
		call save_dms_hdf5_pdf(dm%mb(i)%sbh, group_id, "sbh")
		call save_dms_hdf5_pdf(dm%mb(i)%wd, group_id, "wd")
		call save_dms_hdf5_pdf(dm%mb(i)%ns, group_id, "ns") 
		call save_dms_hdf5_pdf(dm%mb(i)%rg, group_id, "rg")
		call save_dms_hdf5_pdf(dm%mb(i)%dark_matter, group_id, "dark_matter")
		call save_dms_hdf5_pdf(dm%mb(i)%nakedHe, group_id, "Naked_HeliumStar")
		CALL h5gclose_f(group_id, error)
	end do
	!print*, "***********"
	call save_dms_hdf5_pdf(dm%all%all, hdf5_file%file_id, "all")
	call save_dms_hdf5_pdf(dm%all%star, hdf5_file%file_id, "star")
	call write_down_oe_attributions(hdf5_file%file_id,"oe_star",star_type_ms,oe_star)
	
	call save_dms_hdf5_pdf(dm%all%sbh, hdf5_file%file_id, "sbh")
	call write_down_oe_attributions(hdf5_file%file_id,"oe_sbh",star_type_bh,oe_sbh)
	call write_down_oe_attributions(hdf5_file%file_id,"oe_wd",star_type_wd,oe_wd)
	call write_down_oe_attributions(hdf5_file%file_id,"oe_ns",star_type_ns,oe_ns)
	call write_down_oe_attributions(hdf5_file%file_id,"oe_bd",star_type_bd,oe_bd) 
	call write_down_oe_attributions(hdf5_file%file_id,"oe_rg",star_type_rg,oe_rg)

	call save_dms_hdf5_pdf(dm%all%wd, hdf5_file%file_id, "wd")
	call save_dms_hdf5_pdf(dm%all%ns, hdf5_file%file_id, "ns")
	call save_dms_hdf5_pdf(dm%all%bd, hdf5_file%file_id, "bd") 
	call save_dms_hdf5_pdf(dm%all%rg, hdf5_file%file_id, "rg")
	call save_dms_hdf5_pdf(dm%all%nakedHe, hdf5_file%file_id, "Naked_HeliumStar") 

	call h5gcreate_f(hdf5_file%file_id, "dej", sub_group_id, error)
	!print*, "????"
    call output_de_hdf5(dm,sub_group_id)
    call h5gclose_f(sub_group_id, error) 
	call h5gcreate_f(hdf5_file%file_id, "spp", sub_group_id, error)
    call output_spp_hdf5(sub_group_id)
    call h5gclose_f(sub_group_id, error) 
	call save_dms_aux_data_hdf5(dm,spp_new,hdf5_file%file_id) 
	
	CALL hdf5_file%close() 
end subroutine
 
subroutine output_spp_hdf5(group_id)
	use md_hdf5
	use md_star_pot
	implicit none
	integer i
	character*(4) tmpi
	integer(HID_T)::group_id,attr_id
	real(8) mbhmass(n_record)
	!print*, "n_record=",n_record
	if(n_record.eq.0) return
	do i=1, n_record
		write(unit=tmpi,fmt="(I4)") i 
		call spp_record(i)%fphi_star%save_hdf5(group_id,   trim(adjustl(tmpi))//"_fphi")
		call spp_record(i)%frho_star%save_hdf5(group_id,   trim(adjustl(tmpi))//"_frho")
		call spp_record(i)%fgx_ir%save_hdf5(group_id,   trim(adjustl(tmpi))//"_fgx_ir")
	end do
	mbhmass(1:n_record)=spp_record(1:n_record)%mbh_dmless
	call add_attr_dble_arr(group_id,attr_id,"mbh_mass",mbhmass,n_record)
end subroutine
 
subroutine 	output_de_hdf5(dm,group_id)
    use md_hdf5
    use com_main_gw
    implicit none
    type(diffuse_mspec)::dm
    integer(HID_T):: sub_group_id, group_id
    type(s2d_type)::s2d
    integer i
	real(8) smin,smax,tmin,tmax
    if(dm%mb(1)%dc%s2_de_0%nx==0)return
    call dm%mb(1)%dc%s2_de_0%save_hdf5(group_id, "de_0_1")
    select case(ctl%ebin_type)
	case(ebin_type_log)
		smin=log10emin_factor
		smax=log10emax_factor 
	end select
	select case(ctl%jbin_type) 
	case(jbin_type_log)
		tmin=log10jmin_value
		tmax=log10jmax_value
	end select
	call s2d%init(dm%df_coe_bins, dm%df_coe_bins, smin,smax,tmin,tmax, sts_type_dstr) 
	s2d%xcenter=dm%mb(1)%dc%s2_de_110%xcenter
	s2d%ycenter=dm%mb(1)%dc%s2_de_110%ycenter
	s2d%fxy=0
	s2d%fxy=dm%dc0%s2_de_0%fxy
	do i=1, dm%n
		s2d%fxy=s2d%fxy+dm%mb(1)%mc/dm%mb(i)%mc*dm%mb(i)%dc%s2_de_110%fxy
	end do
	call s2d%save_hdf5(group_id, "de1")
	call dm%mb(1)%dc%s2_de_110%save_hdf5(group_id, "de_110_1")

    if(dm%n>1)then 
		call s2d%init(dm%df_coe_bins, dm%df_coe_bins, smin,smax,tmin,tmax, sts_type_dstr)
 		s2d%xcenter=dm%mb(1)%dc%s2_de_110%xcenter
		s2d%ycenter=dm%mb(1)%dc%s2_de_110%ycenter
		s2d%fxy=0
		s2d%fxy=dm%dc0%s2_de_0%fxy
		do i=1, dm%n
			s2d%fxy=s2d%fxy+dm%mb(2)%mc/dm%mb(i)%mc*dm%mb(i)%dc%s2_de_110%fxy
		end do
		call s2d%save_hdf5(group_id, "de2")
        
		call dm%mb(2)%dc%s2_de_110%save_hdf5(group_id, "de_110_2")

    end if
    
	call dm%dc0%s2_dee%save_hdf5(group_id, "dee")
	call dm%dc0%s2_dej%save_hdf5(group_id, "dej")
	call dm%dc0%s2_djj%save_hdf5(group_id, "djj")

        s2d%fxy=dm%dc0%s2_dj_rest%fxy
        do i=1, dm%n
            s2d%fxy=s2d%fxy+(dm%mb(1)%mc+dm%mb(i)%mc)/dm%mb(i)%mc/2d0 &
                *dm%mb(i)%dc%s2_dj_111%fxy
        end do
	call dm%mb(1)%dc%s2_dj_111%save_hdf5(group_id,"dj1_111")      	   
	call s2d%save_hdf5(group_id, "dj1")
    call dm%dc0%s2_dj_rest%save_hdf5(group_id,"dj_rest") 
			!print*, "dfsdfsdf"
    if(dm%n>1)then        
            s2d%fxy=dm%dc0%s2_dj_rest%fxy
            do i=1, dm%n
                s2d%fxy=s2d%fxy+(dm%mb(2)%mc+dm%mb(i)%mc)/dm%mb(i)%mc/2d0 &
                    *dm%mb(i)%dc%s2_dj_111%fxy
            end do
            call s2d%save_hdf5(group_id, "dj2")        
			call dm%mb(2)%dc%s2_dj_111%save_hdf5(group_id,"dj2_111")      
    end if
end subroutine

subroutine read_dms_hdf5_pdf(so, group_id, str_)
	use md_hdf5	
	use com_main_gw
	implicit none
	type(dms_stellar_object)::so
	character*(*) str_
	integer i, error
	character*(4) tmpi 
	INTEGER(HID_T) :: group_id,sub_group_id,s2d_group_id,lapl_id       ! File identifier	
	logical:: exist 
		CALL h5gopen_f(group_id, trim(adjustL(str_)), sub_group_id, error)
		if(error.eq.0)then
			call so%fden%read_hdf5(       sub_group_id,   "fden")
			call so%fmden%read_hdf5(     sub_group_id,"fmden")
			call so%fNa%read_hdf5(      sub_group_id,"fNa")
			call so%fslope%read_hdf5(   sub_group_id,"fslope")
			call so%fMa%read_hdf5(      sub_group_id,"fMa") 
			call so%fden_simu%read_hdf5(sub_group_id, "fden_simu")	 
			so%n_real=1
			so%n=1
		else
			so%n_real=0
			so%n=0 
		end if
		CALL h5gclose_f(sub_group_id, error)
	!end if
end subroutine

subroutine save_dms_hdf5_pdf(so, group_id, str_)
	use md_hdf5	
	use com_main_gw
	implicit none
	type(dms_stellar_object)::so
	character*(*) str_
	integer i, error
	character*(4) tmpi 
	INTEGER(HID_T) :: group_id,sub_group_id       ! File identifier	 
	if(so%n>0)then
		CALL h5gcreate_f(group_id, trim(adjustL(str_)), sub_group_id, error)
		call so%fden%save_hdf5(       sub_group_id,   "fden")
		!print*, "1.1"
		if(ctl%barge_evl_method.eq.barge_evl_method_direct)then
			call so%fden_simu%save_hdf5(  sub_group_id,   "fden_simu")
		end if
		if(ctl%include_stellar_evolution.ge.1)then
			call so%fm_dstr%save_hdf5(sub_group_id, "fm_dstr")
			call so%fr_dstr%save_hdf5(sub_group_id, "fr_dstr")
		end if
		call so%fmden%save_hdf5(     sub_group_id,"fmden")
		call so%fNa%save_hdf5(      sub_group_id,"fNa") 
		call so%fslope%save_hdf5(   sub_group_id,"fslope")
		call so%fMa%save_hdf5(      sub_group_id,"fMa") 
		
		select case(ctl%barge_evl_method)
		case(barge_evl_method_grid_2d) 
			call so%gxj_ir%save_hdf5(sub_group_id, "gxj_ir") 
			call so%nxj_ir%save_hdf5(sub_group_id,"nxj_ir") 
			call so%nx_ir%save_hdf5(sub_group_id,"nx_ir") 
			call so%barge_ir%save_hdf5(    sub_group_id,"fgx_ir")		 
		 
		end select 
		CALL h5gclose_f(sub_group_id, error)
	end if
end subroutine  

subroutine output_event_datas_hdf5(ed,fl,star_type)
	use md_event_datas
	use md_hdf5
	use md_star_pot
	implicit none
	character*(*) fl
	type(event_data)::ed!,em_copy
	integer(HID_T)::attr_id
	type(hdf5_file_type)::hdf5_file
	integer star_type
	! logical::isGasStar

	call hdf5_file%open(trim(adjustl(fl))//".hdf5")
	!em_copy=em
	call save_emris_hdf5(ed%emris,   hdf5_file%file_id,star_type)
	call save_td_hdf5(ed%td,   hdf5_file%file_id,star_type)
	call hdf5_file%close()
end subroutine 
#endif

