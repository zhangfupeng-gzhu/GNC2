

subroutine input_chains_bin(sps, fl)
	use md_chain
	use model_basic
	implicit none
	character*(*) fl
	integer i,n
	type(chain_type)::sps
	type(chain_pointer_type),pointer::pt
    integer flag
    integer,parameter::flag_sg=1 
	open(unit=999,file=trim(adjustl(fl))//".bin", form='unformatted',access='stream',status='old')
	read(unit=999) n
    read(999) ctl%current_version_number
	call input_random_seeds(999)
    print*, "input:chain length=", n
    call sps%init(n)
    pt=>sps%head
    do while(associated(pt))
        read(999) flag
        select case(flag)
        case(flag_sg)
            allocate(particle_sample_type::pt%ob)
            call pt%ob%read_info(999,ctl%current_version_number)
        
        end select        
        pt=>pt%next
    end do
	close(unit=999)
end subroutine
subroutine output_chains_bin(sps, fl)
	use md_chain
	use model_basic
	implicit none
	character*(*) fl
	integer i,n
	type(chain_type)::sps
	type(chain_pointer_type),pointer::pt
    integer,parameter::flag_sg=1 

	open(unit=999,file=trim(adjustl(fl))//".bin", form='unformatted',access='stream')
	
    call sps%get_length(n)
    print*, "output:chain length=",n
    write(999) n
    write(999) ctl%current_version_number
    call output_random_seeds(999)
    !print*, "n=",n
    pt=>sps%head
    do while(associated(pt))
        select type (ca=>pt%ob)
        type is(particle_sample_type)
            write(999) flag_sg
            call ca%write_info(999,ctl%current_version_number)
        end select        
        pt=>pt%next
    end do
	close(unit=999)
end subroutine


subroutine output_particle_sams_arr_bin(bksps, fl)
	use md_bk_species
	use model_basic,only:ctl
	implicit none
	character*(*) fl
	integer i,n
	type(particle_samples_arr_type)::bksps
	open(unit=999,file=trim(adjustl(fl))//".bin", form='unformatted',access='stream')
	write(unit=999) ctl%current_version_number
	call output_random_seeds(999)

	write(unit=999) bksps%n
	!print*, "output_sams_arr_bin:", sps%n
	do i=1, bksps%n
		call bksps%sp(i)%write_info(999,ctl%current_version_number)
	end do
	close(unit=999)
end subroutine

subroutine input_particle_sams_arr_bin(bksps, fl)
	use md_bk_species
	use model_basic,only:ctl
	implicit none
	character*(*) fl
	integer i,n
	type(particle_samples_arr_type)::bksps
	open(unit=999,file=trim(adjustl(fl))//".bin", form='unformatted',access='stream', status='old')
	read(unit=999) ctl%current_version_number
	call input_random_seeds(999)
	read(unit=999) bksps%n
	!print*, "output_sams_arr_bin:", sps%n
	call bksps%init(bksps%n)
	do i=1, bksps%n
		call bksps%sp(i)%read_info(999,ctl%current_version_number)
	end do
	close(unit=999)
end subroutine

 
subroutine gethering_samples(fdir, nsnap, bksar, ex)
	use com_main_gw
	 
	implicit none
	integer nsnap, nvalid,i,n
	type(particle_samples_arr_type)::bksar
	!type(chain_type)::bks
	!type(chain_type)::bys
	character*(*) fdir
	type(particle_samples_arr_type),allocatable::smsa(:)
	type(chain_type),allocatable::sms(:)
	
	character*(6) tmprid, tmpspid
	logical ex
	character*(9) str_
	character*(200) fl,fltpcs

	allocate(smsa(ctl%ntask_total) )
	allocate(sms(ctl%ntask_total) ) 

	write(unit=tmpspid,fmt="(I6)") nsnap

	nvalid=0
	print*, "ctl%ntask_total=",ctl%ntask_total
	do i=1, ctl%ntask_total
		
		write(unit=tmprid,fmt="(I6)") i
		str_=trim(adjustl(tmprid))//"_"//trim(adjustl(tmpspid))
		fl=trim(adjustl(fdir))//"/bin/single/samchn"//trim(adjustl(str_))
		fltpcs=trim(adjustl(fdir))//"/bin/single/event_scol"//trim(adjustl(str_))
		inquire(file=trim(adjustl(fl))//".bin",exist=ex)
		if(ex)then
			nvalid=nvalid+1
			!call smsa(nvalid)%input_bin(trim(adjustl(fl)))
			!print*, "arr nvalid=",nvalid
			call input_chains_bin(sms(nvalid),trim(adjustl(fl))) 
            call all_chain_to_arr_single(sms(nvalid),smsa(nvalid))
			call sms(nvalid)%destory()
		!	print*, "ch nvalid=",nvalid		
		else
			print*, trim(adjustl(fl))//" does not exist"
            return
		endif
	end do
	print*, "sms chain gathering finished"

	call smmerge_arr_single(smsa,nvalid,bksar)
	print*, "bksar%n=",bksar%n
	nvalid=0
    
	deallocate(sms)
	deallocate(smsa)
end subroutine
subroutine gethering_samples_ini(fdir, bksar, ex)
	use com_main_gw
	 
	implicit none
	integer nsnap, nvalid,i,n
	type(particle_samples_arr_type)::bksar
	!type(chain_type)::bks
	!type(chain_type)::bys
	character*(*) fdir
	type(particle_samples_arr_type),allocatable::smsa(:)
	type(chain_type),allocatable::sms(:)
	
	character*(6) tmprid
	logical ex
	character*(9) str_
	character*(200) fl,fltpcs

	allocate(smsa(ctl%ntask_total) )
	allocate(sms(ctl%ntask_total) ) 


	nvalid=0
	print*, "ctl%ntask_total=",ctl%ntask_total
	do i=1, ctl%ntask_total
		
		write(unit=tmprid,fmt="(I6)") i
		str_=trim(adjustl(tmprid))
		fl=trim(adjustl(fdir))//"/bin/single/samchn"//trim(adjustl(str_))
		fltpcs=trim(adjustl(fdir))//"/bin/single/event_scol"//trim(adjustl(str_))
		inquire(file=trim(adjustl(fl))//".bin",exist=ex)
		if(ex)then
			nvalid=nvalid+1
			!call smsa(nvalid)%input_bin(trim(adjustl(fl)))
			!print*, "arr nvalid=",nvalid
			call input_chains_bin(sms(nvalid),trim(adjustl(fl))) 
            call all_chain_to_arr_single(sms(nvalid),smsa(nvalid))
			call sms(nvalid)%destory()
		!	print*, "ch nvalid=",nvalid		
		else
			print*, trim(adjustl(fl))//" does not exist"
            return
		endif
	end do
	print*, "sms chain gathering finished"

	call smmerge_arr_single(smsa,nvalid,bksar)
	print*, "bksar%n=",bksar%n
	nvalid=0
    !bysar%n=0
   
	deallocate(sms)
	deallocate(smsa)
end subroutine

subroutine output_diffuse_mspec_bin(fl)
	use model_basic
	use com_main_gw!,only:id_saver
	use md_star_pot
	 
	use md_mbh_evl_acc
	use md_dms_saving_data
	implicit none
	character*(*) fl
	integer i 
	integer,parameter:: file_unit=209999
	

	open(unit=file_unit,file=trim(adjustl(fl)), access='stream', form='unformatted')
	write(unit=file_unit) ctl%current_version_number
	write(unit=file_unit) id_saver
	write(unit=file_unit) ctl%run_seq_idx, ctl%ts_spshot_dt,ctl%trlx_rh0, ctl%n_spshot_total,ctl%num_step_per_update
	write(unit=file_unit) ctl%run_snap_time_i, ctl%run_snap_time_f
	write(unit=file_unit) ctl%total_energy,ctl%energy_lost_emin
	write(unit=file_unit) ctl%m_bins

	
	write(unit=file_unit) ctl%ini_stellar_tot(1:n_tot_comp), ctl%ini_stellar_each_mass(1:n_tot_comp,1:ctl%m_bins)

	!write(unit=file_unit) ctl%min_mbh_factor,ctl%max_mbh_factor, ctl%max_logrmin, ctl%log10rmin_factor
	write(unit=file_unit) m0_cl,r0_cl

	
	call save_spp_tables_bin(spp_new,file_unit)
	call save_spp_tables_bin(spp_old,file_unit)
	
	write(unit=file_unit) sample_logemin,sample_logemax,sample_logrmin,sample_logrmax, &
	sample_nxgx_logemin,sample_nxgx_logemax		
	
	write(unit=file_unit) dms%n, dms%idx_ref, dms%df_coe_bins, dms%dstr_bins_r, dms%dstr_bins_e, dms%emin, &
		dms%emax, dms%jmin,dms%jmax, dms%mtot, dms%v0, dms%n0, dms%r0_cl, &
		dms%logrmin,dms%logrmax,dms%x_boundary,dms%ebin_type,dms%jbin_type,dms%e_iregular!,&

	write(file_unit) dms%alpha_r, dms%jc,dms%fr_phi, dms%rc
	
	write(file_unit) dms%frc_x
	
	if(ctl%fden_ana_est_method.eq.fden_ana_est_method_2d)then
		write(file_unit) dms%jc_sample_erange
	end if
	!if(ctl%barge_grid_type.ne.barge_grid_type_regular)then
		write(file_unit) dms%barp_ir,dms%dlxb_ir
	!else
	!	write(file_unit) dms%barp
	!end if

	do i=1, dms%n
		call dms%mb(i)%write_mb(file_unit)
	end do
	call dms%dc0%write_grid(file_unit)
	!print*, dms%mb(1)%dc%s2_djj%fxy(1,1:10)
	!print*, dms%dc0%s2_djj%fxy(1,1:10)
	!stop

	!write(file_unit) ctl%ini_nx_log, ctl%ini_ge
	do i=1, ctl%m_bins
		write(file_unit) ctl%ini_frho(i), ctl%ini_ge(i)
	end do
	write(file_unit) ctl%ini_frho_tot,ctl%ini_nx_tot,ctl%ini_ge_tot

	write(file_unit) dms%all%all%nx, dms%all%all%fna, dms%all%all%fna_simu
	write(file_unit) dms%all%all%fmden,dms%all%all%fden
	write(file_unit) dms%all%all%fma, dms%all%all%fma_simu,dms%all%all%barge_ir
	do i=1, n_tot_comp_sg
		write(file_unit) dms%all%dsp(i)%p%fden
	end do
	write(file_unit) dms%rp,dms%ra,dms%pd
	
	if(ctl%enable_evl_mbh.ge.1)then
		call write_mmg(file_unit)
	end if
	write(file_unit) oe_star,oe_rg,oe_sbh,oe_wd,oe_ns,oe_bd
	!print*, "6", oe_bd%se_td%n
	!if(ctl%gw_radiation_otby.ge.1)then
	!	write(file_unit) sample_prestep_within_emri
	!end if

	
	close(file_unit)

end subroutine
subroutine input_random_seeds(file_unit)
	use md_gaussian
	use model_basic,only:ctl
	use mpi_comu,only:rid
	implicit none
	integer file_unit,random_seed_size

	read(unit=file_unit) ctl%random_seed_size
	!print*, "random_seed_size=",ctl%random_seed
	if(allocated(ctl%random_seed))then
		deallocate(ctl%random_seed)
	end if
	allocate(ctl%random_seed(ctl%random_seed_size))
	read(unit=file_unit) ctl%random_seed(1:ctl%random_seed_size)

	call random_seed(size=random_seed_size)
	if(random_seed_size.eq.ctl%random_seed_size)then
		call random_seed(put=ctl%random_seed(1:ctl%random_seed_size))
	else
		print*, "random_seed size not matched, use own random seed"
		call same_random_seed(ctl%seed_value)
	end if
	read(unit=file_unit) flag, v1, v2, s
	print*, "rid, read in seed:=",rid, ctl%random_seed(1:ctl%random_seed_size)
end subroutine
subroutine output_random_seeds(file_unit)
	use md_gaussian
	use model_basic,only:ctl
	use mpi_comu,only:rid
	implicit none
	integer file_unit

	call random_seed(size=ctl%random_seed_size)
	write(unit=file_unit) ctl%random_seed_size
	if(allocated(ctl%random_seed))then
		deallocate(ctl%random_seed)
	end if
	allocate(ctl%random_seed(ctl%random_seed_size))
	call random_seed(get=ctl%random_seed(1:ctl%random_seed_size))
	write(unit=file_unit) ctl%random_seed(1:ctl%random_seed_size)
	write(unit=file_unit) flag, v1, v2, s
	print*, "rid, write down seed:=",rid, ctl%random_seed(1:ctl%random_seed_size)
end subroutine

subroutine input_diffuse_mspec_bin(fl)
	!use model_basic
	!use MPI_comu,only:rid
	use com_main_gw
	use md_star_pot
	 
	use md_mbh_evl_acc
	use md_dms_saving_data
	implicit none
	character*(*) fl
	integer i
	integer,parameter:: file_unit=209999
	real(8) total_time,burn_in_time,update_dt,tnr, ts_spshot_dt
	real(8) mbh_org
	integer n_spshot_total,num_step_per_update

	open(unit=file_unit,file=trim(adjustl(fl)), access='stream', form='unformatted', status='old')
	read(unit=file_unit) ctl%current_version_number
	if(rid.eq.0)then
		print*, "current_version_number=",ctl%current_version_number
	end if
	read(unit=file_unit) id_saver
	read(unit=file_unit) ctl%run_seq_idx, ts_spshot_dt,ctl%trlx_rh0, n_spshot_total,num_step_per_update
	!print*, "trlx_rh0=",ctl%trlx_rh0
	!read(*,*)
	read(unit=file_unit) ctl%run_snap_time_i, ctl%run_snap_time_f
	read(unit=file_unit) ctl%total_energy,ctl%energy_lost_emin
	read(unit=file_unit) ctl%m_bins
	
	if(allocated(ctl%ini_stellar_tot)) then
		deallocate(ctl%ini_stellar_tot,ctl%ini_stellar_each_mass)
	end if
	allocate(ctl%ini_stellar_tot(n_tot_comp),ctl%ini_stellar_each_mass(n_tot_comp,ctl%m_bins))
	read(unit=file_unit) ctl%ini_stellar_tot, ctl%ini_stellar_each_mass
	!read(unit=file_unit) ctl%min_mbh_factor,ctl%max_mbh_factor, ctl%max_logrmin, ctl%log10rmin_factor
	!call get_logrmin_from_mbh(mbh,(ctl%denhenmodel_mtot*m0_cl),  ctl%log10rmin_factor)
	read(unit=file_unit) m0_cl,r0_cl
	
	call read_spp_tables_bin(spp_new,file_unit)
	call read_spp_tables_bin(spp_old,file_unit)
	
	read(unit=file_unit) sample_logemin,sample_logemax,sample_logrmin,sample_logrmax, &
		sample_nxgx_logemin,sample_nxgx_logemax
	if(rid.eq.0)then
		print*, "input_diffuse_mspec_bin: sample_logemin,emax=",sample_logemin,sample_logemax
	end if
	read(unit=file_unit) dms%n, dms%idx_ref, dms%df_coe_bins, dms%dstr_bins_r, dms%dstr_bins_e, dms%emin, &
		dms%emax, dms%jmin,dms%jmax, dms%mtot, dms%v0, dms%n0, dms%r0_cl, &
		dms%logrmin,dms%logrmax,dms%x_boundary,dms%ebin_type,dms%jbin_type,dms%e_iregular!,&

	call get_ini_ebounds()
	! print*, "dminit",dms%n
	call dms%init(dms%n)
	! print*, "stoinit"
	call init_stellar_obj_rtables(dms)
	! print*, "---xxxxxx"
	read(file_unit) dms%alpha_r, dms%jc,dms%fr_phi, dms%rc

	read(file_unit) dms%frc_x
	if(ctl%fden_ana_est_method.eq.fden_ana_est_method_2d)then
		read(file_unit) dms%jc_sample_erange
	end if
	!call dms%alpha_r%print("read: alpha_r")
	!print*, size(dms%mb)
	!print*, "dms%n=",dms%n
	!emin_factor=10**dms%logemin; emax_factor=10**dms%logemax

	! print*, "00xxxxx"
	
	!call init_diffuse_mspec_rtables(dms)
	!print*, 'dc'
	call init_dms_dc(dms)  
	!call init_diffuse_mspec_etables(dms)

	call set_mass_bin_mass_given(dms, ctl%bin_mass, ctl%bin_mass_m1,&
	ctl%bin_mass_m2, ctl%asymptot_ini, ctl%m_bins)
	! print*, "gx"
	call init_diffuse_mspec_gxtables(dms)

	!if(ctl%barge_grid_type.ne.barge_grid_type_regular)then
		read(file_unit) dms%barp_ir,dms%dlxb_ir
	!else
	!	read(file_unit) dms%barp
	!end if
	!if(ctl%barge_grid_type.eq.barge_grid_type_iregular)then
    !    call set_etable_ir_bins()
    !end if

	!print*, dms%all%all%barge_ir%nbin
	!print*, dms%mb(1)%all%barge_ir%nbin
	!print*, dms%mb(1)%star%barge_ir%nbin
	!print*, dms%mb(1)%sbh%barge_ir%nbin
	

	!emax_dstr_factor=10**dms%fphi_star%fx(1)
	!emin_dstr_factor=10**dms%fphi_star%fx(dms%fphi_star%nbin)
	!call dms%mb(1)%star%barge_ir%print("br0")
	!if(ctl%barge_grid_type.eq.barge_grid_type_iregular)then
	!	call set_gx_nx_ranges_ir(dms)
	!end if
	!call dms%mb(1)%star%barge_ir%print("br1")
	!read(*,*)
	!call set_gx_ranges_ir(dms)
	!print*, "11"
	!print*, "mb"
	do i=1, dms%n
		call dms%mb(i)%read_mb(file_unit)
	end do
	! print*, "read 1 finished"
	call dms%dc0%read_grid(file_unit)
	!call dms%dc0_bk%read_grid(999)
	!print*, "read 2 finished"
	!read(file_unit) 
	!print*, dms%mb(1)%dc%s2_djj%fxy(1,1:10)
	!print*, dms%dc0%s2_djj%fxy(1,1:10)
	!stop
	!call dms%mb(1)%star%barge_ir%print("br0")
	!read(*,*)
	
	!read(file_unit) ctl%ini_nx_log, ctl%ini_ge
	do i=1, ctl%m_bins
		read(file_unit) ctl%ini_frho(i), ctl%ini_ge(i)
	end do
	read(file_unit) ctl%ini_frho_tot,ctl%ini_nx_tot,ctl%ini_ge_tot
	!print*, "read 4 finished"
	read(file_unit) dms%all%all%nx, dms%all%all%fna, dms%all%all%fna_simu
	read(file_unit) dms%all%all%fmden,dms%all%all%fden
	read(file_unit) dms%all%all%fma, dms%all%all%fma_simu,dms%all%all%barge_ir
	do i=1, n_tot_comp_sg
		read(file_unit) dms%all%dsp(i)%p%fden
	end do
	read(file_unit) dms%rp,dms%ra,dms%pd
	
	
	!if(ctl%dc_grid_type.eq.dc_grid_irregular)then
	!	read(file_unit) common_jc, common_rp, common_ra, common_pd
	!end if
	

	! print*, "00000"
	if(ctl%enable_evl_mbh.ge.1)then
		call read_mmg(file_unit)
	end if
	! print*, "read 5 finished"
	read(file_unit) oe_star, oe_rg, oe_sbh, oe_wd, oe_ns, oe_bd
	! print*, "----"

	close(unit=file_unit)

end subroutine
