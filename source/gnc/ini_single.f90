module md_ini_data
	use model_basic,only:ctl,log_clone_bd_sep, clone_e0_factor
	use MPI_comu,only:rid
	implicit none
	real(8),allocatable:: cx_tmp(:),fx_tmp(:)
	real(8) fmax
	integer clone_amplifier
contains
	subroutine set_ini_dstr(midx)
		integer i, midx
		integer max_lvl
		real(8) get_clone_deep
		if(allocated(fx_tmp))deallocate(fx_tmp)
		allocate(fx_tmp(ctl%ini_nx_tot%nbin))
		fmax=0
		if(ctl%clone_scheme.ge.1)then
			
			clone_amplifier=ctl%clone_factor(midx)
			! if(allocated(fx_tmp))deallocate(fx_tmp)
			
			do i=1, ctl%ini_nx_tot%nbin
				if(ctl%ini_nper_bin(midx)%fx(i)<0.5)then
					fx_tmp(i)=1d-80
				else
					fx_tmp(i)=10**ctl%ini_nx_log(midx)%fx(i)
				end if
				!print*, "xb,fx=",ctl%ini_nx_tot%xb(i), fx_tmp(i)
			end do
			do i=1, ctl%ini_nx_tot%nbin						
				max_lvl=get_clone_deep(10**ctl%ini_nx_tot%xb(i),log_clone_bd_sep,  clone_e0_factor)
				if(max_lvl>=0)then
					fmax=max(fmax,fx_tmp(i)*dble(clone_amplifier)**dble(max_lvl))
				else
					fmax=max(fmax,fx_tmp(i))
				end if
			end do
		else 
			do i=1, ctl%ini_nx_tot%nbin
				fx_tmp(i)=10**ctl%ini_nx_log(midx)%fx(i)
			end do
			do i=1, ctl%ini_nx_tot%nbin						
				fmax=max(fmax,fx_tmp(i))
			end do
		end if
		
		if(allocated(cx_tmp))deallocate(cx_tmp)
		allocate(cx_tmp(ctl%ini_nx_tot%nbin))
		if(ctl%clone_scheme.ge.1)then
			clone_amplifier=ctl%clone_factor(midx)
			
			! allocate(cx_tmp(ctl%ini_nx_tot%nbin))
			
			cx_tmp(1)=fx_tmp(1)
			do i=2, ctl%ini_nx_tot%nbin
				if(ctl%ini_nper_bin(midx)%fx(i)<0.5)then
					cx_tmp(i)=cx_tmp(i-1)+fx_tmp(i)
				else
					max_lvl=get_clone_deep(10**ctl%ini_nx_tot%xb(i),log_clone_bd_sep,  clone_e0_factor)
					cx_tmp(i)=cx_tmp(i-1)+fx_tmp(i)*dble(clone_amplifier)**dble(max_lvl)
				end if
				!print*, "xb,fx=",ctl%ini_nx_tot%xb(i), fx_tmp(i)
			end do
			cx_tmp=cx_tmp/cx_tmp(ctl%ini_nx_tot%nbin)
		else
			clone_amplifier=ctl%clone_factor(midx)
			! if(allocated(cx_tmp))deallocate(cx_tmp)
			
			cx_tmp(1)=fx_tmp(1)
			do i=2, ctl%ini_nx_tot%nbin
				if(ctl%ini_nper_bin(midx)%fx(i)<0.5)then
					cx_tmp(i)=cx_tmp(i-1)+fx_tmp(i)
				else 
					cx_tmp(i)=cx_tmp(i-1)+fx_tmp(i)
				end if 
			end do
			cx_tmp=cx_tmp/cx_tmp(ctl%ini_nx_tot%nbin)
		end if 
	end subroutine
end module


 
subroutine get_numbers_each_bin(nstellar_tot,nstellar_each_mass, n_comp, n_mass)
	use com_main_gw
	implicit none
	integer n_comp, n_mass
	integer i,j
	integer nstellar_tot(n_comp)
	integer nstellar_each_mass(n_comp, n_mass)
	do i=1, n_mass
		do j=1, n_comp
			nstellar_each_mass(j,i)=ctl%asymptot_ini(j+1,i)*ctl%bin_mass_particle_number(i)
			!print*, "i,j,n=",i,j,nstellar_each_mass(j,i)
		end do
		if(rid.eq.0)then
			print*, "i,n=",i,nstellar_each_mass(:,i)
		end if
	end do
	do j=1, n_comp
		nstellar_tot(j)=sum(nstellar_each_mass(j,:))
	end do
	!stop
end subroutine

subroutine get_init_samples(bksps_arr_ini)
	use com_main_gw
	use md_ini_data
	implicit none
	type(particle_samples_arr_type)::bksps_arr_ini
	integer i, j, k, nsg_tot, nby_tot
	!character*(4) tmpj
	integer nsg0, nby0,stellar_type
	
	call get_numbers_each_bin(ctl%ini_stellar_tot,ctl%ini_stellar_each_mass,n_tot_comp,ctl%m_bins)

	nsg_tot=0;
	nby_tot=0
	do i=1, n_tot_comp_sg
		call get_type_from_ctl_obidx_sg(i,stellar_type)
		nsg_tot=nsg_tot+ctl%ini_stellar_tot(i)
	end do
	
	if(rid.eq.0)then
		write(unit=*,fmt="(A4, 12A15)") "", "NStar", "NSbh", "NNs", "NWd", "NBd", "Nrg"
		write(*, fmt="(A4,10I15)") "TOT=", ctl%ini_stellar_tot(1:6)
		do i=1, ctl%m_bins 
			write(*, fmt="(I4,20I15)") i,ctl%ini_stellar_each_mass(1:6,i)
		end do
		print*, "start bksps_ini", nsg_tot
	end if
	call bksps_arr_ini%init(nsg_tot)

	!print*, "ctl%m_bins=",ctl%m_bins
    nsg0=0; nby0=0;
	do i=1, ctl%m_bins
		call set_ini_dstr(i)
		! print*, "i=",i
		do j=1, n_tot_comp_sg
			! print*, "j=",j
			call get_type_from_ctl_obidx_sg(j,stellar_type)	
			if(ctl%ini_stellar_each_mass(j,i)>0)then			
				do k=1, ctl%ini_stellar_each_mass(j,i)
					bksps_arr_ini%sp(k+nsg0)%obtype=stellar_type
					bksps_arr_ini%sp(k+nsg0)%obidx=j
					! print*, "k=",k
					call init_particle_sample_one_model_rnd(bksps_arr_ini%sp(k+nsg0),i,&
						ctl%ini_sample_sg_mode)
				end do
			end if
			nsg0=nsg0+ctl%ini_stellar_each_mass(j,i)
		end do
	end do

	call set_clone_weight_arr(bksps_arr_ini)
	call set_real_weight_arr(bksps_arr_ini)
end subroutine
logical function test_if_component_exists(component_type)
	use com_main_gw
	implicit none
	integer component_type
	integer j
	call get_obidx_from_type_sg(component_type,j)
	if(j<=0) then
		test_if_component_exists=.false.
		return
	end if
	!if(sum(ctl%ini_stellar_each_mass(j,:))>0)then
	if(ctl%exist_stellar_type(j).ge.1)then
		test_if_component_exists=.true.
	else
		test_if_component_exists=.false.
	end if
end function

subroutine set_chain_samples_single(cbk, bksps_arr)
	use com_main_gw
	IMPLICIT NONE
	integer i,j, flag, obtype, typeidx,n
	type(chain_pointer_type),pointer::ptbk, ptby
	type(particle_samples_arr_type)::bksps_arr
	type(chain_type)::cbk
	!integer,parameter::flag_sg=1,flag_by=2

   ! bksps_arr%n=10
    !bysps_arr%n=10
	call cbk%init(bksps_arr%n)
    ptbk=>cbk%head
	do i=1, bksps_arr%n
		if(.not.allocated(ptbk%ob)) allocate(particle_sample_type::ptbk%ob)
        select type (ca=>ptbk%ob)
        type is(particle_sample_type)
            ca=bksps_arr%sp(i)
            ca%create_time=0d0
            ca%simu_bgtime=0d0 
            ca%en0=ca%en
            ca%jm0=ca%jm
			if(ca%jm0<jmin_value)then
				print*, "jm0=",ca%jm0, "<=jmin_value=",jmin_value
				stop
			end if 
		end select
        ptbk=>ptbk%next
	end do
    
end subroutine

 
subroutine get_sample_ini_stellar_history(sp,tstart)
	use com_main_gw
	use md_stellar_history
	use md_stellar_evolution
	implicit none
	type(particle_sample_type)::sp
	integer kstar,flag
	real(8) tstart

	kstar=get_kstar_integer("MS")
	ctl%metal_z=0.02d0
	call sp%sh%get_history(kstar,sp%m,ctl%metal_z,12000d0,tstart,0) 
	call get_current_particle_stellar_info(sp,0d0,flag,.false.)
end subroutine
subroutine init_particle_sample_one_model_rnd(bkps, midx, ini_sample_sg_mode)
	use com_main_gw
	use md_ini_data
	implicit none
	type(particle_sample_type)::bkps
	real(8) ecmax,logpd, beta, GET_T_GW, ini_ex
	real(8),external:: rnd, gen_ran_from_dstr_consider_clone_cum,get_clone_deep
	real(8),external:: gen_ran_from_dstr_consider_clone
	real(8),external:: gen_ran_from_data, p_EJ_dmless_fast,r_c_iter
	real(8) rmax,rc,jc_dmless,pd_xy,rp_xy,ra_xy,jc_xy,radius 
	integer midx, maxlocation(1),  max_lvl
	integer i,ier,ini_sample_sg_mode 
	if(bkps%obtype.eq.0.or.bkps%obidx.eq.0)then
		print*, "error:particle type not assigned", bkps%obtype, bkps%obidx
		stop
	end if

	! print*, fx_tmp

100 select case(ini_sample_sg_mode)
	case(ini_sample_mode_given)
		bkps%m=ctl%bin_mass(midx)
		bkps%byot%ms%m=bkps%m
	
		bkps%byot%ms%obtype=bkps%obtype
		bkps%byot%ms%obidx=bkps%obidx

		if(ctl%include_stellar_evolution.ge.1.and.bkps%obtype.eq.star_type_ms)then
			call get_sample_ini_stellar_history(bkps,0d0)
		else
			call set_star_radius(bkps%byot%ms)
		end if
 
		if(bkps%byot%ms%radius.eq.0.and.bkps%obtype.eq.star_type_ms) then
			print*, "error!: ini", bkps%byot%ms%radius
			stop
		end if
	case(ini_sample_mode_mobse)
		select case(ctl%ini_mass_bin_mode)
		case(ini_mass_bin_mode_kroupa)
			if(bkps%obtype.ne.star_type_ms.and. bkps%obtype.ne.star_type_bd)then
				print*, "error! kroupa mass bin should all assign to stars or brown dwarf"
				stop
			end if 
			call get_one_kroupa_sample(bkps,ctl%bin_mass_m1(midx),ctl%bin_mass_m2(midx),0d0) 
		case(ini_mass_bin_mode_topheavy)
			if(bkps%obtype.ne.star_type_ms.and. bkps%obtype.ne.star_type_bd)then
				print*, "error! kroupa mass bin should all assign to stars or brown dwarf"
				stop
			end if 
			call get_one_topheavy_sample(bkps,ctl%bin_mass_m1(midx),ctl%bin_mass_m2(midx),0d0) 
		case(ini_mass_bin_mode_pow)
			if(bkps%obtype.ne.star_type_ms)then
				print*, "error! pow mass bin should all assign to stars"
				stop
			end if
			
		end select 
		bkps%byot%ms%m=bkps%m

		bkps%byot%ms%obtype=bkps%obtype
		bkps%byot%ms%obidx=bkps%obidx
		!bkps%byot%ms%radius=radius	
	
	case default
		print*, "error! ini_sample_sg_mode not defined", ini_sample_sg_mode
		stop
	end select
	!call get_mass_idx(bkps%m,midx)
	
	if(ctl%clone_scheme.ge.1)then
		 
		ini_ex=gen_ran_from_dstr_consider_clone(ctl%ini_nx_tot%xb,fx_tmp, &
			ctl%ini_nx_tot%y2, ctl%ini_nx_tot%nbin, ctl%ini_nx_tot%xmin,ctl%ini_nx_tot%xmax,&
		fmax, clone_e0_factor, clone_amplifier,  2)
		
	else
		ini_ex=gen_ran_from_data(ctl%ini_nx_tot%xb, fx_tmp,&
			ctl%ini_nx_tot%nbin, ctl%ini_nx_tot%xmin,ctl%ini_nx_tot%xmax, &
			fmax, 2)
	end if
	
	select case(ctl%ebin_type)
	case(ebin_type_log)
		bkps%en= -ctl%v0**2*10**ini_ex
		bkps%x=10**ini_ex
		call get_rmax_accurate(spp_new,  dms%fr_phi, ini_ex,rmax)
	 
	end select
	rc=r_c_iter(spp_new,bkps%x,ier) 
	jc_xy=jc_dmless(rc,spp_new)
	bkps%jc=jc_xy*(r0_cl*ctl%v0)
	
200		call set_jm_init(bkps)
 
    call get_rpra_dmless(spp_new, bkps%x, bkps%jm, jc_xy, &
                log10(rc), rmax, rp_xy,ra_xy)
	bkps%rp=rp_xy*r0_cl       
    bkps%ra=ra_xy*r0_cl
    pd_xy=p_EJ_dmless_fast(spp_new, bkps%x,bkps%jm,  jc_xy, rp_xy,ra_xy)
    bkps%period=pd_xy*r0_cl/ctl%v0

	bkps%jph=bkps%jm*bkps%jc 
	call init_particle_sample_common(bkps)
	if(ctl%include_loss_cone.ge.1)then 
		call get_sample_r_td(bkps)
		call get_sample_jlc(bkps%x,spp_new%mbh_dmless,bkps%r_lc/r0_cl,bkps%jc/(ctl%v0*r0_cl),spp_new,sample_jlc_dimless,ier)
		if(ier>0.or.ra_xy<bkps%r_lc/r0_cl)then 
			if(ra_xy<bkps%r_lc/r0_cl)then 
				goto 100
			else
				 
				goto 200
			end if
		end if 
	end if
	call get_mass_idx(bkps%m,midx)
	bkps%weight_n=ctl%Weight_n(midx)
	
end subroutine
subroutine get_one_kroupa_sample(sp,m1,m2,tstart)
	use com_main_gw
	use md_stellar_history
	use md_stellar_evolution
	implicit none
	type(particle_sample_type)::sp
	real(8) m, canonical_IMF,radius,m1,m2,tstart
	integer kstar,kwtype,flag
	
	m=canonical_IMF(imf_para_nt,imf_para_nc,imf_para_nq,m1,m2)
	
	sp%m=m
	kstar=get_kstar_integer("MS")
	 
	call sp%sh%get_history(kstar,m,ctl%metal_z,12000d0,tstart,0) 
	if(sp%sh%n>8)then 
		call simplfy_history(sp%sh) 
	end if 

	call get_current_particle_stellar_info(sp,0d0,flag,.false.)
 
end subroutine
subroutine simplfy_history(history)
	use md_stellar_history
	implicit none
	type(type_stellar_history)history, history_simplfy,history_tmp
	integer n, i,n_sm
	n=history%n
	history_tmp=history
	n_sm=1
	do i=1, n-1
		if(history%kwtype(i).ne.history%kwtype(i+1))then
			n_sm=n_sm+1
			goto 100
		end if
		if(abs(history%mass(i+1)-history_tmp%mass(n_sm))/history_tmp%mass(n_sm)>0.1)then
			n_sm=n_sm+1
			goto 100
		end if
		if(abs(history%radius(i+1)-history_tmp%radius(n_sm))/history_tmp%radius(n_sm)>0.1)then
			n_sm=n_sm+1
			goto 100
		end if
		if(history%ktype(i).ne.history%ktype(i+1))then
			n_sm=n_sm+1
			goto 100
		end if
		if(i.eq.n-1)then
			n_sm=n_sm+1
			goto 100
		end if
		cycle
100     history_tmp%kwtype(n_sm)=history%kwtype(i+1)
		history_tmp%mass(n_sm)=history%mass(i+1)
		history_tmp%ktype(n_sm)=history%ktype(i+1)
		history_tmp%radius(n_sm)=history%radius(i+1)
	end do
	! print*, "n, n_sm=",n,n_sm
	call history_simplfy%init(n_sm)
	n_sm=1
	history_simplfy%kwtype(1)=history%kwtype(1)
	history_simplfy%mass(1)=history%mass(1)
	history_simplfy%ktype(1)=history%ktype(1)
	history_simplfy%radius(1)=history%radius(1)

	do i=1, n-1
		if(history%kwtype(i).ne.history%kwtype(i+1))then
			n_sm=n_sm+1
			goto 200
		end if
		if(abs(history%mass(i+1)-history_simplfy%mass(n_sm))/history_simplfy%mass(n_sm)>0.1)then
			n_sm=n_sm+1
			goto 200
		end if
		if(abs(history%radius(i+1)-history_simplfy%radius(n_sm))/history_simplfy%radius(n_sm)>0.1)then
			n_sm=n_sm+1
			goto 200
		end if
		if(history%ktype(i).ne.history%ktype(i+1))then
			n_sm=n_sm+1
			goto 200
		end if
		if(i.eq.n-1)then
			n_sm=n_sm+1
			goto 200
		end if
		cycle
200     history_simplfy%kwtype(n_sm)=history%kwtype(i+1)
		history_simplfy%mass(n_sm)=history%mass(i+1)
		history_simplfy%ktype(n_sm)=history%ktype(i+1)
		history_simplfy%radius(n_sm)=history%radius(i+1)
		history_simplfy%time(n_sm)=history%time(i+1)
	end do

	history=history_simplfy
end subroutine
subroutine get_one_topheavy_sample(sp,m1,m2,tstart)
	use com_main_gw
	use md_stellar_history
	use md_stellar_evolution
	implicit none
	type(particle_sample_type)::sp
	real(8) m, triple_break_IMF,radius,m1,m2,tstart
	integer kstar,kwtype,flag
	
	m=triple_break_IMF(imf_para_nt,imf_para_nc,imf_para_nq,m1,m2,&
		topheavy_alpha,topheavy_xb)
	
	sp%m=m
	kstar=get_kstar_integer("MS")
	
	call sp%sh%get_history(kstar,m,ctl%metal_z,12000d0,tstart,0)
	 

	call get_current_particle_stellar_info(sp,0d0,flag,.false.)
 
end subroutine

subroutine get_jm_idx(jm, idx, rdx,evjum)
	use model_basic
	use md_coeff
	implicit none
	real(8) jm, evjum,jmin,  jmax,rdx
	integer idx 
	if(jm<jmin_value)then
		jm=2*jmin_value-jm 
	end if

	if(jm>jmax_value)then
		jm=jmax_value 
	end if 

	select case(ctl%jbin_type) 
	case(jbin_type_log) 
		evjum=log10(jm) 
	case default
		print*, "error! define jbtype", ctl%jbin_type
	end select

	select case(ctl%jbin_type) 
	case(Jbin_type_log)
		jmin=log10(jmin_value)
		jmax=log10(jmax_value) 
	end select 
	call return_idx(evjum, jmin,jmax, dms%df_coe_bins, idx, &
		coeff_sts_type_dc)

	if(idx<-1)then
		idx=1
		print*, "*****************evjum=",evjum, jmin, jmin_value, jmax_value, jm
		stop
	end if
	if(idx>dms%df_coe_bins)then
		idx=dms%df_coe_bins
	end if 
	rdx=(evjum-log10jmin_value)/dc_grid_ystep 
end subroutine
subroutine get_ex_idx_ir(ex,idx,rdx, even)
	use model_basic
	use md_coeff
	implicit none
	real(8) ex, even,rdx
	integer idx

	even=log10(ex)
	if(even>log10emax_factor) then
		even=log10emax_factor
	end if
	if(even<log10emin_factor)then
		even=log10emin_factor
	end if 

	call return_idx_ir(dms%dc0%s2_dee%xcenter,dms%df_coe_bins,even, idx,sts_type_dstr)
		 
	if(idx<-1)then
		idx=1
		print*, "even, ex=",even, ex
	end if
	if(idx>dms%df_coe_bins)then
		idx=dms%df_coe_bins
		print*, "even=",even
	end if
	if(ctl%method_interpolate.eq.method_int_linear)then 
		rdx=idx-1+(even-(dms%dlxb_ir%xb(idx)-dms%dlxb_ir%xsteps(idx)/2d0))&
			/dms%dlxb_ir%xsteps(idx)
		 
	else
		rdx=idx
	end if
end subroutine
subroutine get_ex_idx(ex,idx,rdx,even)
	use model_basic
	use md_coeff
	implicit none
	real(8) ex, even, smin, smax,rdx
	integer idx
	select case(ctl%ebin_type)
	case(ebin_type_log)
		even=log10(ex)
		smin=log10emin_factor; smax=log10emax_factor 
	case default
		print*, "get_ex_idx define"
		stop
	end select
	
	if(even>smax) then
		even=smax
	end if
	if(even<smin)then
		even=smin
	end if 
	call return_idx(even, smin, smax,dms%df_coe_bins, idx,coeff_sts_type_dc) 
	if(idx<-1)then
		idx=1
		print*, "even, ex=",even, ex
	end if
	if(idx>dms%df_coe_bins)then
		idx=dms%df_coe_bins
		print*, "even=",even
	end if


	select case(ctl%ebin_type)
	case(ebin_type_log)
		rdx=(even-log10emin_factor)/dc_grid_xstep 
	end select
	

end subroutine
subroutine init_particle_sample_common(bkps)
	use com_main_gw
	implicit none
	type(particle_sample_type)::bkps
	real(8),external:: rnd, gen_ran_from_data 

	call set_star_spin_random(bkps%byot%ms) 
	!end if
	bkps%en0=bkps%en 
    bkps%jm0=bkps%jm 
	bkps%id=gen_id()
	bkps%byot%ms%id=bkps%id
	call track_init(bkps,0)
	bkps%state_flag_last=state_ae_evl
	bkps%state_emri_last=0
	bkps%state_emri_current=0
	
	bkps%byot_ini=bkps%byot
	bkps%exit_flag=exit_normal
	bkps%rid=rid
    bkps%source=source_bk 
	bkps%N_gene=1 
	call set_particle_sample_other(bkps)
	
end subroutine
subroutine set_particle_sample_other(bkps)
	use com_main_gw
	implicit none
	type(particle_sample_type)::bkps
	real(8) rnd, period
	 
	bkps%byot%ms%m=bkps%m; bkps%byot%mm%m=spp_new%mbh; bkps%byot%mtot=bkps%m+spp_new%mbh
	bkps%byot%Inc=acos(rnd(-1d0,1d0)); bkps%byot%Om=rnd(0d0,2*pi); 
    bkps%byot%pe=rnd(0d0, 2*pi); 
	bkps%byot%bname='byot'
	bkps%byot_bf%bname='byot_bf'

	bkps%byot%ms%N_gene=bkps%N_gene

	bkps%byot%ms%id=bkps%id
	bkps%byot%mm%obtype=star_type_bh
	bkps%byot%mm%radius=spp_new%mbh/1d8
    bkps%byot%me=rnd(0d0,2*pi)
    bkps%byot%an_in_mode=an_in_mode_mean
    call by_em2st(bkps%byot)
    call by_split_from_rd(bkps%byot)
end subroutine
    

subroutine init_output()
	use com_main_gw
	use md_star_pot
	use md_event_datas
	implicit none
	character*(4) tmpi

	write(unit=tmpi,fmt="(I4)") rid+ctl%ntask_bg+1
	call output_chains_bin(bksams,"output/ini/bin/single/samchn"//trim(adjustl(tmpi)))
	
	ctl%run_snap_time_f=0
	ctl%run_snap_time_i=0
		
	select case(ctl%timestep_mode)
	case(timestep_mode_trh)
		call get_trlx_time_at_rh(ctl%trlx_rh0)
	case(timestep_mode_tnr)
		call get_trlx_min_across_cluster(ctl%trlx_rh0)
	end select
	if(rid.eq.0)then
		print*, "TNR0=",ctl%trlx_rh0
	end if
	 
	if(rid.eq.0)then 
		call output_dms_hdf5_pdf(dms, "./output/ini/hdf5/dms_0") 
		call output_ctl_ini_dstr_hdf5("./output/ini/hdf5/ini") 

		call output_diffuse_mspec_bin("output/ini/hdf5/dms_0")
		write(*,*) "init finished"	 
	end if 
	
end subroutine