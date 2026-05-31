subroutine run_one_sample(pt, run_time)
	use com_main_gw
	use md_coeff
	implicit none
	type(chain_pointer_type),target::pt
	real(8) elp, en0,en1,  steps, run_time, GET_T_GW
    real(8) tgw, period, npi,npf
	real(8) rp_out, rtmax, total_time, time, time_next,time_dt, time_create
	integer out_flag_boundary
	character*(100) str_flag
	real(8) ipdi,ipdf, wsi
    integer flag, flag_dedj,out_flag_clone,ier
    integer,parameter::flag_sg=1,flag_by=2
	integer(8),parameter:: ouput_freq=1,pause_freq=1
	type(coeff_type)::coeNr, coeRR, coeGW
    logical:: init_out, particle_cloned
	integer,save::num_out=0
	interface 
		subroutine run_one_sample_inside_cluster(pt, time, total_time)
			use com_main_gw
			implicit none
			type(chain_pointer_type),target::pt
			real(8) time, total_time
		end subroutine
	end interface
	
	associate(sample=>pt%ob)
		time=sample%simu_bgtime*1d6*2*pi
		total_time=run_time*1d6*2*pi

		if(sample%weight_real.eq.0d0.or.sample%weight_n.eq.0d0)then
			print*, "ac_ec_evl_single:error, sample%weight_real or n=0"
			call sample%print("ac_ec_evl single")
			stop
		end if
		call show_sample_init(sample, total_time, time)
 
		out_flag_boundary=0

		!print*, "sample%en, ctl%energy_boundary=", sample%en, ctl%energy_boundary
		
100  	if(sample%en>ctl%energy_boundary)then
			print*, "run_bd: time, total_time=",time/2d6/pi, total_time/2d6/pi, sample%id
			print*, "??", sample%en, ctl%energy_boundary
			print*, "sample%x, x2=",sample%x,sample%en/ctl%energy0, ctl%x_boundary
			call sample%print("run_one_sample")
			stop 
		else
			select type(ca=>sample)
			type is(particle_sample_type)
				call run_one_sample_particle_inside_cluster(pt,time, total_time)
			end select
			if(pt%ob%exit_flag.eq.exit_boundary_min)then
				ctl%num_boundary_elim=ctl%num_boundary_elim+1
			end if
		end if 
		if(ctl%chattery.ge.1)then 
			if(ctl%chattery.eq.1)then
				if(sample%exit_flag.ne.exit_boundary_min)then
					select type(sample)
					type is(particle_sample_type)
						call print_results_single(pt,  pt%idx, pt%ed%idx, sample)
					end select
				end if
			else
				select type(sample)
				type is(particle_sample_type)
					call print_results_single(pt,  pt%idx, pt%ed%idx, sample)
				end select
			end if
			if(ctl%chattery.ge.4.or.ctl%debug.ge.1)then
				if(ctl%chattery.eq.4)then
					if(sample%exit_flag.ne.exit_normal.and.sample%exit_flag.ne.exit_invtransit)then
						call get_exit_flag_str(sample%exit_flag, str_flag)
						print*, "exit=", trim(adjustl(str_flag))
						read(*,*)
					end if
				else
					read(*,*)
				end if
			end if

		end if
	end associate

end subroutine

subroutine run_one_sample_particle_inside_cluster(pt, time, total_time)
	use com_main_gw
	use md_coeff
	use md_stellar_evolution
	implicit none
	type(chain_pointer_type)::pt
	type(particle_sample_type)::sp_test
	!type(particle_sample_type),pointer::sample
	real(8) elp, en0,en1,jm0, jm1,  steps, run_time
	real(8),external:: GET_T_GW
    real(8) tgw, period
	real(8) total_time, time, time_next,time_dt, time_create
	integer(8) j
	real(8) rnd,wsi
    integer flag, out_flag_clone,flag_pass_rp,ier, flag_stellar_evl
    integer,parameter::flag_sg=1,flag_by=2
	integer(8),parameter:: ouput_freq=1,pause_freq=1
	type(coeff_type)::coeNr, coeRR, coeGW
	integer,save::num_out=0
	real(8) dt_block, tp, jlc,ratio,c0,mass_change
	logical::tidal_condition
 
    associate(sample=>pt%ob)
		!if(sample%en>ctl%energy_boundary.and.sample%jm<0.04d0) then
		select type(sample)
			type is (particle_sample_type)
				j=0
				!print*, "run_one_sample_particle_inside_cluster"
				
				! if(sample%m>20)then
				!  	ctl%chattery=5
				! else
				!  	ctl%chattery=0
				! end if
				! if(rid.eq.8.and.time/1d6/2d0/pi>200.7)then
				! 	ctl%chattery=2
				! end if

				if(time.ge.total_time) time_next=time
				
				call get_mass_idx(sample%m,sample_mass_idx)
				
				if(sample_mass_idx==-1)then
					print*, "sample source, weight_real=", sample%source, sample%weight_real
				end if
				if(sample%weight_real<1d-30)then
					print*, "error!, sample%weight_real=",sample%weight_real, sample%weight_n, sample%weight_clone,sample%Lvl_clone
					stop
				end if 
				loop1:do while(time<total_time)
 
					ctl%num_of_loops=ctl%num_of_loops+1

					if(ctl%include_stellar_evolution.ge.1)then
						call get_current_particle_stellar_info(sample,time/1d6/2d0/pi, flag_stellar_evl,ctl%chattery.ge.3) 
						if(flag_stellar_evl.eq.-1)then
							sample%exit_flag=exit_stellar_evl_supnov
							exit loop1
						end if
					end if

					call update_sample_para(sample,spp_new)
					 
					if(ctl%include_loss_cone.ge.1 )then
						call get_sample_r_td(sample)
						if(sample%jc.le.0d0)then
							print*, "error! sp%jc<=0", sample%jc
							stop
						end if
						call get_sample_jlc(sample%x,spp_new%mbh_dmless,sample%r_lc/r0_cl,sample%jc/(ctl%v0*r0_cl),spp_new,&
							sample_jlc_dimless,ier)
						if(ier.eq.1)then 
							print*, "rtd, sp%x,  sp%r_lc, star_type_str(sp%obtype)"
							print*, sample%r_lc/r0_cl,  sample%x, sample%r_lc, star_type_str(sample%obtype),&
								sample%byot%ms%radius
							call if_sample_pass_rp(sample, steps,flag_pass_rp)
							print*, "flag_pass_rp=",flag_pass_rp
							print*, "mbh_dmless,mbh=",spp_new%mbh_dmless,spp_new%mbh 
							print*, "sample%rp,ra=",sample%rp,sample%ra
							!stop
							sample%exit_flag=exit_lc
							if(ctl%trace_all_sample.ge.record_track_nes.or.&
							sample%write_down_track.ge.record_track_detail)then
								call add_track(time/2d6/pi,sample,state_plunge)
							end if
							exit loop1
						end if
						sample_jlc=sample_jlc_dimless*sample%jc
					end if
					call get_coeff(sample,coeNr, coeRR, coeGW) 

					if(ctl%include_loss_cone.ge.1)then
						call if_sample_within_lc(sample)
					end if
					call get_step(sample,coeNr, coeRR, coeGW,steps, total_time, time)
					if(ctl%include_stellar_evolution.ge.1)then
						call get_step_stellar_evl(sample,steps,time,ctl%chattery.ge.3)
					end if
 
					en0=sample%en
					jm0=sample%jm 
					if(steps>1d99)then
						print*, "af get_steps steps=",steps,ieee_is_finite(steps)
						call sample%print("sample")
						print*, "sample%x,jm=",sample%x,sample%jm
						print*, "keplerian a,e=",sample%byot%a_bin,sample%byot%e_bin
						print*, "sample%period=",sample%period
						print*, "sample_jlc_dimless=",sample_jlc_dimless
						print*, "sample%rp, sample%ra=", sample%rp, sample%ra
						print*, "sample%jc=",sample%jc
						stop
					end if

					!period=P(sample%byot%a_bin)
					period=sample%period
					time_dt=steps*period 
					call get_de_dj(sample, coeNR, coeRR,coeGW, time,time_dt,steps, period) 
					
					call get_move_result(sample,sample_den,sample_djp,steps,&
						sample_enf, sample_jf, sample_mef, sample_af) 
					if(ctl%chattery.ge.5)then
						print*, "======================================================================="
						print*, "sample_rlx_e_time,sample_rlx_j_time,sample_tgw_time=",sample_rlx_e_time,sample_rlx_j_time, sample_tgw_time, &
							sample%state_emri_current,sample%state_emri_last
						print*, "======================================================================="
					end if
					
					if(ctl%include_loss_cone.ge.1 .and.&
						sample%en<ctl%energy_boundary)then
						call if_sample_pass_rp(sample, steps,flag_pass_rp)
						if(ctl%chattery.ge.3)then
							print*, "flag_pass_rp=",flag_pass_rp
						endif
						if(flag_pass_rp.ge.1.or.sample%ra<sample%r_lc)then
							if(ctl%chattery.ge.3)then
								print*, "within_jt=",sample%within_jt
							end if
							if(sample%within_jt.eq.1)then	
								 
								if(ctl%chattery.ge.3)then
									print*, "========================exiting=================="
								end if
								if(ctl%gw_radiation_otby.ge.1)then
									if(sample%state_emri_current.ge.1.and.sample%state_emri_last.ge.1)then
										sample%exit_flag=exit_emri_single 
										if(ctl%trace_all_sample.ge.record_track_detail)then
											call add_track(time/1d6/(2*pi), sample,state_emri)
										end if
										exit loop1 
									end if 
								end if

								tidal_condition=.false.
								select case(sample%obtype)
									case(star_type_ms,star_type_rg,star_type_nakedHe, star_type_bd)
										if(sample%rp>2*mbh_radius)then
											tidal_condition=.true.
										end if
									case default
								end select
								
								if(tidal_condition)then
									if(abs(sample_djp0)<sample_jlc_dimless)then
										sample%exit_flag=exit_tidal_empty
									else
										sample%exit_flag=exit_tidal_full
									end if
									if(ctl%trace_all_sample.ge.record_track_nes.or.&
									sample%write_down_track.ge.record_track_detail)then
										call add_track(time/1d6/(2*pi), sample,state_td)
									end if
								else   
									sample%exit_flag=exit_lc
									if(ctl%trace_all_sample.ge.record_track_nes.or.&
									sample%write_down_track.ge.record_track_detail)then
										call add_track(time/2d6/pi,sample,state_plunge)
									end if
									ctl%num_of_lc=ctl%num_of_lc+1
								end if							
								exit loop1
								 
							else
								!! test if the particle dominate by GW radiation EACH TIME IT PASSED PARICENTER
								if(ctl%gw_radiation_otby.ge.1)then
								sample%state_emri_last=sample%state_emri_current
								if(sample_tgw_time<ctl%emri_cri*min(sample_rlx_e_time,sample_rlx_j_time))then
									sample%state_emri_current=1										
								else
									sample%state_emri_current=0
								end if
							end if
							end if		
						end if				
					end if

					call update_track(sample, j)
					!print*, ctl%trace_all_sample.ge.record_track_detail
					if(sample%write_down_track.ge.record_track_detail&
						.or.ctl%trace_all_sample.ge.record_track_detail)then
						call add_track(time/1d6/(2*pi), sample,state_ae_evl)
					end if

					j=j+1 
					if(j>MAX_RUN_LENGTH)then
						sample%exit_flag=exit_max_reach
						print*, "j=",j
					exit loop1
					end if
					if((j.gt.MAX_RUN_LENGTH/10.and.j.lt.MAX_RUN_LENGTH/10+10).or.j.eq.MAX_RUN_LENGTH/2)then
						print*, "single:warning, j, rid,sp%id=", j, rid, sample%id
						print*, "step,ao,eo=", steps
						print*, "within_jt=",sample%within_jt
						print*, "sample%obtype=",star_type_str(sample%obtype)
						 
						print*, "sample_jlc, jlc_dmless=",sample_jlc, sample_jlc/sample%jc
						print*, "sample%x,jm,jc,jcdm=",sample%x, sample%jm, sample%jc, sample%jc/(ctl%v0*r0_cl)
						print*, "time, period=", time, period
						print*, "sample_step_nr_e, nr_j, lc, gw=", sample_step_nr_e, sample_step_nr_j, &
							sample_step_lc, sample_step_gw
						print*, "coenr%ee,coenr%jj=",coenr%ee, coenr%jj, coenr%j
						print*, sample_table_idx,sample_table_idy
						print*, "sample_den,sample_djp=",sample_den,sample_djp
						print*, "sample_enf, sample_jf=",sample_enf, sample_jf
						print*, "sample_xf,rp,r_td=",sample_enf/ctl%energy0,sample%rp,sample%r_lc
						print*, "sample_mef=",sample_mef
						!stop
					end if
 

					if(steps>1d99.or.ieee_is_nan(steps).or.steps<0)then
						print*, "bf get_dedj steps=",steps,ieee_is_finite(steps)
						call sample%print("ac_ec_evl_single")
						print*, "time, rid=", time, rid
						stop
					end if 
					
					call move_de_dj_one(spp_new,sample,sample_eni,sample_enf, sample_jf, sample_mef,sample_af)
					 
					en1=sample%en
					 
					time_next=time+time_dt 
					if(en1<ctl%energy_max.and.spp_new%mbh_dmless.ne.0)then 
						sample%exit_flag=exit_boundary_max
						boundary_sts_emax_cros=boundary_sts_emax_cros+1
						 
						if(ctl%trace_all_sample.ge.record_track_nes)then
							sample%x=sample%en/ctl%energy0 
							call add_track(time_next/2d6/pi,sample,state_emax)
						end if
						exit loop1
					end if
					
					if(ctl%clone_scheme.ge.1)then
 
						if(.not.(ctl%gw_radiation_otby.ge.1.and. (sample%state_emri_last.eq.1.and. sample%state_emri_current.ge.1)))then
							if(sample_mass_idx.ne.-1)then
								call clone_scheme(pt, en0, en1, ctl%clone_factor(sample_mass_idx),&
									time_next/1d6/2d0/pi, out_flag_clone) 
								if(out_flag_clone.eq.100)then
									sample%exit_flag=exit_invtransit
									ctl%num_clone_elim=ctl%num_clone_elim+1
									exit loop1
								end if
							end if
						end if				
					end if 

					if(en1>ctl%energy_boundary)then 
						sample%exit_flag=exit_boundary_min
						ctl%mass_move_out_of_emin=ctl%mass_move_out_of_emin+sample%weight_real*sample%m
						exit loop1
					end if
					
					time=time_next
					if(ctl%chattery.ge.3)then
						print*, "=================================end"
						print*, "current time, time_tot=",time/1d6/2d0/pi, total_time/1d6/2d0/pi
						print*, "===================================="
					end if 
				end do loop1
				if((sample%exit_flag.ne.exit_boundary_min))then
 
					if(spp_new%mbh_dmless.ne.0)then
						call update_sample_para(sample,spp_new)
					end if
 
				else
					sample%x=sample%en/ctl%energy0
				end if
				if(ctl%chattery.ge.3) then
					print*, "exit time=", time/1d6/2d0/pi, sample%exit_flag
					if(sample%obtype.eq.star_type_bh.and.sample%exit_flag.eq.exit_boundary_max)then
						read(*,*)
					end if
				end if
				sample%exit_time=time_next/1d6/2d0/pi

				if(ctl%chattery.ge.5) then
						write(*,*) "finished, rid=",rid
						write(*,fmt="(A25, 1PE11.4)") "-------time exit:",sample%exit_time
						write(*,fmt="(A25, I7)") "-------exit flag:",sample%exit_flag
						print*, "Enter to go"	
						read(*,*)
				end if
		end select
		
	end associate
end subroutine

subroutine print_results_single(pt,  id, eid,sample)
	use com_main_gw
	implicit none
	type(chain_pointer_type):: pt
	type(particle_sample_type)::sample
	integer id, eid,  spid
	character*(3) str_type

	!call get_binary_types(sample%bytype,str_type)
	spid=sample%id
	call get_star_type(sample%obtype,str_type)
	select case(sample%exit_flag)
		case(exit_normal) 
			if(ctl%chattery.ge.2)then
				write(chattery_out_unit,fmt="(A25, A5, 5I8, 2I10)") "-------time out",&
                str_type,id, sample%source,rid,eid, sample%n_gene,spid, pt%ed%ob%id 
			end if 
		case(exit_tidal_empty)
			write(chattery_out_unit,fmt="(A25, A5, 5I8, 2I10)") "-------td empty",&
                str_type,id, sample%source,rid,eid, sample%n_gene, spid, pt%ed%ob%id
		case(exit_tidal_full)
			write(chattery_out_unit,fmt="(A25, A5, 5I8, 2I10)") "-------td full",&
                str_type,id, sample%source,rid,eid, sample%n_gene, spid, pt%ed%ob%id
		case(exit_max_reach)
			write(chattery_out_unit,fmt="(A25, A5, 5I8, 2I10)") "-------MORE STEPS NEEDED",str_type,id, sample%source,rid,eid,&
                     sample%n_gene, spid, pt%ed%ob%id
			write(chattery_out_unit,fmt="(A25, 1P4E10.3)") "--ac,ec,ain,ein=",&
                sample%byot%a_bin,sample%byot%e_bin
		case (exit_lc)
			write(chattery_out_unit,fmt="(A25, A5, 5I8, 2I10)") "----emri_plunge",&
                str_type,id, sample%source,rid,eid, sample%n_gene, spid, pt%ed%ob%id
		case (exit_emri_single)
			write(chattery_out_unit,fmt="(A25, A5, 5I8, 2I10)") "----emri_single",str_type,id, &
                sample%source,rid,eid, sample%n_gene, spid, pt%ed%ob%id
		case(exit_boundary_min)
			write(chattery_out_unit,fmt="(A25, A5, 5I8, 2I10)") "-------exit_to_emin",str_type,id, &
                sample%source,rid,eid, sample%n_gene, spid, pt%ed%ob%id
		case(exit_boundary_max)
			write(chattery_out_unit,fmt="(A25, A5, 5I8, 2I10)") "-------exit_to_emax",str_type,id, &
                sample%source,rid,eid, sample%n_gene, spid, pt%ed%ob%id
 
		case(exit_invtransit)
		case(exit_other)
			write(chattery_out_unit,fmt="(A25, A5, 5I8, 2I10)") "-------other",str_type,id, &
                sample%source,rid,eid, spid, sample%n_gene, pt%ed%ob%id 
		case(exit_stellar_evl_supnov)
			write(chattery_out_unit,fmt="(A25, A5, 5I8, 2I10)") "-------supernovea",str_type,id, &
                sample%source,rid,eid, sample%n_gene, spid, pt%ed%ob%id
		case default
			write(chattery_out_unit,*) "define state:", sample%exit_flag
			!stop
	end select
end subroutine
