module com_main_gw
	use model_basic 
    use md_coeff
    use MPI_comu
    use md_star_pot
    use md_bk_species
    !use md_star_pot
	implicit none      
	real(8)::sample_jlc, sample_jlc_dimless
	real(8)::sample_enf, sample_jmf, sample_mef, sample_af, sample_eni
    real(8)::sample_jf ,sample_alpha
    real(8)::sample_even, sample_evjum
    real(8)::sample_djp,sample_elp,sample_den,sample_djp0
    real(8)::sample_step_nr_e, sample_step_nr_j, sample_step_lc, sample_step_gw
    real(8)::sample_tgw_time, sample_rlx_j_time, sample_rlx_e_time 
    real(8)::sample_table_rdx, sample_table_rdy 

    integer,save::id_saver
    integer::sample_table_idx, sample_table_idy
	integer::sample_mass_idx    

contains
      
    subroutine get_coeff(sample, coeNR, coeRR, coeGW) 
		implicit none
		real(8) drrjj,drrj, ac, ec, sample_mass
		integer i,idx,idy
		type(coeff_type)::coerr, coenr, coegw
        !type(diffuse_coeffient_type),pointer::dc0
		class(particle_sample_type)::sample
	
		!print*, "even,evjum=", even, evjum
        if(ctl%two_body_relaxation_on.ge.1)then
		    call get_coenr(sample_even, sample_evjum, sample%m, sample%en, sample%jc,&
            coenr,sample_table_idx,sample_table_idy)
        else
            coenr%e=0;coenr%ee=0;coenr%j=0;coenr%jj=0;coenr%ej=0
            sample_rlx_e_time=1d30; sample_rlx_j_time=1d30
        end if

		if(ctl%gw_radiation_otby.ge.1 )then
            !if(sample%obtype.eq.star_type_bh)then
                !print*, "sample%en,jm, rp, a(1-e)=", sample%x, sample%jm, sample%rp, sample%byot%a_bin, &
                !    sample%byot%e_bin, sample%byot%a_bin*(1-sample%byot%e_bin)
                if(sample%byot%e_bin<1d0.and.sample%byot%e_bin>0d0)then
                !if(sample%rp<)then
                    !print*, sample%byot%a_bin, sample%byot%e_bin
			        call get_coegw(sample%rp, sample%byot%e_bin, sample%m, sample%period,coegw)
                    if(ctl%chattery.ge.4)then
                        print*, "==coegw bound===================================="
                        print*, "coegw%e,coegw%j=",coegw%e,coegw%j
                        print*, "ac,ec=",sample%byot%a_bin,sample%byot%e_bin
                    end if
                elseif(sample%byot%e_bin>=1d0)then
                    call get_coegw_unbound(sample%rp, sample%byot%e_bin, sample%m, sample%period, coegw)
                    if(ctl%chattery.ge.4)then
                        print*, "==coegw unbound===================================="
                        print*, "coegw%e,coegw%j=",coegw%e,coegw%j
                        print*, "ac,ec=",sample%byot%a_bin,sample%byot%e_bin
                    end if
                else
                    print*, "error! ec<0,ebin,jm=", sample%byot%e_bin, sample%jm
                    stop
                    coegw%e=0; coegw%j=0
                endif
                sample_tgw_time=abs(sample%en/coegw%e)

		end if
		
	end subroutine


	integer function gen_id()
        implicit none
        !print*, "id_saver=",id_saver
		gen_id=id_saver+rid*10000000
        id_saver=id_saver+1
        !print*, "id=", gen_id
	end function
    subroutine create_clone_particle(pt,lvl,amplifier,time)
        !use com_main_gw
        implicit none
        type(chain_pointer_type)::pt
    !	integer,parameter::number_of_clone=9
        type(chain_pointer_type),pointer::ps, pe
        integer lvl, i
        real(8) time
        integer amplifier
        pe=>pt%ed
        call pt%ed%create_chain(amplifier-1) 
        pt%ob%weight_real=pt%ob%weight_real/amplifier
        pt%ob%weight_clone=pt%ob%weight_clone/amplifier
        pt%ob%Lvl_clone=pt%ob%Lvl_clone+1
        ps=>pe%next
        do i=1, amplifier-1
           ! print*, allocated(pe%ob)
            call pt%copy_to(ps)
            ps%ob%create_time=time
            ps%ob%simu_bgtime=time
            !print*, "clone: id_saver=",id_saver
            ps%ob%id=gen_id()
            !print*, "ps%ob%id, pt%ob%id=",ps%ob%id,pt%ob%id
            !read(*,*)
            !ps%ob%source=lvl
            !call ps%ob%get_weight_clone(ctl%clone_scheme, &
            !            amplifier, ctl%clone_e0)
            !ps%ob%weight_clone=pt%ob%weight_clone
            !ps%ob%weight_real=pt%ob%weight_real
      !call particle_sample_get_weight_clone(ps%ob%en, ctl%clone_scheme, &
!                        amplifier,ctl%clone_e0,ps%ob%weight_clone)                        
            ps=>ps%next
        end do
        
    end subroutine
    subroutine show_sample_init(sample, total_time, time)
        !use com_main_gw
        implicit none
        class(particle_sample_type)::sample
        real(8) total_time, time
        integer ntrack_esti

            if(ctl%chattery.ge.3) then
                if(ctl%ntasks.gt.1)then
                    write(unit=chattery_out_unit,fmt="(A25)") "init", star_type_str(sample%obtype), rid
                    write(unit=chattery_out_unit,fmt="(A25, 1PE10.3, I10)") "time=",time/(2*pi)/1d6
                    write(unit=chattery_out_unit,fmt="(A25, 1P3E10.3)") "ac,ec,en=", sample%byot%a_bin,&
                        sample%byot%e_bin,sample%en/ctl%energy0
                else
                    write(*,fmt="(A25)") "init", star_type_str(sample%obtype), rid
                    write(*,fmt="(A25, 1PE10.3, I10)") "time=",time/(2*pi)/1d6
                    write(*,fmt="(A25, 1P3E10.3)") "ac,ec,en=", sample%byot%a_bin,sample%byot%e_bin,sample%en/ctl%energy0 
                end if
                
            end if
            if(sample%exit_flag.ne.exit_normal)then
                print*, "error! sample%exit ne exit_normal"
                print*, "sample%en>ctl%energy_boundary=", sample%en>ctl%energy_boundary
                call sample%print("sample")
                stop
            end if
            !sample%exit_flag=exit_normal
    end subroutine
 
	subroutine if_sample_within_lc(sample)
		implicit none
		class(particle_sample_type)::sample
		!sample%rp=sample%byot%a_bin*(1-sample%byot%e_bin)
        !sample%rp=dms%rp%fxy(sample_table_idx, sample_table_idy)*r0_cl

		if(sample%rp<sample%r_lc)then
			sample%within_jt=1
		else
			sample%within_jt=0
		end if
        if(ctl%chattery.ge.3)then
            print*, "==sample_within_jt=", sample%within_jt
        end if
	end subroutine
	subroutine if_sample_pass_rp(sample,steps, flag)
		implicit none
		class(particle_sample_type)::sample
		integer flag
		real(8) ipdi,ipdf
		real(8) steps, npi,npf

		npi=sample%byot%me
		npf=npi+steps*pi*2
		!print*, "loss cone: sample%me=",npi
		ipdi=npi/2d0/pi-int(npi/2d0/pi)
		ipdf=npf/2d0/pi-int(npf/2d0/pi)
		!print*, "ipdi,ipdf=",ipdi,ipdf
		if(ipdi<0.5.and.ipdf>0.5.or.steps.ge.1d0)then
			flag=1
            if(ctl%chattery.ge.3)then
                print*, "passed pericenter"
            end if
		else
			flag=0
		end if
        if(ctl%chattery.ge.3)then
            print*, "npi, npf=", npi, npf
        end if
	end subroutine
    subroutine get_step(sample,coeNr, coeRR, coeGW, steps, ttot, tnow)
        implicit none
        integer(8) nstemp
        integer num
        real(8) tgw_otby,get_t_gw, r_mean
        real(8)  period, steps
        real(8),external:: rnd
        real(8) NminNr, NminRR, ttot, tnow, ac, na 
        real(8) e2,e1, jnr2, jrr2, jgwr, j1,time_dt_nr
        real(8) nmjl, nmjmr, nmjrr, jm,  nmendrift,time_dt_e, time_dt_j
        real(8) nmende2, nmendee2,deb2, njdrift,net
        real(8) ntmax, nmingwe, nmingwj!,lambda(ctl%m_bins)
        real(8),parameter::avvr=0.31
        real(8) tvrr_step, ratio,wsi
        class(particle_sample_type)::sample
        type(coeff_type)::coerr, coenr, coegw
        period=sample%period

        
        jm=sample%jph
        if(ctl%two_body_relaxation_on.ge.1)then
            select case(ctl%dejmodel)
            case(dejmodel_EJ)
                call get_steps_nr_EJ(sample%en, sample%jm, coenr, sample%jc,  time_dt_e, time_dt_j)
                time_dt_nr=min(time_dt_e, time_dt_j)
                sample_step_nr_e=time_dt_e/period
                sample_step_nr_j=time_dt_j/period
            case(dejmodel_xj)
                call get_steps_nr_xj(sample%en, sample%jm, coenr,  time_dt_nr)
            case default
                print*, "error! define dejmodel", ctl%dejmodel
                stop
            end select
        else
            sample_step_nr_e=1d8
            sample_step_nr_j=1d8
            time_dt_nr=1d8
        end if
        
       	steps=time_dt_nr/period 

		ntmax=(ttot-tnow)/Period
        steps=min(steps,ntmax)

        if(steps.eq.1d6) then
            print*, "warnning steps=", steps
            print*, " ntmax=",    ntmax
            print*, "nmendrift,steps=", nmendrift,steps
            print*, "e1,e2=",e1,e2, coeNr%ee, coeNr%e
            print*, "P,ac=", Period, sample%byot%a_bin, sample%byot%e_bin, sample%en
            print*, "id, rid=",sample%id, rid
            stop
        end if
        
        !print*, "steps=", steps
        if(ieee_is_nan(steps))then
            print*,"1:steps is nan:steps=", steps
            print*, "jnr2, e1, e2, jm=", jnr2, e1, e2, jm
            print*, "sample%id=", sample%id
            print*, "sample%byot%a_bin,ec, period=",sample%byot%a_bin,sample%byot%e_bin, period
            read(*,*)
        end if
        !if(steps>1d99)then
        !    print*, "1: get_steps steps=",steps,ieee_is_finite(steps), &
        !    jrr2, nmjmr, nmjrr, nmjl, jm, jmin
        !    stop
        !end if
        if(sample%r_lc<0)then
            print*, "get_step:error, sample%r_lc<0", sample%r_lc
            print*, "sample%id=",sample%id
            call print_binary(sample%byot)
            select type(sample)
            type is(particle_sample_type)
                print*, "particle"
            
            end select
            stop
        end if

        if(ctl%include_loss_cone.ge.1 .and.&
		sample%en<ctl%energy_boundary)then
			!if(sample%rp<sample%r_lc)then
			!if(sample%within_jt.eq.1)then
			!	steps=min(steps, 1d0)
			!else				
				jnr2=coenr%jj*period
				nmjl=(max(0.05*sample_jlc, 0.25*abs(jm-sample_jlc)))**2/jnr2
				steps=min(steps, nmjl)
                sample_step_lc=nmjl
			!end if
        end if
        if(steps.le.0)then
            print*, "error! steps<0", steps
            print*, "period=",period
            print*, "nmjl,  ntmax=",nmjl,  ntmax
            print*, "ex,jm=", sample_even, sample_evjum
            print*, "coenr%jj,coenr%ee,period=",coenr%jj,coenr%ee,period
            print*, "ttot, tnow=", ttot, tnow
            print*, "sample%jc=",sample%jc,sample%x,sample%jm
            print*, sample_step_nr_e, sample_step_nr_j, sample_step_gw, ntmax
            if(ctl%include_loss_cone.ge.1)then
                print*, "sample_step_lc=", sample_step_lc
            end if
            stop
        end if
 
    
        if(ctl%gw_radiation_otby.ge.1 .and.&
        sample%en<ctl%energy_boundary)then
        	 !if(sample%within_jt.eq.0)then
            call get_step_gw(coegw, sample%en, sample%jph, sample%jm, period,&
                sample_step_gw)
			steps=min(steps, sample_step_gw)
            !else
			!	steps=min(steps,1d0)
			!end if
            
        !	print*, steps
        endif
    !	if(sample%within_kl)then  ! forget the reason for this??
    !		steps=min(steps, 100d0)
    !	end if
    
		
        if(ieee_is_nan(steps).or.steps<0)then
            print*,"2:steps is nan:steps,  nmjl =",  steps,  nmjl
            print*, "jnr2,jmin,jmax,rtd/ac=", jnr2, sample_jlc,sample%jc, sample%r_lc/sample%byot%a_bin
            print*, "coeNr%jj,Period=", coeNr%jj,Period
            print*, "sample%jm, jph=",sample%jm,sample%jph
             
            print*, "sample%byot%a_bin,ec, period=",sample%byot%a_bin,sample%byot%e_bin, period 
            read(*,*)
        end if
        if(ctl%chattery.ge.4)then
			print*, "====get steps=============================="
			print*, "steps, steps_nr_e, step_nr_j, steps_gw, step_ntmax, sample_step_lc=",steps, &
            sample_step_nr_e, sample_step_nr_j, sample_step_gw, ntmax,sample_step_lc
        end if
    end subroutine
	
	subroutine get_sample_r_td(sp)
		implicit none
		class(particle_sample_type)::sp
		real(8) wsi, ac,rtd,jphlc2,jc_dm

        !---first get the traditional tidal radius
        
		select type(sp)
		type is(particle_sample_type)
			select case(sp%obtype)
			case(star_type_ms,star_type_rg, star_type_bd,star_type_nakedHe)
				if(sp%byot%ms%radius.eq.0)then
					print*, "ms radius=0?? check",sp%obtype
                    call sp%sh%print()
					stop
				end if
				rtd=(spp_new%mbh/sp%m*ctl%eta_td_star)**(1/3d0)*sp%byot%ms%radius
			case(star_type_bh,star_type_ns,star_type_wd,star_type_dark_matter)				
                rtd=0d0
			case default
				print*, "error! star type not defined",sp%obtype
				stop
			end select
		 
		class default
			print*, "?? which type??"
			stop
		end select

        call get_riso_compact_obj_given_spin_inc(mbh_spin,sp%byot%inc,mbh_iso_radius)
        sp%r_lc=max(mbh_iso_radius, rtd)        

	end subroutine
     
    subroutine set_jm_init(bkps)
		implicit none
		class(particle_sample_type)::bkps
		real(8) rtd, jlc,tmp
		real(8),external::fpowerlaw_rnd, rnd
		
		select case(ctl%boundary_fj)
		case(boundary_fj_iso)
			bkps%jm=fpowerlaw_rnd(1d0,0.0044d0,0.99999d0)  
		case default
			print*, "define flag INI", ctl%boundary_fj
			stop
		end select
	end subroutine
    
    subroutine get_type_idx(sample, idxob)
        implicit none
        class(particle_sample_type)::sample
        integer idxob
        select type(sample)
            type is(particle_sample_type)
                select case(sample%obtype)
                case(star_type_ms)
                    idxob=1
                case(star_type_bh)
                    idxob=2
				case(star_type_ns)
					idxob=3
				case(star_type_wd)
					idxob=4
				case(star_type_bd)
					idxob=5
                case default
                    print*, "error in get_type_idx:sg",sample%id
                    stop
                end select
             
        end select
    end subroutine
    subroutine update_track(sp,j)
        implicit none
        class(particle_sample_type)::sp
        integer(8) j
        
        if((sp%write_down_track.ge.record_track_nes&
        .or.ctl%trace_all_sample.ge.record_track_all))then
            !print*, "sp%track_step=",sp%track_step
            !read(*,*)
            if(mod(j,sp%track_step).eq.0)then
                if(sp%length_to_expand>MAX_LENGTH)then
                    sp%track_step=sp%track_step*10 
                    call track_compress(sp, 10)
                    print*,"track compressed"
                end if
            end if
    end if
    end subroutine
    subroutine track_compress(sp, ns)
        use model_basic
        implicit none
        class(particle_sample_type)::sp
        type(track_type),allocatable::tk(:)
        integer ns,j,i
    
        allocate(tk(sp%length))
        tk=sp%track
        j=0
        do i=1, sp%length, ns
            j=j+1
            sp%track(j)=tk(i)
        end do
        sp%length=j
    end subroutine 
	
    subroutine get_sample_weight_real(sp)
        implicit none
        class(particle_sample_type)::sp
        sp%weight_real=sp%weight_clone*sp%weight_n*ctl%n_basic
    end subroutine
 

    subroutine output_sample_track_txt(sp,fl)
        class(particle_sample_type)::sp
        character*(*) fl
        select type(ca=>sp)
        type is (particle_sample_type)
            call output_sg_sample_track_txt(ca,fl)
         
        end select
    end subroutine
    subroutine add_track(t,sp, state_flag)
        !use model_basic
        implicit none
        class(particle_sample_type) ::sp
        type(track_type),allocatable:: tk(:)
        real(8) t
        integer i, state_flag
        i=sp%length
       ! print*, "i=",i, sp%length_to_expand
        if(i.eq.sp%length_to_expand)then
          if(ctl%chattery.ge.4) then
              print*, "track expanding"
          end if
          sp%length_to_expand=sp%length_to_expand+track_length_expand_block
          allocate(tk(i))
          tk(1:i)=sp%track(1:i)
          if(allocated(sp%track)) deallocate(sp%track)
          allocate(sp%track(sp%length_to_expand))
          sp%track(1:i)=tk(1:i)
        end if
        i=i+1
        sp%length=i
        sp%track(i)%time=t
        sp%track(i)%r_lc=sp%r_lc
        sp%track(i)%meout=sp%byot%me
        sp%track(i)%ac=sp%byot%a_bin
        sp%track(i)%ec=sp%byot%e_bin
        sp%track(i)%rp=sp%rp
        sp%track(i)%jm=sp%jm
        sp%track(i)%x=sp%x
        sp%track(i)%incout=sp%byot%inc
        sp%track(i)%omout=sp%byot%om
        sp%track(i)%state_flag=state_flag


        select type(ca=>sp)
		type is(particle_sample_type)
			!do i=1, sp%length
			!	print*, sp%track(i)%time, sp%track(i)%state_flag
			!end do
			!read(*,*)
         
        end select

		
     !	print*, sp%length
     end subroutine

	 subroutine get_move_result(sample,den,djp,steps, &
			enf, jf, mef, af)
		implicit none
		class(particle_sample_type)::sample
		real(8) den, djp,steps
		real(8) ai,ei,Eni, Enf, Ji, Jf, af
		real(8) jmf, mef

		!sample%byot%a_bin=-mbh/(2*sample%en)
		!ai=sample%byot%a_bin
		!ei=sample%byot%e_bin		
		select case(ctl%dejmodel)
		case(dejmodel_EJ)
			Eni=sample%en
			Ji=sample%jm*sample%jc
			Enf=Eni+den
			jf=Ji+djp

		case(dejmodel_xj)
			Eni=sample%en/ctl%energy0
			Enf=(Eni+den)*ctl%energy0
			Ji=sample%jm
			Jf=ji+djp
			af=-spp_new%mbh/(2*Enf)
			!jmf=jf
		
			
		case default
			print*, "error! define dejmodel", ctl%dejmodel
			stop
		end select
		!print*, "af,jf=",af,jf
		!stop
		mef=sample%byot%me+(steps-int(steps))*2*pi
		!print*, "me=", sample%byot%me

        
        if(ctl%chattery.ge.4)then
            print*, "==========get_move_result=================="
            print*, "Eni, den, enf, xf=", Eni, den, enf, enf/ctl%energy0
            print*, "ji,  djp, jf=", ji, djp, jf
        end if
	 end subroutine
     subroutine set_ebound_for_samples(sample)
        implicit none
		class(particle_sample_type)::sample
        !type(star_pot_para)::spp

        ! the maximum energy can not be crossed, if there is no black hole in the center.
        if(sample%en<ctl%energy_max)then
            !print*, "sample%x,jm:i=", sample%en/ctl%energy0, sample%jm, sample%jph
            !sample%en=sample%byot_bf%e
            !sample%en=2*ctl%energy_max-sample%en
            !print*, "sample%en,emax=",sample%en,ctl%energy_max
            sample%en=ctl%energy_max*0.9999999
            running_correction_emax=running_correction_emax+1
            !call get_ex_idx(sample%en/ctl%energy0, idx,even)
            !sample%jph=sample%jm*dms%jc%fx(idx)*ctl%v0*r0_cl
            !print*, "sample%x,jm:f=", sample%en/ctl%energy0, sample%jm, sample%jph
            !read(*,*)
        end if
     end subroutine
	 

end module 
 
subroutine set_star_spin_random(pr)
    use com_main_gw
    implicit none
    type(particle)::pr
    real(8) rnd
    real(8) theta ,phi
    real(8) pd
    select case(pr%obtype)
    case(star_type_ms)
        theta=rnd(0d0,pi)
        phi=rnd(0d0,2*pi)
        ! assuming rotation period of 30 days
        pd=30d0/(365.2425/2d0/pi)
        pr%spin=(/cos(theta)*sin(phi), cos(theta)*cos(phi), sin(theta)/)*2*pi/pd
    case default
        pr%spin=0d0
    end select
end subroutine
    