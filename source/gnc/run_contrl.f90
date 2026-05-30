
subroutine run_one_snap(cur_time_i, cur_time_f,  smsa, n,update_dms)
    use com_main_gw
    implicit none
    real(8) cur_time_f, cur_time_i!, cur_time_f_istep
    integer n, ierr, istep
	type(particle_samples_arr_type) smsa(n)
    integer,save::norm=5
    real(8) t1, t2
	logical::update_dms
    
    call show_memory_usage()
    if(rid.eq.0) call get_collection_memory_usage(smsa, n)
    call show_total_memory_usage()
    !   call check_rp_ra("run_one_snap: start")
    if(n.eq.0.or.size(smsa).eq.0)then
        print*, "error! smsa not initialized!", n
        stop
    endif
    !print*,"size of smsa:", size(smsa),rid
    if(rid.eq.0)then
		call cpu_time(t1)
	end if
    
    call update_arrays_single(.true.)
    ctl%run_snap_time_inner_steps_i=cur_time_i
    ctl%run_snap_time_inner_steps_f=cur_time_i
    do istep=1, ctl%num_step_per_update
        ctl%run_snap_time_inner_steps_i=ctl%run_snap_time_inner_steps_f
        ctl%run_snap_time_inner_steps_f=cur_time_i+&
        (cur_time_f-cur_time_i)/dble(ctl%num_step_per_update)*dble(istep)
        
        
        write(unit=chattery_out_unit,fmt=*)  "sg:istep,rid, cur_time_f_istep=",&
            istep,rid, ctl%run_snap_time_inner_steps_f
        !call check_rp_ra("before RR_mpi")
 

        call mpi_BARRIER(mpi_comm_world, ierr)
        !call check_boundary("bg")
        call RR_mpi(bksams, ctl%run_snap_time_inner_steps_f)
        !call check_boundary("ed")
        call mpi_BARRIER(mpi_comm_world, ierr) 

        if(rid.eq.0)then
            print*, "reset time"
        end if
        call reset_create_time(bksams)
        !call init_chattery()

        !if(ctl%barge_evl_method.eq.barge_evl_method_direct)then
        !    call get_sample_para(dms,bksams_arr_norm,.false.,spp_new)
        !end if
        if(ctl%trace_all_sample.gt.0.or.ctl%insnapmode.eq.snap_mode_one)then
            call update_arrays_single(.false.)
            cycle
        else
            call update_arrays_single(.true.)
        end if 
        call get_system_dstr_adb_cor_one_time() 
        call mpi_barrier(mpi_comm_world,ierr) 

        call update_sample_energy_indvd(ctl%replace_sample_eceed_emax)
        print*, "stpt:update energy finished", rid        
        call mpi_barrier(mpi_comm_world,ierr)

    end do
    
    if(rid.eq.0)then
        call print_norm_dms(dms,chatterY_out_unit) 
    end if

    ! if(ctl%gw_radiation_inby.ge.1)then
    !     write(unit=chattery_out_unit,fmt=*)  "rid, n_colld,n_colld_bh,n_colld_bh_2g, nexchange_tot, nexchange_2gene=",rid,&
    !     ctl%num_collide_tot, ctl%num_bh_collide_tot, ctl%num_bh_2genecollide_tot, ctl%num_exchange_tot, &
    !     ctl%num_exchange_2gene
    ! end if 

    if(rid.eq.0)then
		call cpu_time(t2)
        write(unit=chattery_out_unit,fmt=*)  "one snap running time:", t2-t1
	end if
     
    if(ctl%two_body_relaxation_on.ge.1)then
        if(update_dms)then
            call init_dms_dc(dms)   
            call get_dms(dms) 

        end if
    end if
    if(rid.eq.0)then 
        call show_memory_usage()
    end if 
    call mpi_barrier(mpi_comm_world, ierr) 

    call check_weightings("end of snapshot")
end subroutine
subroutine check_weightings(str_)
    use com_main_gw
    implicit none
    character*(*) str_
    type(chain_pointer_type),pointer::ps
    real(8) weight_c 
    ps=>bksams%head
    do while(associated(ps))
        select type(ca=>ps%ob)
        type is (particle_sample_type)
            weight_c=ca%weight_clone*ca%weight_n*ctl%n_basic
            if(abs(weight_c-ca%weight_real)>0.1)then
                print*, "warnning! weighting not matched!"
                print*, "weight_correct,weight_now,weight_clone,weight_n,exit_flag,source=",&
                    weight_c,ca%weight_real, ca%weight_clone, ca%weight_n,ca%exit_flag,ca%source
                stop
            end if
        end select
        ps=>ps%next
    end do
	if(rid.eq.0)then
    	print*, str_, " check finished!"
	end if

end subroutine
  
 subroutine get_v_dispersion_one_ir(gx, phi_r, rho_r, vdisp)
    use com_main_gw
    implicit none
    integer n, i
    type(s1d_ird_type)::gx
    type(s1d_ird_type)::common_gx_ir
    real(8),intent(in)::phi_r, rho_r
    real(8),intent(out)::vdisp
    real(8) int_out
    integer idid
    common_gx_ir=gx
    !call get_none_zero_s1d_ir(gx,common_gx_ir)

    int_out=0
    call my_integral_acc(0d0, (phi_r*2)**0.5,int_out,1d-20,1d-12,fcn, idid)
    vdisp=(int_out*4*pi/rho_r*(2*pi)**(-1.5d0)*ctl%v0**2/3d0)**0.5  ! vdisp^2=<v^2>/3
contains
	subroutine fcn(n, x, y, f, par, ipar)
		use, intrinsic :: ieee_arithmetic
		implicit none
		integer n, ipar(100)
		real(8) x, y(n), f(n), par(100)
		real(8) logx, barge_tmp
        if(phi_r-0.5d0*x**2>10**common_gx_ir%xmin)then
            logx=log10(phi_r-0.5d0*x**2)
            call common_gx_ir%get_value_l(logx,barge_tmp)
            if(barge_tmp<0) barge_tmp=0
        else
            barge_tmp=0d0
        end if
		f(1)=barge_tmp*x**4
        !print*, "f(1), barge, phi_r-0.5*x**2=", f(1), barge_tmp, phi_r-0.5*x**2
	end subroutine    
 end subroutine
 
 subroutine get_trlx_time_at_rh(trlx_rh)
    use com_main_gw
    use md_star_pot
    implicit none
    !type(s1d_type)::s1d_trlx
    real(8) trlx_rh,phi_tmp,rho_tmp
    real(8) rh_dmless,rh,vh,nm2_tot,lglambda, rho_tot,rh_now
    integer k
    if(spp_new%mbh.eq.0)then
        call get_lambda(lglambda)
    else
        !call get_lambda(lglambda)
        lglambda=log(spp_new%mbh*0.5)
    end if
    call get_rh_now(rh_now)
    !call get_rh_in_pc(mbh,rh_now)
    rh=10**rh_now*r0_cl/pc
    rh_dmless=rh*pc/r0_cl
    nm2_tot=0
    rho_tot=0
    do k=1, dms%n
        !!=================
        !print*, "here, spt_rho_rmin should be corrected"
        call get_rho_full_range(dms%mb(k)%all%fmden,dms%mb(k)%all%spt_rho_rmin,&
        log10(rh_dmless),rho_tmp,spp_new) 
        !!=================
        nm2_tot=nm2_tot+RHO_TMP*ctl%n0*dms%mb(k)%mc
        rho_tot=rho_tot+rho_tmp
        !if(rid.eq.0)then
        !    print*, "k, ctl%n0, rho_tmp=",k,ctl%n0, rho_tmp,dms%mb(k)%mc
        !end if
    end do

    call get_phi_star_full_range(spp_new,log10(rh_dmless),phi_tmp)
    phi_tmp=10**phi_tmp+spp_new%mbh_dmless/rh_dmless
    !select case(ctl%barge_grid_type)
    !case(barge_grid_type_iregular_barp,barge_grid_type_iregular_phi,barge_grid_type_regular)

        call get_v_dispersion_one_ir(dms%all%all%barge_ir,phi_tmp,rho_tot,vh)
 

    if(nm2_tot.ne.0.and.vh.ne.0)then
        call get_two_body_relaxation_time(vh, nm2_tot,lglambda,trlx_rh)
 
    else
        trlx_rh=1d99
    end if
    if(rid.eq.0)then
        print*, "trlx_rh,rho_tot,vh,rh=", trlx_rh,rho_tot,vh,rh_dmless
    end if

 end subroutine
 subroutine get_trlx_time_across_cluster(s1d_trlx)
    use com_main_gw
    use md_star_pot
    ! use md_lagrange
    implicit none
    !type(s1d_type)::fmden
    real(8) trlx, lglambda, nm2_tot, vh!, rhalf_mass
    real(8) thf, mtot, r0, rhoin, rmin,  rho_tmp, phi_tmp
    real(8) rcmin,r_c_iter
    type(s1d_type)::s1d_trlx
    integer i,k,ier
    
    call get_lambda(lglambda)
    select case(ctl%fden_ana_est_method)
    case(fden_ana_est_method_1d_iso)
        s1d_trlx=spp_new%fphi_star
    case(fden_ana_est_method_2d)
 
        call dms%rc%get_value_l(sample_logemax,rcmin)
        !print*, "rcmin, sample_logemax=",rcmin,sample_logemax
        call s1d_trlx%init(log10(rcmin),sample_logrmax,ctl%dstr_bins_r,sts_type_dstr)
        call s1d_trlx%set_range()
    end select
    !print*, "11111111"
    if(s1d_trlx%nbin.eq.0)then
        print*, "error! get_trlx_time_across_cluster: s1d_trlx%nbin=",s1d_trlx%nbin
        stop
    end if
    do i=1, s1d_trlx%nbin
        nm2_tot=0
        do k=1, dms%n
            !!=================
            !print*, "here, spt_rho_rmin should be corrected"
            call get_rho_full_range(dms%mb(k)%all%fmden,dms%mb(k)%all%spt_rho_rmin,&
                s1d_trlx%xb(i),rho_tmp,spp_new)
 
            !!=================                
            nm2_tot=nm2_tot+RHO_TMP*ctl%n0*dms%mb(k)%mc
        end do

        call get_phi_star_full_range(spp_new,s1d_trlx%xb(i),phi_tmp)
        phi_tmp=10**phi_tmp+spp_new%mbh_dmless/10**s1d_trlx%xb(i) 
		call get_v_dispersion_one_ir(dms%mb(1)%all%barge_ir,phi_tmp,rho_tmp,vh) 

        if(nm2_tot>0.and.vh.gt.1d-5)then
            
            call get_two_body_relaxation_time(vh, nm2_tot,lglambda,trlx)
             
            s1d_trlx%fx(i)=trlx!/2d0/pi/1d6
        else
            
            if(nm2_tot<=0.or.vh<=0)then
                if(rid.eq.0)then
                    print*, " get trlx2: x,vh, nm2_tot=", s1d_trlx%xb(i), vh , nm2_tot
                end if
            end if
            s1d_trlx%fx(i)=1d99
        end if
    end do 
end subroutine 
subroutine get_two_body_relaxation_time(vh, nm2_tot, lglambda, trlx) ! in unit of Myr
    use constant
    implicit none
    !type(s1d_type)::fmden
    real(8) trlx,vh, nm2_tot, lglambda
    trlx=0.34*vh**3/nm2_tot/lglambda/2d0/pi/1d6
end subroutine
subroutine get_trlx_min_across_cluster(trlxmin)
    use com_main_gw
    implicit none
    type(s1d_type)::s1d_trlx
    real(8) trlxmin
    call get_trlx_time_across_cluster(s1d_trlx)
    trlxmin=minval(s1d_trlx%fx)
    if(rid.eq.0)then
        call s1d_trlx%print("s1d_trlx")
    end if
end subroutine
subroutine get_simu_time_step(dt)
    use com_main_gw
     
    use md_mbh_evl_acc
    implicit none
    real(8) trlx,dt, tmax_collision,tmin_acc
    real(8) tmin_star_formation
    
    !real(8),parameter::tfractor_collision=0.02

    tmax_collision=1d99
    
    select case(trim(ctl%time_unit))
    case("TNR")
        if(ctl%update_2bdrlx_per_shapshot.ge.1)then
            call get_trlx_min_across_cluster(ctl%trlx_rh0)
        end if
        trlx=ctl%trlx_rh0
       
    case("TRH")
        if(spp_new%mbh.eq.0)then
            print*, "error! TRH0 mode need mbh.ne.0"
            stop
        endif
        if(ctl%update_2bdrlx_per_shapshot.ge.1)then
            call get_trlx_time_at_rh(ctl%trlx_rh0)
        end if
        trlx=ctl%trlx_rh0
    end select
    
    dt=ctl%trlx_rh0*ctl%ts_snap_dt_per_snap

    if(rid.eq.0)then
        print*, "trlx,tfrac=",ctl%trlx_rh0,ctl%ts_snap_dt_per_snap
    end if

    if(spp_new%mbh_dmless>0.and.ctl%enable_evl_mbh.ge.1)then
        if(mbh_mmg%gas_reservior_left>0.1*spp_new%mbh)then
            call get_min_time_due_to_accretion(spp_new%mbh,tmin_acc)
        else
            tmin_acc=1d99
        end if
        dt=min(dt, tmin_acc)
        if(rid.eq.0)then
            print*, "tmin_acc=",tmin_acc
        end if

    end if 

    dt=min(dt,ctl%tmax_timestep)

    if(rid.eq.0)then
        print*, "dt=",dt
    end if
end subroutine
subroutine get_min_time_due_to_accretion(mbh,tmin_acc)
    implicit none
    real(8) tmin_acc,medd,mbh
    real(8),external::mdot_edd_msun_yr
    real(8),parameter::frac=0.1
    medd=mdot_edd_msun_yr(0.1d0,mbh)
    tmin_acc=frac*mbh/medd/1d6 
    ! in Myr, minimum time required to increase more than 10% of current MBH mass.
end subroutine


subroutine main_run()
    use com_main_gw
    implicit none
    integer i
    character*(4) tmpssnapid
    real(8) cur_time_i, cur_time_f,dt
    integer num_update
    real(8) dt0, dt1
    real(8) rho0, rho1
    !print*,  ctl%time_run_mode
    select case(ctl%time_run_mode)
	case(time_run_mode_snap)
        ctl%run_snap_time_i=ctl%run_snap_time_0
        
        do i=ctl%n_spshot_bg+1, ctl%n_spshot_total
            ctl%run_seq_idx=i 
			call get_simu_time_step(dt) 
            ctl%ts_spshot_dt=dt
            ctl%run_snap_time_f=ctl%run_snap_time_i+dt		
			call run_snapshots_sequence(i,ctl%run_snap_time_i,ctl%run_snap_time_f,dt)
            ctl%run_snap_time_i=ctl%run_snap_time_f

		end do
	case(time_run_mode_ttot)
        ctl%run_snap_time_i=ctl%run_snap_time_0
		ctl%run_seq_idx=ctl%n_spshot_bg
			
        do while(ctl%run_snap_time_f<ctl%total_time)
            ctl%run_seq_idx=ctl%run_seq_idx+1

			call get_simu_time_step(dt)
			ctl%ts_spshot_dt=min(dt,ctl%total_time-ctl%run_snap_time_i)
            if(rid.eq.0)then
                print*, "dt,tot-ti=",dt,ctl%total_time,ctl%run_snap_time_i
            end if
            ctl%run_snap_time_f=ctl%run_snap_time_i+ctl%ts_spshot_dt	
            if(ctl%run_snap_time_f.ge.ctl%total_time)then
                ctl%n_spshot_total=ctl%run_seq_idx
            end if
			call run_snapshots_sequence(ctl%run_seq_idx,ctl%run_snap_time_i,ctl%run_snap_time_f,ctl%ts_spshot_dt)
            ctl%run_snap_time_i=ctl%run_snap_time_f
        end do
    case default
        print*, "readin error! define snap mode", ctl%time_run_mode
	end select
    open(unit=199999,file="./output/run_summary.txt",form="formatted")
    write(unit=199999,fmt=*) ctl%run_seq_idx
    close(unit=199999)
end subroutine
subroutine run_snapshots_sequence(i,ti,tf,dt)
    use com_main_gw
     
    use md_dms_saving_data
    use md_event_datas
    implicit none
    integer i    
    character*(6) tmpssnapid, tmprid
    real(8) ti,tf,dt
    integer k
	character*(100) str_, str_update
	type(particle_samples_arr_type)::smstot,smstot_sel
	type(sts_fc_type) ::fc
	type(particle_samples_arr_type),allocatable::smsa(:)
	logical update_dms, output_data_condition
	integer num_update

   	
	allocate(smsa(ctl%ntasks)) 
		 

    write(chattery_out_unit,fmt=*) "start.. snap, cur_time_i, f, dt=", i, ti , tf,dt
    
    if(ctl%trace_all_sample.ge.1)then
        update_dms=.false.
    else
        update_dms=.true.
    end if
    call run_one_snap(ti, tf,  smsa, ctl%ntasks,update_dms)
            
    write(chattery_out_unit,fmt=*) "finished snapshot i, rid, update_j=", i, rid 
    
    if(ctl%trace_all_sample.eq.0)then
        call get_dms_saving_data_all() 
    end if

    if((mod(i,ctl%output_dms_freq).eq.1.or.(ctl%output_dms_freq.eq.1)).and.ctl%trace_all_sample.eq.0)then
        if(rid.eq.mpi_master_id)then
            print*, "start output", rid
            !write(unit=tmprid,fmt="(I4)") rid+1+ctl%ntask_bg
            write(unit=tmpssnapid,fmt="(I6)") i
            call output_diffuse_mspec_bin("output/ecev/dms/dms_"//trim(adjustl(tmpssnapid)))
            call output_dms_hdf5_pdf(dms, "output/ecev/dms/dms_"//trim(adjustl(tmpssnapid)))
            print*, "output finished", rid
        endif
    end if
    if(ctl%output_data_freq.ge.1)then
        output_data_condition=(mod(i,ctl%output_data_freq).eq.1.or.(ctl%output_data_freq.eq.1))&
        .and.ctl%trace_all_sample.eq.0.and.(i.ne.1)
        if(output_data_condition.or.i.eq.ctl%n_spshot_total)then
            print*, "start output bins"
        !===================================
            write(unit=tmprid,fmt="(I6)") rid+1+ctl%ntask_bg
            write(unit=tmpssnapid,fmt="(I6)") i
            str_=trim(adjustl(tmprid))//"_"//trim(adjustl(tmpssnapid))
            call output_chains_bin(bksams,"output/ecev/bin/single/samchn"//trim(adjustl(str_)))
        end if
    end if
    if(ctl%gw_radiation_otby.ge.1.and.ctl%trace_all_sample.eq.0.and.ctl%output_gw_emri_data_bin.ge.1)then
        call collection_all_emris_series()
        if(rid.eq.0)then
            write(unit=tmprid,fmt="(I6)") rid+1+ctl%ntask_bg
            write(unit=tmpssnapid,fmt="(I6)") i
            str_=trim(adjustl(tmpssnapid))
            call output_all_emris_datas("output/ecev/bin/event_data/"//trim(adjustl(str_)))
        end if
        !call output_emris_event_data_series(emris_datas,"output/ecev/bin/"//trim(adjustl(str_)))
    end if
    if(ctl%trace_all_sample.gt.0)then
        if(ctl%trace_all_sample.ge.1)then 
            ctl%output_track_norm=0
            print*, "start output tracks"
            call all_chain_to_arr_single(bksams,bksams_arr)
            call output_sams_sg_track_txt(Bksams_arr, "output/indvd/")
        end if
    end if
    if (.not.(i.eq.ctl%n_spshot_total)) then
        print*, "create_sams"
        call convert_sams()	
    end if
end subroutine
 

subroutine get_num_of_runs(fdir,na)
    implicit none
    character*(*) fdir
    character*(1024) output
    integer na
    integer exitstat
    call execute_command_line('find "'//trim(adjustl(fdir))//'" -maxdepth 1 -type f -name "*.hdf5" | wc -l > output/tmp', &
        exitstat=exitstat, cmdmsg =output, wait=.true.)
        ! print*, "output=",output
    if(exitstat ==0)then
        open(unit=1,file="output/tmp")
        read(unit=1,fmt=*) na
        close(1)
        print*, "number of files:", na
    else
        print*, "error! can not read dir"
        stop
    end if
end subroutine