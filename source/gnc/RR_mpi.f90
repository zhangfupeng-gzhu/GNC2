subroutine RR_mpi(ch,total_time)
use com_main_gw
use md_mbh_evl_acc
! timeunit=0.158yr
IMPLICIT NONE
	integer n
	integer steps,i,del_flag
    integer,save::try=1
	!integer,parameter:: chunk=4
	logical :: flag_normal
	type(chain_type)::ch
	type(chain_pointer_type),pointer::pt,ps
	real(8) t1,t2,total_time
	interface 
	subroutine run_one_sample_particle(pt, run_time)
		use com_main_gw
		implicit none
		type(chain_pointer_type),target::pt
		real(8) run_time
	end subroutine
	subroutine delete_some_samples(ps,pt,ch,flag)
		use com_main_gw
		implicit none
		type(chain_pointer_type),pointer::pt,ps
		type(chain_type)::ch
		integer flag
	end subroutine
	end interface 
		call cpu_time(t1)
	!end if
    if(ctl%chattery.ge.1)then
		write(chattery_out_unit,fmt=*) 'proc',rid,"starting..."
    end if
	if(rid.eq.mpi_master_id)then
		write(chattery_out_unit,*) "simu begin"
		if(ctl%chattery.ge.1) write(chattery_out_unit,fmt=*) "total number of samples:", bksams%n 
		if(ctl%chattery.ge.1) write(chattery_out_unit,fmt=*) "total time:", total_time,"Myr"
		if(ctl%chattery.ge.1) write(chattery_out_unit,fmt=*) "total number of procs:", ctl%ntasks
	end if
    if(ctl%chattery.ge.1)then
		write(chattery_out_unit,fmt="(A25, A5, 10A7)") "-------result","type","idx", "nhiar","cpuid","eid", "ngene"
	endif
	dc_grid_xstep=dms%dc0%s2_dee%xstep
	dc_grid_ystep=dms%dc0%s2_dee%ystep
	if(ctl%method_interpolate.eq.method_int_linear)then
		call prepare_common_s2ds()
	end if
	running_correction_emax=0
	ctl%num_get_sample_para_exact=0
	ctl%num_get_sample_para_kpl=0
	ctl%num_of_loops=0
	ctl%num_clone_created=0
	ctl%num_clone_elim=0
	ctl%num_of_lc=0
	spp_new%phi_star0=spp_new%phi_r1r2_s+0.5d0*10**(dms%logrmin*2)*spp_new%spt_rho_rmin
	!ctl%bin_mass_Nbalance=0
	call mbh_mmg%init_one_snap()
    !ctl%bin_mass_Nbalance_ot=0
    pt=>ch%head
    !if(rid.eq.7)then
    !	print*, "rid,i, ac,ec=",rid,i, sp(449)%ob%ac, sp(449)%ob%ec
    !end if
    !call pt%bg%output()
    !call bksams%output_screen(2)

	
loop1:	do while(associated(pt))        
		if(pt%ob%exit_flag.eq.exit_normal)then
			select case (ctl%trace_all_sample)
			case(-1)
				if(pt%ob%obtype.eq.star_type_MS)then
					pt%ob%write_down_track=record_track_detail
				else
					pt%ob%write_down_track=0
				end if
			case(-2)
				if(pt%ob%obtype.eq.star_type_BH)then
					pt%ob%write_down_track=record_track_detail
				else
					pt%ob%write_down_track=0
				end if
				
			end select 
			call run_one_sample(pt,total_time) 
		end if 
        ps=>pt
        pt=>pt%next
        !print*, "1"
        if(ctl%clone_scheme.ge.1)then
			if(ps%ob%exit_flag.eq.exit_invtransit &
			.and.ctl%del_cross_clone.ge.1.or.&
			(ctl%del_exit_min.ge.1.and.(ps%ob%exit_flag.eq.exit_boundary_min)))then
				call delete_some_samples(ps,pt,ch,del_flag)
				if(del_flag.eq.-1) exit loop1   ! if head is deleted, chain is empty
				cycle
			end if
        end if
        !call check("after delete in RRMPI")
    end do loop1
	!if(rid.eq.0)then
		call cpu_time(t2)
		write(*,fmt="(A22,F15.4,A3, I4,A40,6I14)")  "run_mpi used time=",t2-t1, " s", rid,&
			" n_exact,kpl, loops,clone,invclone,lc=",&
		 ctl%num_get_sample_para_exact, ctl%num_get_sample_para_kpl, ctl%num_of_loops,&
			 ctl%num_clone_created,ctl%num_clone_elim,ctl%num_of_lc
		call collection_int(running_correction_emax)
		!write(*,fmt="(A22,F15.4,A3, I4,A40,5I20)"), "run_mpi used time=",t2-t1, " s", rid,&
		!	" n_exact,kpl, loops,clone,invclone=",&
		!	 ctl%num_get_sample_para_exact
		!write(*,fmt="(5I20)"),  ctl%num_get_sample_para_exact
	!end if

end subroutine
subroutine prepare_common_s2ds()
	use com_main_gw
	implicit none
	! call dms%dc0%s2_dee%print("s2_dee")
	common_dee_log=dms%dc0%s2_dee
	common_dee_log%fxy=log10(dms%dc0%s2_dee%fxy)
	common_djj_log=dms%dc0%s2_djj
	common_djj_log%fxy=log10(dms%dc0%s2_djj%fxy)
	common_pd_log=dms%pd
	common_pd_log%fxy=log10(dms%pd%fxy)
end subroutine
subroutine delete_some_samples(ps,pt,ch,flag)
	use com_main_gw
	implicit none
	type(chain_pointer_type),pointer::pt,ps
	integer flag
	type(chain_type)::ch
	 
	if(associated(ps%prev))then
		if(allocated(ps%ob))then
			deallocate(ps%ob)
		end if
		!call bksams%output_screen()
		call chain_pointer_delete_item_chain_type(ps)
		!stop
	else
		if(.not. associated(pt)) then
			flag=-1
			return
		end if
			
		!print*, "delete head node cycle"
		!cycle 
		!call ch%output_screen(2,10)
		call destroy_attach_pointer_chain_type(ch%head)
		ch%head=>pt
		call pt%set_head()     
		!print*, "delete the fisrt node",rid
		pt%prev=>null()
		 
		
	end if

end subroutine