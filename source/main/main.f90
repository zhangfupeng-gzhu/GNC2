program main
    use com_main_gw
!	use md_coeff
	implicit none
	integer i,j,ierror
	real(8) t1,t2
	!character*(4) tmprid, tmpssnapid
	character*(9) str_
	logical ex
    INTEGER*4 getpid, pid
	integer num_update
	character*(20) pid_root_str, pid_cld_str

    call head("model.in","mfrac.in")
     
	select case(ctl%insnapmode)
	case(snap_mode_new,snap_mode_one)
		call prepare_ini_data()
		
	case(snap_mode_append)
		call prepare_ini_data_append()
	end select 

	if(rid.eq.0)then
		call cpu_time(t1)
	end if	
	!call set_seed(ctl%same_rseed_evl, ctl%seed_value+rid)
    !if(rid.eq.0)then
	!	print*, "total simu time, tnr=", ctl%total_time,ctl%tnr
	!end if
	!stop
	!print*, ctl%time_run_mode
	call main_run()
	if(rid.eq.0)then
		call cpu_time(t2)
		print*, "total time:", t2-t1
	end if
	print*, "finished main", rid
!	call mpi_barrier(mpi_comm_world,ierror)
	call end()

end
