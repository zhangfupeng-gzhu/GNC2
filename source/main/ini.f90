program ini
	use com_main_gw
	use md_star_pot 
	use md_event_datas
	implicit none
	integer i
	character*(4) tmpi
	logical ex  
	integer ierr, cal_dms
    real(8) t1,t2,t3,t4
	character*(4) arg1
	!type(particle_samples_arr_type)::sps_arr
	
	call get_command_argument(1,arg1)
	if(arg1.ne.'')then
		read(unit=arg1,fmt=*) cal_dms
		cal_dms=0
	else
		cal_dms=1
	end if

	call head("model.in","mfrac.in")

	call set_dm_init(dms) 
	call prepare_tables() 
	if(ctl%include_loss_cone.ge.1)then
		call input_rp_iso(iso_kerr, ctl%rv_nw_conv_dir)
	endif

    if(rid.eq.0)then
	    !call print_model_par()
        call print_current_code_version()
        call cpu_time(t1)
    end if

	if (ctl%intaskmode.eq.task_mode_append) then
		do i=1, ctl%ntask_bg
			write(unit=tmpi,fmt="(I4)") i 
		end do
		do i=ctl%ntask_bg+1, ctl%ntasks+ctl%ntask_bg
			write(unit=tmpi,fmt="(I4)") i
			inquire(file="output/ini/bin/samchn"//trim(adjustl(tmpi))//".bin",exist=ex)
			if(ex) then
				print*, "output/ini/bin/samchn"//trim(adjustl(tmpi))//".bin", " exists"
				stop
			end if
		end do
	end if

	print*, "proc rid start", rid

	call set_seed(ctl%same_rseed_ini, ctl%seed_value+rid) 
	call get_bin_number() 
	if(rid.eq.0)then
		call cpu_time(t3)
	end if
	call get_init_samples(bksams_arr_ini)
	if(rid.eq.0)then
		call cpu_time(t4)
		print*, "MC sample gen consumed time:", t4-t3
	end if
	do i=1, bksams_arr_ini%n
        if(bksams_arr_ini%sp(i)%exit_flag.ne.exit_normal)then
            print*, "i=",i,bksams_arr_ini%sp(i)%exit_flag
            stop
        end if 
    end do 
	call check_weightings("init")
	call set_chain_samples_single(bksams, bksams_arr_ini) 
	
	call get_ini_self_consistent_solution() 
	call init_dms_dc(dms)  
	if(cal_dms.eq.1)then
		call get_dms(dms)
	end if
	
	call MPI_barrier(MPI_COMM_WORLD, ierr)
	call init_output()
	if(rid.eq.0)then
		call cpu_time(t2)		
		print*, "used time=",t2-t1
	end if 
	call stop_mpi()
end
   