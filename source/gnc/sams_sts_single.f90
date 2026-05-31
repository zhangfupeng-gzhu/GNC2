
 
subroutine get_sample_select_numbers(sa,flag,n,n_w)
	use com_main_gw
	implicit none
	integer flag
	type(particle_samples_arr_type)::sams_sel,sa
	real(8),allocatable::weight_real(:)
	integer n
	real(8) n_w
	call sa%select(sams_sel,flag,-1,-1d0,-1d0)

	if(sams_sel%n>0)then
		allocate(weight_real(sams_sel%n))
		weight_real(1:sams_sel%n)=sams_sel%sp(1:sams_sel%n)%weight_real
		n=sams_sel%n
		n_w=sum(weight_real)
	else
		n=0
		n_w=0
	end if

end subroutine

subroutine collect_to_root_emri_single(em_send,em, n)
	use mpi_comu
	use md_event_datas
	implicit none
	integer i,n
	type(event_data_series)::em(n),em_send
	
	do i=1, ctl%ntasks
		if(i.ne.mpi_master_id+1)then
			!print*, "rid,i=",rid,i
			call send_emri_data_arr_mpi(em_send,em(i), i-1, mpi_master_id)
		else
			em(mpi_master_id+1)=em_send
		end if
	end do
end subroutine

subroutine send_emri_data_arr_mpi(em_send, em_recv, proc_id_source,proc_id_dest)
	use mpi_comu
	use md_event_datas
	implicit none
	type(event_data_series)::em_send,em_recv
	integer ierr, proc_id_source,proc_id_dest, i, nintarr, nrealarr
	integer status(MPI_Status_size)
	integer,allocatable::intarr(:,:)
	real(8),allocatable::realarr(:,:)

	if(rid.eq.proc_id_source)then
		allocate(intarr(nint_em,em_send%n))
		allocate(realarr(nreal_em,em_send%n))
		call conv_em_int_real_arrays(em_send, intarr(1:nint_em,1:em_send%n), realarr(1:nreal_em,1:em_send%n))

		call mpi_send(em_send%n,1, MPI_INTEGER,proc_id_dest,0,MPI_COMM_WORLD,ierr)
		call mpi_send(intarr, nint_em*em_send%n, MPI_INTEGER, proc_id_dest, 0, MPI_COMM_WORLD, ierr)
		call mpi_send(realarr, nreal_em*em_send%n, MPI_DOUBLE_PRECISION, proc_id_dest, 0, MPI_COMM_WORLD, ierr)
	elseif(rid.eq.proc_id_dest)then
		call mpi_recv(em_recv%n,1, MPI_INTEGER,proc_id_source,0,MPI_COMM_WORLD,status,ierr)
		call em_recv%init(em_recv%n)
		allocate(intarr(nint_em,em_recv%n))
		allocate(realarr(nreal_em,em_recv%n))
		nintarr=nint_em*em_recv%n; nrealarr=nreal_em*em_recv%n
		call mpi_recv(intarr, nintarr, MPI_INTEGER, proc_id_source, 0, mpi_comm_world,status, ierr)
		call mpi_recv(realarr, nrealarr, MPI_DOUBLE_PRECISION, proc_id_source, 0, mpi_comm_world, status,ierr)
		call conv_int_real_arrays_em(em_recv, intarr(1:nint_em,1:em_recv%n), realarr(1:nreal_em,1:em_recv%n))

	end if
end subroutine

subroutine collection_all_emris_series()
	use com_main_gw
	use md_event_datas
	implicit none
	logical,external::test_if_component_exists
	!print*, "test if sbh exist=", test_if_component_exists(star_type_bh)
	if(test_if_component_exists(star_type_ms))then
		call collection_a_event_data(data_star)
	end if

	if(test_if_component_exists(star_type_bh))then
		call collection_a_event_data(data_sbh)
	end if
	if(test_if_component_exists(star_type_ns))then
		call collection_a_event_data(data_ns)
	end if
	if(test_if_component_exists(star_type_wd))then
		call collection_a_event_data(data_wd)
	end if
	if(test_if_component_exists(star_type_bd))then
		call collection_a_event_data(data_bd)
	end if
	if(test_if_component_exists(star_type_rg))then
		call collection_a_event_data(data_rg)
	end if
end subroutine
subroutine collection_a_event_data(ed)
	use md_event_datas
	implicit none
	type(event_data)::ed
	call collection_a_event_data_series(ed%emris)
	call collection_a_event_data_series(ed%td)
end subroutine
subroutine collection_a_event_data_series(em)
	use com_main_gw
	use md_event_datas
	implicit none
	type(event_data_series)::em, em_collect(ctl%ntasks)!, em_tot
	integer i, n_tot,j,k


	call collect_to_root_emri_single(em,em_collect,ctl%ntasks)
	if(rid.eq.0)then
		n_tot=sum(em_collect(:)%n)
				
		call em%init(n_tot)
		j=0
		do i=1, ctl%ntasks
			do k=1, em_collect(i)%n
				j=j+1
				call copy_i_event_data_series(em_collect(i),k,em,j)
			end do
		end do
		!call output_emris_event_data_series(em_tot, fl)
	end if
end subroutine

subroutine dms_saving_data_get_sts(sars,se)
	use com_main_gw
	use md_dms_saving_data
	implicit none
	type(particle_samples_arr_type)::sars
	type(snap_event)::se
	real(8),allocatable::x(:),w(:)
	real(8) nr,nw,m1,m2
	character*(10) star_type_string
	integer n
	integer i,ier

	nr=sars%n
	n=sars%n
	if(allocated(x))then
		deallocate(x,w)
	end if
	allocate(x(n),w(n))
	
	w(1:n)=sars%sp(1:n)%weight_real
	nw=sum(w)

	call collection_event_numbers(nr,nw)
	
	se%n=nr
	se%nw=nw
	se%rate=nw/(ctl%run_snap_time_f-ctl%run_snap_time_i)
	call collection_snap_events(se,ctl%run_snap_time_f-ctl%run_snap_time_i)

	if(nw>0)then
		associate(fdstr=>se%fdstr_x)
			x(1:n)=sars%sp(1:n)%x
			call fdstr%init(sample_logemin,sample_logemax,ctl%dstr_bins_e,use_weight=.true.)
			call fdstr%set_range()
			if(n.ne.0)then
				call fdstr%get_hst(log10(x),w,n)
			else
				fdstr%fxw=0
			end if
			call collection_and_avg_fx(fdstr%nbw,fdstr%nbin)
			call collection_and_avg_fx(fdstr%fxw,fdstr%nbin)
		end associate

		associate(fdstr=>se%fdstr_m)
			x(1:n)=sars%sp(1:n)%m
			m1=minval(x(1:n))
			m2=maxval(x(1:n))		
			
			call collection_and_get_min_real(m1)
			call collection_and_get_max_real(m2)

			if(rid.eq.0)then
				if(n>0)then
					call get_star_type(sars%sp(1)%obtype,star_type_string)
					print*, "type,min max mass=",star_type_string, m1,m2
				end if
			end if
			m1=m1*0.9999
			m2=m2*1.0001
			call fdstr%init(m1,m2,&
				ctl%dstr_bins_e,use_weight=.true.)
			call fdstr%set_range()
			if(n>0)then
				call fdstr%get_hst(x,w,n)
			else
				fdstr%fxw=0
			end if
 
			call collection_fx_int(fdstr%nb,fdstr%nbin)
			call collection_and_avg_fx(fdstr%nbw,fdstr%nbin)
			call collection_and_avg_fx(fdstr%fxw,fdstr%nbin)
		end associate
	end if 
end subroutine

subroutine dms_saving_data_get_sts_td(sars,se)
	use com_main_gw
	use md_dms_saving_data
	implicit none
	type(particle_samples_arr_type)::sars
	type(snap_event)::se
	real(8),allocatable::x(:),w(:)
	real(8) nr,nw,m1,m2,xmin,xmax
	character*(10) star_type_string
	integer n
	integer i,ier

	nr=sars%n
	n=sars%n
	if(allocated(x))then
		deallocate(x,w)
	end if
	allocate(x(n),w(n))
	
	w(1:n)=sars%sp(1:n)%weight_real
	nw=sum(w)

	call collection_event_numbers(nr,nw)
	
	se%n=nr
	se%nw=nw
	se%rate=nw/(ctl%run_snap_time_f-ctl%run_snap_time_i)
	call collection_snap_events(se,ctl%run_snap_time_f-ctl%run_snap_time_i)

	if(nw>0)then
		associate(fdstr=>se%fdstr_x)
			x(1:n)=sars%sp(1:n)%rp/sars%sp(1:n)%r_lc
			xmin=0d0
			xmax=1d0
			call fdstr%init(xmin,xmax,20,use_weight=.true.)
			call fdstr%set_range()
			if(n.ne.0)then
				call fdstr%get_hst(x,w,n)
			else
				fdstr%fxw=0
			end if
			call collection_and_avg_fx(fdstr%nbw,fdstr%nbin)
			call collection_and_avg_fx(fdstr%fxw,fdstr%nbin)
		end associate

		associate(fdstr=>se%fdstr_m)
			x(1:n)=sars%sp(1:n)%m
			m1=minval(x(1:n))
			m2=maxval(x(1:n))		
			
			call collection_and_get_min_real(m1)
			call collection_and_get_max_real(m2)

			if(rid.eq.0)then
				if(n>0)then
					call get_star_type(sars%sp(1)%obtype,star_type_string)
					print*, "type,min max mass=",star_type_string, m1,m2
				end if
			end if
			m1=m1*0.9999
			m2=m2*1.0001
			call fdstr%init(m1,m2,&
				ctl%dstr_bins_e,use_weight=.true.)
			call fdstr%set_range()
			if(n>0)then
				call fdstr%get_hst(x,w,n)
			else
				fdstr%fxw=0
			end if 
			call collection_fx_int(fdstr%nb,fdstr%nbin)
			call collection_and_avg_fx(fdstr%nbw,fdstr%nbin)
			call collection_and_avg_fx(fdstr%fxw,fdstr%nbin)
		end associate
	end if 
end subroutine

subroutine get_dms_saving_data(star_type_number,oe,ed)
	use com_main_gw
	use md_dms_saving_data
	use md_event_datas
	implicit none
	integer n,ierr
	integer star_type_number
	type(particle_samples_arr_type)::bksma_arr_sel,bksma_arr_emri
	type(particle_samples_arr_type)::bksma_arr_ms_td,bksma_arr_lc
	real(8) nr, nw
	real(8),allocatable::x(:),w(:)
	logical,external::selection_emri,selection_td
	type(obj_events)::oe
	type(event_data)::ed

	call sams_arr_select_type_single(bksams_arr, bksma_arr_sel,star_type_number) 
	!print*, "selection"
	call get_sample_select_numbers(bksma_arr_sel,exit_boundary_max,n, nw)
	nr=n
	!print*, "nr, nw=", nr, nw
	
	call collection_event_numbers(nr,nw)
	
	oe%se_emax%n=nr
	oe%se_emax%nw=nw
	oe%se_emax%rate=nw/(ctl%run_snap_time_f-ctl%run_snap_time_i)

	call collection_snap_events(oe%se_emax,ctl%run_snap_time_f-ctl%run_snap_time_i)

	if(ctl%include_loss_cone.ge.1)then
		call sams_arr_select_condition_single(bksma_arr_sel, bksma_arr_ms_td,selection_td)
		call dms_saving_data_get_sts_td(bksma_arr_ms_td,oe%se_td)
		call dms_saving_data_record(bksma_arr_ms_td,ed%td)
		call bksma_arr_sel%select(bksma_arr_lc,exit_lc,-1,-1d0,-1d0)
		call dms_saving_data_get_sts(bksma_arr_lc,oe%se_lc)
	end if

	if(ctl%gw_radiation_otby.ge.1)then
		!call get_sample_select_numbers(bksma_arr_sel,exit_emri_single,n, nw)
		!nr=n
		
		call sams_arr_select_condition_single(bksma_arr_sel,bksma_arr_emri,selection_emri)
		print*, "emris:spsel%n, spobj%n, type,rid=",bksma_arr_sel%n, bksma_arr_emri%n, &
			star_type_str(star_type_number),rid
		!call bksma_arr_sel%select(bksma_arr_emri,exit_emri_single,-1,-1d0,-1d0)
		call dms_saving_data_get_sts(bksma_arr_emri,oe%se_emris)
		call dms_saving_data_record(bksma_arr_emri,ed%emris)
		call dms_saving_data_get_emris_sts(bksma_arr_emri,oe)

	end if

end subroutine

logical function selection_emri(sp)
	use md_particle_sample
	implicit none
	type(particle_sample_type)::sp
	integer flag
	!real(8) tgw, GET_T_GW
	!real(8) ac,ec

	flag=sp%exit_flag
	selection_emri=.false.
	if(flag.eq.exit_emri_single)then
		selection_emri=.true.
	end if
	if(flag.eq.exit_boundary_max)then		
		if(sp%state_emri_last.ge.1.and. sp%state_emri_current.ge.1)then	
			selection_emri=.true.
		end if
		!tgw=get_t_gw(sp%m,spp_new%mbh,ac, ec)
	end if
end function

subroutine dms_saving_data_record(sp,em)
	use com_main_gw
	!use md_dms_saving_data
    use md_event_datas
	implicit none
	type(particle_samples_arr_type)::sp
	type(event_data_series)::em
	integer i
	call em%init(sp%n)

	do i=1, em%n
		em%x(i)=sp%sp(i)%x
		em%jm(i)=sp%sp(i)%jm
		em%ac(i)=sp%sp(i)%byot%a_bin
		em%ec(i)=sp%sp(i)%byot%e_bin
		em%rp(i)=sp%sp(i)%rp
		em%ra(i)=sp%sp(i)%ra
		em%r_ls(i)=sp%sp(i)%r_lc
		em%mass(i)=sp%sp(i)%m
		em%w(i)=sp%sp(i)%weight_real
		em%inc(i)=sp%sp(i)%byot%Inc
		em%radius(i)=sp%sp(i)%byot%ms%radius
		em%exit_flag(i)=sp%sp(i)%exit_flag
		em%nlvl(i)=sp%sp(i)%Lvl_clone
	end do
end subroutine


subroutine dms_saving_data_get_emris_sts(sp,oe)
	use com_main_gw
	use md_dms_saving_data
	implicit none
	type(particle_samples_arr_type)::sp
	type(obj_events)::oe
	real(8) weights(sp%n), x(sp%n), m1(sp%n), m2(sp%n),mtot(sp%n),enx(sp%n)
	real(8) jm(sp%n),rp(sp%n), jc(sp%n), ac(sp%n),ec(sp%n)
	real(8) agw(sp%n),egw(sp%n),rangemax,rangemin!,get_clone_deep
	integer n, ns, i,nl(sp%n),flag(sp%n)
	real(8),parameter::freqob=0.1 !mHz
	real(8),parameter::freq_simu=freqob/1d3*86400*365.2425/2d0/pi
	real(8) acrit(sp%n)
	n=sp%n

	m1=sp%sp(1:n)%m; m2=spp_new%mbh
	jm=sp%sp(1:n)%jm; 
	enx=log10(sp%sp(1:n)%x)
	rp=sp%sp(1:n)%rp/r0_cl; 
	jc=sp%sp(1:n)%jc/(r0_cl*ctl%v0) 

	do i=1, n
		acrit(i)=(1/freq_simu/(2d0/pi)*(m1(i)+m2(i)))**(1d0/3d0)
	end do
	do i=1, n
		call get_kpl_ae(jm(i),jc(i),spp_new%mbh_dmless, rp(i),ac(i),ec(i))
	end do
	ac=ac*r0_cl
	
	call get_ae_given_fgw(ac, ec, m1 ,m2, n, freq_simu, agw,egw)
	!print*, "1111"
	ns=0
	do i=1, n
		if(egw(i)<1)then
			if(sp%sp(i)%exit_flag.eq.exit_lc.and.ac(i)>acrit(i)) cycle
			ns=ns+1
			x(ns)=1-egw(i)
			weights(ns)=sp%sp(ns)%weight_real	
		end if
	end do

	associate(fdstr=>oe%fd_emris_ecc)
		rangemin=-2
		rangemax=0d0!.5
		!if(maxval(x)>1.5) print*, "egw=",egw(maxloc(x,dim=1))
		call fdstr%init(rangemin,rangemax,50,fc_spacing_linear,use_weight=.true.)
		call fdstr%set_range()
		!call set_range_fc
		if(ns.ne.0)then
			call get_fc_weight(log10(x(1:ns)),weights(1:ns),ns,fdstr)
			!call fdstr%get_hst(log10(x),weights,n)
		else
			fdstr%fxw=0
			fdstr%nbw=0
		end if
		call collection_and_avg_fx(fdstr%nbw,fdstr%nbin)
		call collection_and_avg_fx(fdstr%fxw,fdstr%nbin) 
	end associate
	
	associate(f2dstr=>oe%fd_emris_nxj_ir) 
		f2dstr=dms%mb(1)%dsp(1)%p%nxj_ir
		call f2dstr%get_hst(enx,log10(jm),weights,n)
		call collection_and_avg_s2d(f2dstr%nxyw, f2dstr%nx,f2dstr%ny)
	end associate
end subroutine  

subroutine get_dms_saving_data_all()
	use com_main_gw
	use md_dms_saving_data
	use md_event_datas
	implicit none
	integer n,ierr
	type(particle_samples_arr_type)::bksma_arr_sel
	real(8) nr, nw
	logical,external::test_if_component_exists

	!print*, "get_dms_saving_data"
	call all_chain_to_arr_single(bksams,bksams_arr)
	!call sams_get_weight_clone_single(bksams_arr)
	!call set_real_weight_arr_single(bksams_arr)

	call get_dms_saving_data(star_type_ms,oe_star,data_star)
	if(test_if_component_exists(star_type_rg))then
		call get_dms_saving_data(star_type_rg,oe_rg,data_rg)
	end if
	if(test_if_component_exists(star_type_bh))then
		call get_dms_saving_data(star_type_bh,oe_sbh,data_sbh)
	end if
	if(test_if_component_exists(star_type_wd))then
		call get_dms_saving_data(star_type_wd,oe_wd,data_wd)
	end if
	if(test_if_component_exists(star_type_bd))then
		call get_dms_saving_data(star_type_bd,oe_bd,data_bd)
	end if
	if(test_if_component_exists(star_type_ns))then
		call get_dms_saving_data(star_type_ns,oe_ns,data_ns)
	end if

end subroutine 
 
subroutine get_ae_given_fgw(a,e,m1,m2,n,fgw,agw,egw)
	implicit none
	integer n,i
	real(8) fgw
	real(8) a(n), e(n), c0(n),  agw(n), egw(n)
	real(8) m1(n), m2(n)
	do i=1, n
		call get_c0(a(i),e(i),m1(i),m2(i),c0(i))
		call get_aeob_given_fgwc0(m1(i),m2(i),c0(i), fgw, agw(i),egw(i))
	end do
end subroutine
subroutine get_ae_given_fgw_one(a,e,m1,m2,fgw,agw,egw)
	implicit none
	real(8) fgw
	real(8) a, e, c0,  agw, egw
	real(8) m1, m2
	
	call get_c0(a,e,m1,m2,c0)
	call get_aeob_given_fgwc0(m1,m2,c0, fgw, agw,egw)

end subroutine
logical function selection_td(sp)
	use com_main_gw
	implicit none
	type(particle_sample_type)::sp
	integer flag
	flag=sp%exit_flag

	selection_td=.false.
	if(flag.eq.exit_tidal_full.or.flag.eq.exit_tidal_empty.and.sp%m>0.01d0)then
		selection_td=.true.
	end if
end function



subroutine output_all_emris_datas(fl)
	use md_event_datas
	implicit none
	character*(*) fl
	integer,parameter::funit=99123123
	open(unit=funit,file=trim(adjustl(fl))//"_emdata.bin",form='unformatted', &
		access='stream')
		call data_sbh%write(funit)
		call data_bd%write(funit)
		call data_ns%write(funit)
		call data_wd%write(funit)
		call data_star%write(funit)
		call data_rg%write(funit)
	close(funit)
end subroutine

subroutine input_all_emris_datas(fl)
	use md_event_datas
	implicit none
	character*(*) fl
	integer,parameter::funit=99123123
	open(unit=funit,file=trim(adjustl(fl))//"_emdata.bin",form='unformatted', &
		access='stream',status="old")
		
		call data_sbh%read(funit)
		! print*, "sbh"
		call data_bd%read(funit)
		! print*, "bd"
		call data_ns%read(funit)
		! print*, "ns"
		call data_wd%read(funit)
		! print*, "wd"
		call data_star%read(funit)
		call data_rg%read(funit)
		! print*, "star"
	close(funit)
end subroutine

subroutine save_td_hdf5(em, group_id,star_type)
	use constant
	use md_event_datas
	use md_hdf5
	use md_star_pot
	use md_particle_sample
	implicit none
	type(event_data_series)::em,em_copy
	integer(HID_T) group_id
	! character*(10) tablename
	type(hdf5_table_type)::h5table
	integer::nfields
	real(8) mbhmass, m2(em%n)
	real(8)::pd, acmin, r_lc,mmin_in,mmax_in
	integer n_sel,i,j,star_type
	logical::uselog
	if(em%n>0)then
		mbhmass=spp_new%mbh
		m2=mbhmass
		
		n_sel=em%n
		print*, "td:em%n=",em%n
		call em_copy%init(n_sel)
		j=0
		do i=1, em%n
			j=j+1
			call copy_i_event_data_series(em,i,em_copy,j)
		end do

		select case(star_type)
		case(star_type_bh,star_type_ns,star_type_wd)
			return
		end select
		if(n_sel>0)then
			! print*, "n_sel=", n_sel
			nfields=11
			call h5table%init_table(nfields, em_copy%n, "Data_td")
			h5table%field_names=(/"  x", " jm", "  w", " ac", " ec", " rp", " ls", "  m",  "rad", "inc", "flg"/)
			h5table%field_types(1:10)=H5T_NATIVE_DOUBLE
			h5table%field_types(11)=H5T_NATIVE_INTEGER
			call h5table%prepare_write_table(group_id)
			call h5table%write_column_real(em_copy%x(:))
			call h5table%write_column_real(em_copy%jm(:))
			call h5table%write_column_real(em_copy%w(:))
			call h5table%write_column_real(em_copy%ac(:))
			call h5table%write_column_real(em_copy%ec(:))
			call h5table%write_column_real(em_copy%rp(:))
			call h5table%write_column_real(em_copy%r_ls(:))
			call h5table%write_column_real(em_copy%mass(:))
			call h5table%write_column_real(em_copy%radius(:))
			call h5table%write_column_real(em_copy%inc(:))
			call h5table%write_column_int(em_copy%exit_flag(:))
		end if
		
		call save_td_sts_hdf5(em,group_id )
	end if
end subroutine
subroutine save_td_sts_hdf5(em, group_id )
	use md_event_datas
	use md_hdf5
	use com_sts_type
	use constant
	use model_basic
	implicit none
	type(event_data_series)::em
	integer(HID_T) group_id
	! character*(*) tablename
	type(hdf5_table_type)::h5table
	integer::nfields
	real(8) mass(em%n), radius(em%n),ecc(em%n),weight(em%n)
	type(s1d_hst_type)::fr, fm
	type(s2d_hst_type)::fmr2d
	real(8) xmin,xmax,ymin,ymax
	logical::uselog=.true.
	integer,parameter::nr=20,nm=20

	if(em%n>0)then
		xmin=0.01d0; xmax=150d0
		if(uselog)then
			xmin=log10(xmin); xmax=log10(xmax)
			mass=log10(em%mass(:))
		else
			xmin=xmin; xmax=xmax
			mass=em%mass(1:em%n)
		end if
		weight=em%w(1:em%n)/dble(ctl%ntask_total)
		call fm%init(xmin,xmax,nm,use_weight=.True.)
		call fm%set_range()
		call fm%get_hst(mass(1:em%n),weight(1:em%n),em%n)
		call fm%save_hdf5(group_id,"fmtd")

		ymin=0.05; ymax=1d4
		ymin=log10(ymin); ymax=log10(ymax)
		radius=log10(em%radius(:)/rd_sun)
		call fr%init(ymin,ymax,nr,use_weight=.True.)
		call fr%set_range()
		call fr%get_hst(radius(1:em%n),weight(1:em%n),em%n)
		call fr%save_hdf5(group_id,"frtd")

		call fmr2d%init(nm,nr,xmin,xmax,ymin,ymax,use_weight=.true.)
		call fmr2d%set_range()
		call fmr2d%get_stats_weight(mass(1:em%n),radius(1:em%n),weight(1:em%n),em%n)
		call fmr2d%save_hdf5(group_id,"fmrtd")
	end if
end subroutine
subroutine save_emris_hdf5(em, group_id,star_type)
	use constant
	use md_event_datas
	use md_hdf5
	use md_star_pot
	use md_particle_sample
	use com_main_gw,only:jmax_value
	implicit none
	type(event_data_series)::em,em_copy
	integer(HID_T) group_id,attr_id
	! character*(10) tablename
	type(hdf5_table_type)::h5table
	integer::nfields
	real(8) agw(em%n),egw(em%n), agwtmp, egwtmp
	real(8) freqmhz,emin_possible, rnd
	real(8) mbhmass, m2(em%n),dwe(em%n),tgw(em%n),dwj(em%n),gwf(em%n), see(em%n),sjj(em%n)
	real(8)::freq_simu,pd, acmin, r_lc,mmin_in,mmax_in
	integer n_sel,i,j,star_type
	logical::uselog
	! character*(3) star_type_str

	if(em%n>0)then
		freqmhz=0.1d0
		freq_simu=freqmhz/1d3*86400*365.2425/2d0/pi
		mbhmass=spp_new%mbh
		m2=mbhmass
		! print*,"jmax=", maxval(em%jm)
		j=0
		do i=1, em%n
			call get_ae_given_fgw_one(em%ac(i), em%ec(i), em%mass(i), &
			m2(i), freq_simu, agwtmp,egwtmp)
			if(agwtmp*(1-egwtmp)>em%r_ls(i))then
				j=j+1
			end if
			if(em%w(i)<1d-10)then
				print*, "warnning, i, em%w(i)=",i, em%w(i)
			end if
		end do
		n_sel=j
		print*, "em%n,n_sel=",em%n,n_sel
		call em_copy%init(n_sel)

		j=0
		do i=1, em%n
			call get_ae_given_fgw_one(em%ac(i), em%ec(i), em%mass(i), &
			m2(i), freq_simu, agwtmp,egwtmp)
			if(agwtmp*(1-egwtmp)>em%r_ls(i))then
				j=j+1
				call copy_i_event_data_series(em,i,em_copy,j)
				! emin_possible=(1-0.99**2)**0.5   ! to correct the results is not very accurate when e is very close to zero and a is very small.
				! if(em%jm(i)>0.99)then
				! 	egwtmp=rnd(0d0,emin_possible)
				! end if
				agw(j)=agwtmp
				egw(j)=egwtmp
			end if
		end do
		
		! call em_copy%init(n_sel)
		call get_dwf( agw(1:n_sel), egw(1:n_sel),em_copy%rp(1:n_sel), em_copy%mass(1:n_sel), n_sel, dwe,dwj,1d0,gwf,see,sjj,tgw)
		! em_copy=em
		! n_sel=em%n
		
		
		! call em%get_at_same_freq(spp_new%mbh,0.1d0, acout,ecout)
		if(n_sel>0)then
			! call get_star_type(star_type,star_type_str)
			call add_attr_dble(group_id,attr_id, "frac(0.1mHz)", dble(n_sel)/dble(em%n))
			! print*, "n_sel=", n_sel
			nfields=18
			call h5table%init_table(nfields, em_copy%n, "Data_emris")
			h5table%field_names=(/"  x", " jm", "  w", " ac", " ec", " rp", " ls", "  m", "rad", "agw",&
			                      "egw" , "inc", "dwe", "dwj", "gwf","tgw","lvl", "flg"/)
			h5table%field_types(1:16)=H5T_NATIVE_DOUBLE
			h5table%field_types(17:18)=H5T_NATIVE_INTEGER
			call h5table%prepare_write_table(group_id)
			call h5table%write_column_real(em_copy%x(:))
			call h5table%write_column_real(em_copy%jm(:))
			call h5table%write_column_real(em_copy%w(:))
			call h5table%write_column_real(em_copy%ac(:))
			call h5table%write_column_real(em_copy%ec(:))
			call h5table%write_column_real(em_copy%rp(:))
			call h5table%write_column_real(em_copy%r_ls(:))
			call h5table%write_column_real(em_copy%mass(:))
			call h5table%write_column_real(em_copy%radius(:))
			call h5table%write_column_real(agw(1:n_sel))
			call h5table%write_column_real(egw(1:n_sel))
			call h5table%write_column_real(em_copy%inc(:))
			call h5table%write_column_real(dwe(:))
			call h5table%write_column_real(dwj(:))
			call h5table%write_column_real(gwf(:))
			call h5table%write_column_real(tgw(:))
			call h5table%write_column_int(em_copy%nlvl(:))
			call h5table%write_column_int(em_copy%exit_flag(:))
			
		end if
		select case(star_type)
		case(star_type_bh)
			mmin_in=5d0;mmax_in=25d0
			uselog=.false.
		case(star_type_wd,star_type_ns)
			mmin_in=0.3d0;mmax_in=2.0d0
			uselog=.false.
		! case(star_type_ms)
			! mmin_in=0.1d0;mmax_in=0.1d0
			! uselog=.true.
		case(star_type_bd,star_type_ms)
			mmin_in=0.01d0;mmax_in=1.2d0
			uselog=.true.
		end select
		
		call save_emris_sts_hdf5(em_copy,agw(1:n_sel),egw(1:n_sel),dwe(1:n_sel),dwj(1:n_sel),&
			group_id,mmin_in,mmax_in,uselog)
	end if
end subroutine
subroutine get_dwf(aw,ew,rpi,mass,n,dwe,dwj,obtime,gwf,see,sjj,tgw)
	use com_main_gw
	implicit none
	integer n,i,j
	real(8) aw(n),ew(n),mass(n),dwe(n),dwj(n),tgw(n),rpi(n),gwf(n),see(n),sjj(n)
	real(8) x,jm,rdx,rdy,even,evjm
	integer idx,idy
	real(8) de_110(dms%n),dj_111(dms%n),de_0,dee,dj_rest,djj,dej
	type(coeff_type)::coenr,coegw
	real(8) rp, pd,en,jc
	real(8) obtime ! in unit of yr

	call prepare_common_s2ds()
	do i=1, n
		x=spp_new%mbh/aw(i)/2d0/ctl%v0**2
		if(ew(i)>1)then
			print*, "?? ew(i)>1", ew(i)
			stop
		end if
		jm=(1-ew(i)**2)**(0.5d0)
		jc=(aw(i)*spp_new%mbh)**0.5/(ctl%v0*r0_cl)
		en=spp_new%mbh/2d0/aw(i)
		
		dc_grid_xstep=dms%dc0%s2_dee%xstep
		dc_grid_ystep=dms%dc0%s2_dee%ystep

		call get_ex_idx_ir(x, idx,rdx,even)
		call get_jm_idx(jm, idy,rdy,evjm)
		! print*, "emax,emin=",dms%dc0%s2_dee%xcenter(96),dms%dc0%s2_dee%xcenter(1)
		! print*, "x,jm,idx, rdx, idy, rdy=", x,jm, idx, rdx, idy, rdy


		do j=1, dms%n
            !print*, "dms%df_coe_bins=",dms%df_coe_bins
            call linear_int_2d_xy(idx,idy,rdx,rdy,&
                dms%mb(j)%dc%s2_de_110%fxy,dms%df_coe_bins,dms%df_coe_bins,de_110(j))
            call linear_int_2d_xy(idx,idy,rdx,rdy,&
                dms%mb(j)%dc%s2_dj_111%fxy,dms%df_coe_bins,dms%df_coe_bins,dj_111(j))
        end do

		associate(dc0=>dms%dc0)
            call linear_int_2d_xy(idx,idy,rdx,rdy,&
                dc0%s2_de_0%fxy,dms%df_coe_bins,dms%df_coe_bins,de_0)

            call linear_int_2d_xy(idx,idy,rdx,rdy,&
                common_dee_log%fxy,dms%df_coe_bins,dms%df_coe_bins,dee)
            !print*, "idx,idy,rdx,rdy,dee=",idx,idy,rdx,rdy,dee

            dee=10**dee
            call linear_int_2d_xy(idx,idy,rdx,rdy,&
                dc0%s2_dj_rest%fxy,dms%df_coe_bins,dms%df_coe_bins,dj_rest)
        
            call linear_int_2d_xy(idx,idy,rdx,rdy,&
                common_djj_log%fxy,dms%df_coe_bins,dms%df_coe_bins,djj)
            djj=10**djj
        
            call linear_int_2d_xy(idx,idy,rdx,rdy,&
                dc0%s2_dej%fxy,dms%df_coe_bins,dms%df_coe_bins,dej)    
        end associate

		coeNr%jj=djj*jc*jc; 
        coeNr%e=de_0; coeNr%j=dj_rest
        do j=1, dms%n
            coeNr%e=coeNr%e+mass(i)/dms%mb(j)%mc*de_110(j)
            coeNr%j=coeNr%j+dj_111(j)*(mass(i)+dms%mb(j)%mc)/dms%mb(j)%mc/2d0; 
        end do
        coeNr%ee=dee*en*en;
        !coeNr%ee=dee*(10**even*ctl%energy0)**2

        coeNr%e=coeNr%e*en; 
        !coeNr%e=coeNr%e*(10**even*ctl%energy0)

        coeNR%j=  coeNR%j*jc
        coeNr%ej=  dej*en*jc
        if(ctl%gw_radiation_otby.ge.1)then
            sample_rlx_e_time=1d0/dee
            sample_rlx_j_time=10**(evjm*2)/djj
		
			! read(*,*)
        end if

		rp=aw(i)*(1-ew(i))
		pd=(aw(i)**3/spp_new%mbh)**0.5*2*pi
		call get_coegw(rp, ew(i), mass(i), pd,coegw)

		sample_tgw_time=abs(en/coegw%e)
		
		dwe(i)=(coenr%ee*obtime*2*pi)**0.5/abs(coegw%e*obtime*2*pi) !min(sample_rlx_e_time,sample_rlx_j_time)/sample_tgw_time
		dwj(i)=(coenr%jj*obtime*2*pi)**0.5/abs(coegw%j*obtime*2*pi) !min(sample_rlx_e_time,sample_rlx_j_time)/sample_tgw_time
		! print*, "rp,aw,ew,mbh=",rp,aw(i),ew(i),rpi(i),spp_new%mbh
		!print*, "dwe,dwj, coenr%ee, dee, en=", dwe(i),dwj(i), coenr%ee, coenr%jj, dee, djj, coegw%j

		! print*, "pd, coegw%e, sample_tgw_time=",pd/2/pi, coegw%e, mass(i), sample_tgw_time
		! print*, "sample_tgw_time/min(sample_rlx_e_time,sample_rlx_j_time)=",sample_tgw_time/min(sample_rlx_e_time,sample_rlx_j_time)
		!read(*,*)
		tgw(i)=sample_tgw_time/2/pi
		gwf(i)=sample_tgw_time/min(sample_rlx_e_time,sample_rlx_j_time)
		see(i)=dee
		sjj(i)=djj
	end do

end subroutine

subroutine save_emris_sts_hdf5(em, agw,egw,dwe,dwj, group_id, mmin_in,mmax_in,uselog)
	use md_event_datas
	use md_hdf5
	use com_sts_type
	use model_basic
	implicit none
	type(event_data_series)::em
	integer(HID_T) group_id
	! character*(*) tablename
	type(hdf5_table_type)::h5table
	integer::nfields
	real(8) agw(em%n),egw(em%n),mass(em%n),weight(em%n),dwe(em%n),dwj(em%n)
	type(s1d_hst_type)::fegw, fm, finc
	type(s2d_hst_type)::fme2d,fdwfm2d,fdwfe2d, fdwfm2dj,fdwfe2dj
	real(8) mmin,mmax,mmin_in,mmax_in,dwfmax,dwfmin
	real(8) emin,emax
	logical::uselog

	if(em%n>0)then
		if(uselog)then
			mmin=log10(mmin_in); mmax=log10(mmax_in)
			mass=log10(em%mass(:))
		else
			mmin=mmin_in;mmax=mmax_in
			mass=em%mass(1:em%n)
		end if
		weight=em%w(1:em%n)
		call fm%init(mmin,mmax,15,use_weight=.True.)
		call fm%set_range()
		call fm%get_hst(mass(1:em%n),weight(1:em%n),em%n)
		call fm%save_hdf5(group_id,"fmgw")

		call finc%init(0d0,pi,10,use_weight=.True.)
		call finc%set_range()
		call finc%get_hst(em%inc(1:em%n),weight(1:em%n),em%n)
		call finc%save_hdf5(group_id,"fincgw")

		emin=0d0; emax=1.0d0
		call fegw%init(emin,emax,20,use_weight=.True.)
		call fegw%set_range()
		call fegw%get_hst(egw(1:em%n),weight(1:em%n),em%n)
		!call fegw%get_hst(egw(:),em%w(:),em%n)
		call fegw%save_hdf5(group_id,"fegw")

		! emin=minval(egw); emax=maxval(egw)
		emin=-1.1d0; emax=0.0d0
		call fegw%init(emin,emax,20,use_weight=.True.)
		call fegw%set_range()
		call fegw%get_hst(log10(1-egw(1:em%n)),weight(1:em%n),em%n)
		!call fegw%get_hst(egw(:),em%w(:),em%n)
		call fegw%save_hdf5(group_id,"f1-egw")
		
		emin=0.0d0; emax=0.99d0
		call fme2d%init(30,15,emin,emax,mmin,mmax,use_weight=.true.)
		! print*, "2"
		call fme2d%set_range()
		! print*, "3"
		call fme2d%get_stats_weight(egw(1:em%n),mass(1:em%n),weight(1:em%n),em%n)
		! print*, "4"
		call fme2d%save_hdf5(group_id,"fmegw")
		
		emin=-1.1d0; emax=0.0d0
		
		call fme2d%init(20,15,emin,emax,mmin,mmax,use_weight=.true.)
		! print*, "2"
		call fme2d%set_range()
		! print*, "3"
		call fme2d%get_stats_weight(log10(1-egw(1:em%n)),mass(1:em%n),weight(1:em%n),em%n)
		! print*, "4"
		call fme2d%save_hdf5(group_id,"fm1-egw")

		dwfmin=-5;dwfmax=1.2
		call fdwfm2d%init(25,15,dwfmin,dwfmax,mmin,mmax,use_weight=.true.)
		! print*, "2"
		call fdwfm2d%set_range()
		! print*, "3"
		call fdwfm2d%get_stats_weight(log10(dwe(1:em%n)),mass(1:em%n),weight(1:em%n),em%n)
		! print*, "4"
		call fdwfm2d%save_hdf5(group_id,"fdwfm2de")

		dwfmin=-5;dwfmax=1.2
		call fdwfe2d%init(25,20,dwfmin,dwfmax,emin,emax,use_weight=.true.)
		! print*, "2"
		call fdwfe2d%set_range()
		! print*, "3"
		call fdwfe2d%get_stats_weight(log10(dwe(1:em%n)),log10(1-egw(1:em%n)),weight(1:em%n),em%n)
		! print*, "4"
		call fdwfe2d%save_hdf5(group_id,"fdwfe2de")

		dwfmin=-8;dwfmax=-2
		call fdwfm2dj%init(25,15,dwfmin,dwfmax,mmin,mmax,use_weight=.true.)
		! print*, "2"
		call fdwfm2dj%set_range()
		! print*, "3"
		call fdwfm2dj%get_stats_weight(log10(dwj(1:em%n)),mass(1:em%n),weight(1:em%n),em%n)
		! print*, "4"
		call fdwfm2dj%save_hdf5(group_id,"fdwfm2dj")

		dwfmin=-8;dwfmax=-2
		call fdwfe2dj%init(25,20,dwfmin,dwfmax,emin,emax,use_weight=.true.)
		! print*, "2"
		call fdwfe2dj%set_range()
		! print*, "3"
		call fdwfe2dj%get_stats_weight(log10(dwj(1:em%n)),log10(1-egw(1:em%n)),weight(1:em%n),em%n)
		! print*, "4"
		call fdwfe2dj%save_hdf5(group_id,"fdwfe2dj")

	end if
end subroutine