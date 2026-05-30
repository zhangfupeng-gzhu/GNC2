 
subroutine reset_create_time(sps)
	use com_main_gw
	implicit none
	type(chain_type)::sps
	type(chain_pointer_type),pointer::sp

    sp=>sps%head
    do while(associated(sp))
        if(.not.allocated(sp%ob)) then
            print*, "warnning, why it is not allocated?"
            print*, allocated(sp%ob)
            call sps%output_screen()
            stop
        end if
        select type(ca=>sp%ob)
        type is(particle_sample_type)
            if(sp%ob%exit_flag.eq.exit_normal)then
                sp%ob%simu_bgtime=sp%ob%exit_time
                !print*, "exit_time=",sp%ob%exit_time, sp%ob%m, sp%ob%source, sp%ob%exit_flag
            end if
        ! type is(sample_type)
        !     if(sp%ob%exit_flag.eq.exit_normal)then
        !         sp%ob%simu_bgtime=sp%ob%exit_time
        !     end if
        !     !print*, "exit_time=",sp%ob%exit_time, sp%ob%m, sp%ob%source, sp%ob%exit_flag
        end select
        sp=>sp%next
    end do
    !read(*,*)
end subroutine
   
subroutine set_clone_weight_arr(sms_single)
	use com_main_gw
	implicit none
	type(particle_samples_arr_type)::sms_single
	call sams_get_weight_clone_single(sms_single)
end subroutine
subroutine set_real_weight_arr(sms_single)
	use com_main_gw
	implicit none
	integer i
	type(particle_samples_arr_type)::sms_single
    real(8) weight_tot
    weight_tot=0
	do i=1, sms_single%n
        
        call get_sample_weight_real(sms_single%sp(i))
        if(sms_single%sp(i)%obtype.eq.star_Type_ms)then
            weight_tot=weight_tot+sms_single%sp(i)%weight_real
        end if
        ! print*, "sample%weight_real=",sms_single%sp(i)%weight_real
	end do
    
end subroutine

subroutine set_real_weight_arr_single(sms_single)
	use com_main_gw
	implicit none
	integer i
	type(particle_samples_arr_type)::sms_single

	do i=1, sms_single%n
        !sms_single%sp(i)%weight_asym=dms%weight_asym
        call get_sample_weight_real(sms_single%sp(i))
	end do
end subroutine

 
subroutine set_real_weight(sms)
	use com_main_gw
	implicit none
	integer i
	type(chain_type)::sms
	type(chain_pointer_type),pointer::pt

    pt=>sms%head
    do while(associated(pt).and.allocated(pt%ob))
        !pt%ob%weight_asym=dms%weight_asym
        call get_sample_weight_real(pt%ob)
        pt=>pt%next
    end do
	
end subroutine 

subroutine arr_to_chain_single(bkarr, chain)
    use com_main_gw
    implicit none
    type(chain_type)::chain   
    type(chain_pointer_type),pointer::pt
    integer i
    type(particle_samples_arr_type) ::bkarr
    
    call chain%init(bkarr%n)
    !print*, "arr_to_chain_single:init finished"
    pt=>chain%head
    do i=1, bkarr%n
        allocate(particle_sample_type::pt%ob)
        select type(ca=>pt%ob)
        type is(particle_sample_type)
            ca=bkarr%sp(i)
        end select
        pt%idx=i
        pt=>pt%next
    end do
end subroutine

subroutine convert_sams()
	use com_main_gw
	implicit none
    integer n
    !call bksams%get_length(n)
    !print*, "original chain length=", n
    call refine_chain(bksams)
    !call bksams%get_length(n)
    !print*, "after refine chain length=", n
    !read(*,*)
    call convert_sams_pointer_arr(bksams, bksams_pointer_arr,type=1)
    call reset_create_time(bksams)
    
end subroutine
   
subroutine get_collection_memory_usage(smsa,n)
    use com_main_gw
    implicit none
    integer n,i, nsize_smsa 
    type(particle_samples_arr_type)::smsa(n)
    nsize_smsa=0; 
    do i=1, n
        nsize_smsa=nsize_smsa+sizeof(smsa(i)%sp)/1024 
    end do
    nsize_tot=nsize_tot+nsize_smsa 
end subroutine
subroutine show_memory_usage()
		use com_main_gw
		implicit none
		integer i, j,  nsize_other

		!call get_memo_usage(proc_id)
        nsize_chain_bk=bksams%head%get_sizeof()
		nsize_chain_bk=nsize_chain_bk/1024
		nsize_arr_bk=sizeof(bksams_arr%sp)/1024
		nsize_arr_bk_norm=sizeof(bksams_arr_norm%sp)/1024
		nsize_arr_bk_pointer=sizeof(bksams_pointer_arr%pt)/1024
		nsize_tot_bk=nsize_arr_bk+nsize_chain_bk+nsize_arr_bk_norm+nsize_arr_bk_pointer
        
        nsize_other=nsize_other+sizeof(dms)/1024
        nsize_tot=nsize_tot_bk
    
end subroutine
subroutine show_total_memory_usage()
        use com_main_gw
        implicit none
        integer ierr
        integer nsizetot_all(ctl%ntasks), nsizetot_chain(ctl%ntasks)
        
        call mpi_gather(nsize_tot, 1, MPI_INTEGER,nsizetot_all, 1,MPI_INteger,0,&
             MPI_comm_world,ierr )
        call mpi_gather(nsize_chain_bk, 1, MPI_INTEGER, nsizetot_chain, 1, MPI_INteger, 0, &
             MPI_comm_world,ierr )
        if(rid.eq.0)then
            write(*,fmt=*) "tot memory usage:", sum(nsizetot_all)/1024, " Mb", &
                "where chain total:", sum(nsizetot_chain)/1025," Mb"
        end if
end subroutine

  
subroutine refine_chain(chain)
    use com_main_gw
    implicit none
    type(chain_type)::chain
    type(chain_pointer_type),pointer::pt,ps
    ps=>chain%head
    do while(associated(ps))
        pt=>ps%next
        !! check
        if(ps%ob%exit_flag.eq.exit_normal.and.ps%ob%en>ctl%energy_boundary)then
            print*, "error, ps%ob%exit_flag.ne.exit_boundary_min.and.ps%ob%en>ctl%energy_boundary"
            print*, "ps%ob%exit_flag,ps%ob%en,ctl%energy_boundary=",ps%ob%exit_flag,ps%ob%en,ctl%energy_boundary
            print*, "x,xmin=",ps%ob%x,ps%ob%en/ctl%energy0,ctl%x_boundary
            call ps%ob%print("in refine_chain")
            stop
        end if
        if(ps%ob%exit_flag.ne.exit_normal)then
            if(allocated(ps%ob))then
                deallocate(ps%ob)
            end if
            if(.not.associated(ps%prev))then
                !ps point to the head of the chain
                print*, "ps point to the head of the chain"
                if(associated(pt))then
                !call chain%output_screen()
                    chain%head=>pt
                    if((.not.associated(pt%prev))) then
                        print*, "??? not associted pt%prev", associated(pt%prev)
                        stop
                    end if
                    pt%prev=>null()
                    call destroy_attach_pointer_chain_type(ps)
                    call pt%set_head()
                else
                    chain%head=>null()
                end if
            else
                call chain_pointer_delete_item_chain_type(ps)
            end if
        end if
        ps=>pt
    end do
    !pt=>null()
    !ps=>null()
end subroutine
subroutine particle_sample_get_weight_clone(en, clone,  amplifier, e0, weight_clone,nlvl)
    use com_main_gw
	implicit none
	real(8) e0, en
	integer nlvl, clone
	integer amplifier
    real(8) weight_clone,get_clone_deep
    nlvl=0
	if(clone.ge.1)then
		if(en/e0<0)then
			!sp%weight=1d0
			!obidx=sp%obidx
			weight_clone=1d0 !sp%weight_asym
		else
            !print*, "en, e0=", en, e0, en/ctl%energy0, e0/ctl%energy0
            
            nlvl=int(get_clone_deep(en,log_clone_bd_sep, e0))
            !if(nlvl>log10clone_emax) nlvl=int(log10clone_emax)-1
            !print*, log10clone_emax
            !read(*,*)
            !print*, "nlvl=", nlvl, en/e0
            if(nlvl<0) nlvl=0
			!obidx=sp%obidx
			!sp%weight_clone=dble()**(-dble(nlvl))
			weight_clone=dble(amplifier)**(-dble(nlvl))
			
            
            !if(.not.ieee_is_finite(weight_clone))then
                !print*, "weight_clone=",weight_clone, amplifier, nlvl, en, e0
                !read(*,*)
            !end if
            !print*, "nlvl, weight_clone, en/e0=", nlvl, weight_clone, en/e0,&
            !    en/ctl%energy0, e0/ctl%energy0
           ! read(*,*)
            if(ieee_is_nan(weight_clone).or..not.ieee_is_finite(weight_clone))then
                print*, "error!weight_clone=", weight_clone
                print*, "amplifier, nlvl, sp%en, e0=", amplifier, nlvl, en, e0
                stop
            end if
		end if
	else
		!obidx=sp%obidx
		weight_clone=1d0
		!sp%weight_real=sp%weight_asym
	end if
end subroutine

subroutine prepare_ini_data()
	use com_main_gw
    use md_star_pot
	implicit none
	character*(4) tmprid
    integer j,ier

	write(unit=tmprid,fmt="(I4)") rid+1+ctl%ntask_bg

    call input_chains_bin(bksams,"output/ini/bin/single/samchn"//trim(adjustl(tmprid)))
    ! print*, bksams%head%ob%exit_flag
    ! stop
    call input_diffuse_mspec_bin("output/ini/hdf5/dms_0")
    ctl%run_snap_time_0=0d0
	if(rid.eq.0)then
    	print*, "readin dms finished!"
	end if
    !print*, "before run....................."
    !call show_memory_usage()
    !call show_total_memory_usage()
    call update_arrays_single(.true.) 
    if(rid.eq.1)then
	!	print*, bksams_arr_norm%sp(1:10)%x
        print*, "emin,emax=",emin_factor,emax_factor
        print*, "logrmin,logrmax=",dms%logrmin, dms%logrmax
!        print*, "emindst,emaxdst=", emin_dstr_factor,emax_dstr_factor
	end if
    
    
end subroutine  
subroutine prepare_data(rsid,snapid)
    use com_main_gw
	implicit none
	character*(6) tmprid, tmpsnapid,tmpsupdate_max
    character*(200) str_
    logical ex
    integer rsid, snapid
    integer j, ier

	write(unit=tmprid,fmt="(I6)") rsid
    write(unit=tmpsnapid,fmt="(I6)") snapid

   
    str_=trim(adjustl(tmprid))//"_"//trim(adjustl(tmpsnapid))
    print*, "str=", trim(adjustl(str_))

    inquire(file="output/ecev/bin/single/samchn"//trim(adjustl(str_))//".bin",exist=ex)
    
    if(.not.ex)then
        print*, "output/ecev/bin/single/samchn"//trim(adjustl(str_))//".bin"," not exist"
        stop
    endif

    call input_chains_bin(bksams,"output/ecev/bin/single/samchn"//trim(adjustl(str_)))
    
    call convert_sams()	
    !write(unit=tmpsupdate_max,fmt="(I4)") ctl%i_spshot_bg_last
    call input_diffuse_mspec_bin("output/ecev/dms/dms_"//trim(adjustl(tmpsnapid)))
    print*, "readin dms finished!"

    ctl%run_snap_time_0=ctl%run_snap_time_f
    
    ! do j=1, ctl%ntasks
    !     if(rid.eq.j-1)then
    !         print*, "rid=",rid
    !         print*, "rho_rmax,rmin=", dms%mb(7)%dsp(2)%p%spt_rho_rmax,&
    !             dms%mb(7)%dsp(2)%p%spt_rho_rmin
    !     end if
    !     call mpi_barrier(mpi_comm_world,ier)
    ! end do
    ! stop

    !print*, "the start time of append mode is not correct, stoped"
    !stop

    !print*, "before run....................."
    !call show_memory_usage()
    !call show_total_memory_usage()
    call update_arrays_single(.true.) 
   
end subroutine
subroutine prepare_ini_data_append()
	use com_main_gw
	implicit none
	character*(6) tmprid, tmpsnapid,tmpsupdate_max
    character*(200) str_
    logical ex
    call prepare_data(rid+1+ctl%ntask_bg,ctl%n_spshot_bg)
	
end subroutine

 subroutine get_lambda(lambda)
    use com_main_gw
    use md_star_pot
    implicit none
    real(8) lambda
    real(8),parameter::gamma0=0.4
    !print*, "m0_cl,spp_new%N_r_within_max,gamma0", m0_cl, spp_new%N_r_within_max,gamma0, spp_new%M_r_within_max
    lambda=log(m0_cl*spp_new%N_r_within_max*gamma0)
    !print*, "lambda=",lambda
   ! read(*,*)
 end subroutine
