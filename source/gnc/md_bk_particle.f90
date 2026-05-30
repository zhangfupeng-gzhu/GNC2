module md_bk_species
    use md_coeff
    use com_sts_type
    use md_particle_sample
	!use md_chain_pointer
    implicit none
    
    type particle_samples_arr_type
        integer n 
        type(particle_sample_type),allocatable::sp(:)
        contains
			procedure::init=>init_particle_sample_arr
			!procedure::add_member=>add_member_particle_sample_arr
			!procedure::output_bin=>output_particle_sams_arr_bin
			!procedure::input_bin=>input_particle_sams_arr_bin
			procedure::select=>sams_arr_select_single
			procedure::output_txt=>output_particle_sams_txt
			procedure::deallocation=>deallocation_particle_sams
			!procedure::get_int_total
			!procedure::get_real_total
			!procedure::convert_to_int_arrays
			!procedure::convert_to_real_arrays
    end type	
   
   private::init_particle_sample_arr
  ! private::input_particle_sams_arr_bin, output_particle_sams_arr_bin
   private::sams_arr_select_single,sams_selection_function
   private::output_particle_sams_txt,deallocation_particle_sams
contains
subroutine deallocation_particle_sams(bksps)
	implicit none
	class(particle_samples_arr_type)::bksps
	if(allocated(bksps%sp)) then
	!	print*, "deallocating arr...", bksps%n
		deallocate(bksps%sp)
	!	print*, "deallocate arr success!"
	end if
end subroutine
subroutine init_particle_sample_arr(bksps,n)
    implicit none
    integer n,i
    class(particle_samples_arr_type)::bksps
    !print*, "start arr ini", n
    !print*, associated(bksps%sp)
    if(allocated(bksps%sp)) then
    !	print*, "deallocating arr...", bksps%n
        deallocate(bksps%sp)
    !	print*, "deallocate arr success!"
    end if
    !print*, "start allocation"
    allocate(bksps%sp(n))
    bksps%n=n
    !print*, "allocate success"
    do i=1, n
        call bksps%sp(i)%init()
    end do
    !print*, "finished init"
end subroutine
    subroutine sams_arr_select_condition_single(sps, sps_out,selection_func,ipar,rpar)
		implicit none	
		type(particle_samples_arr_type)::sps,sps_out
		integer nsel,i, exitflag, nhiar
		real(8) timebg,timeed
		logical::selection_func
		logical::selected
		integer,optional:: ipar(10)
		real(8),optional:: rpar(10)
		nsel=0
		do i=1, sps%n
			if(present(ipar))then
				selected=selection_func(sps%sp(i),ipar,rpar)
			else
				selected=selection_func(sps%sp(i))
			end if
			if(selected)	 nsel=nsel+1
		end do
		call sps_out%init( nsel)
		nsel=0
		do i=1, sps%n
			if(present(ipar))then
				selected=selection_func(sps%sp(i),ipar,rpar)
			else
				selected=selection_func(sps%sp(i))
			end if
			if(selected) then
				nsel=nsel+1
	!			PRINT*, "NSEL=",NSEL
				sps_out%sp(nsel)=sps%sp(i)
	!			PRINT*, SPS_OUT%SP(NSEL)%EXIT_FLAG, sps_out%sp(nsel)%agw, sps%sp(nsel)%agw
			end if		
		end do
	end subroutine
	subroutine add_member_particle_sample_arr(sps,idx, sp, sps_tmp)
		!insert after idx
		implicit none	
		type(particle_samples_arr_type)::sps
		type(particle_samples_arr_type)::sps_tmp
		type(particle_sample_type)::sp
		integer i, idx
		
		call sps_tmp%init(sps%n+1)
		do i=1, idx
			sps_tmp%sp(i)=sps%sp(i)
		end do
		sps_tmp%sp(idx)=sp
		do i=idx+1, sps%n+1
			sps_tmp%sp(i)=sps%sp(i-1)
		end do
		call copy_particle_sample_arr(sps_tmp, sps)
	end subroutine
	subroutine copy_particle_sample_arr(scopy, sp)
		!scopy=>sp
		implicit none
		type(particle_samples_arr_type)::scopy, sp
		integer i
		call sp%init(scopy%n)
		do i=1, scopy%n
			sp%sp(i)=scopy%sp(i)
		end do
	end subroutine
    subroutine smmerge_arr_single(sma,n,smam)
		implicit none
		integer n,i,j, nsam
		type(particle_samples_arr_type)::sma(n), smam
		nsam=0
		do i=1, n
			nsam=nsam+sma(i)%n
		end do
		!print*, "nsam=",nsam
		call smam%init(nsam)
		nsam=0
		do i=1, n
			do j=1, sma(i)%n
				nsam=nsam+1
				!print*, nsam
				smam%sp(nsam)=sma(i)%sp(j)
			end do
		end do
	end subroutine

	subroutine set_sample_arr_indexs_rid_particle(sams_arr,rid)
		!use com_main_gw
		implicit none
		type(particle_samples_arr_type)::sams_arr
		integer i,rid
		do i=1, sams_arr%n
			sams_arr%sp(i)%idx=i
			sams_arr%sp(i)%rid=rid
		end do
	end subroutine

   
	subroutine sams_arr_select_single(sps, sps_out, exitflag, source, timebg, timeed)
		implicit none
		class(particle_samples_arr_type)::sps, sps_out
		integer nsel,i, exitflag, source
		real(8) timebg,timeed
		!logical sams_selection_function
		nsel=0
		do i=1, sps%n
			if(sams_selection_function(sps,i,exitflag,source,timebg,timeed)) nsel=nsel+1
		end do
		!print*, "nsel=",nsel
		call sps_out%init(nsel)
		nsel=0
		do i=1, sps%n
			if(sams_selection_function(sps,i,exitflag,source,timebg,timeed)) then
				nsel=nsel+1
				sps_out%sp(nsel)=sps%sp(i)
			end if		
		end do
	!contains

	end subroutine
	logical function sams_selection_function(sps,i,exitflag,source,timebg,timeed)
		implicit none
		class(particle_samples_arr_type)::sps
		integer i,exitflag,source
		real(8) timebg,timeed
		sams_selection_function=.false.
		if(ieee_is_nan(sps%sp(i)%weight_clone).or.ieee_is_nan(sps%sp(i)%weight_N))then
			print*, "selection error:", sps%sp(i)%weight_clone, sps%sp(i)%weight_N, &
			sps%sp(i)%id
			stop
		end if
		!if(sps%sp(i)%id.eq.98172504)then
		!	print*, "98172504%exit_flag=",sps%sp(i)%exit_flag
		!end if
		!if(abs(sps%sp(i)%weight_real-sps%sp(i)%weight_asym*sps%sp(i)%weight_n*&
		!	sps%sp(i)%weight_clone)>0.1)then
		!	print*, "error!, id=", sps%sp(i)%id, sps%sp(i)%weight_real, &
		!	sps%sp(i)%weight_asym, sps%sp(i)%weight_n,  &
		!	sps%sp(i)%weight_clone
		!	stop
		!end if
		! print*, sps%sp(i)%exit_flag, exitflag
		if(sps%sp(i)%exit_flag.eq.exitflag.or.(exitflag.eq.-1) &
			.or.(exitflag.eq.-2.and.sps%sp(i)%exit_flag.ne.exit_invtransit))then
			if((sps%sp(i)%source.eq.source).or.(source.eq.-1))then
				if(sps%sp(i)%exit_time>timebg.or.timebg<0)then
					if(sps%sp(i)%exit_time<timeed.or.timeed<0)then
						sams_selection_function=.true.
					end if
				end if
			end if
		end if
	end function

	subroutine output_particle_sams_txt(sps, fl)
		implicit none
		character*(*) fl
		integer i,n
		class(particle_samples_arr_type)::sps
		real(8) x,y,z
        real(8),external::rnd
        
		!type(particle_sample_type),pointer::sp
		open(unit=999,file=trim(adjustl(fl))//".txt")
		write(unit=999,fmt="(13A20, 3A10)") "m", "aout", "eout", "inc","om","pe", "rp", "en", "jm","weight", "x","y","z", &
			"exit_flag", "source", "obtype"
		do i=1, sps%n
			associate(sp=>sps%sp(i))		
				sp%byot%an_in_mode=an_in_mode_mean
                sp%byot%me=rnd(0d0,2*pi)
				call by_em2st(sp%byot)
				call by_split_from_rd(sp%byot)
				write(unit=999,fmt="(13E20.10, 3I10)") sp%m, sp%byot%a_bin,  sp%byot%e_bin, &
			sp%byot%inc, sp%byot%om, sp%byot%pe ,sp%rp, sp%en, sp%jm,sp%weight_real,sp%byot%rd%x,&
                 sp%exit_flag, sp%source, sp%obtype
			end associate
		end do
		close(unit=999)
	end subroutine

end module