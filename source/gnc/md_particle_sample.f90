module md_particle_sample
    use,intrinsic::ieee_arithmetic
	use md_stellar_history
    use md_binary 
    real(8)::log10clone_emax
    type track_type
        real(8) time, r_lc, rp, ain, ein, x, jm, ac,ec
        real(8) Incin, incout, omin, omout, meout
        integer state_flag, ms_star_type, mm_star_type
    end type
    type particle_sample_type
		integer id                    ! a unque id for the particle
		integer rid, idx  !rid=the proc id, idx=the index in array
		integer obtype, obidx
		integer state_flag_last         ! flag of the last dynamic, use to determine the type of exit_flag
		integer exit_flag
		integer length, length_to_expand
		integer source, track_step
		integer write_down_track, within_jt
		integer N_gene  ! number of generation
		integer state_emri_last, state_emri_current, Lvl_clone ! the level of clone
        type(type_stellar_history)::sh
		real(8) r_lc, m, en0, jm0
		!type(particle)::p
		real(8) period, rp, tgw, simu_bgtime
		real(8) En, x, create_time,exit_time
		real(8) jph, jm, jc,ra, raq
		!real(8) djp, elp, den, djp0
		real(8) weight_clone  ! the weight factor due to clone
		real(8) weight_N      ! the weight factor due to number of particles set by simulation
							  ! weight_N=1 by default, and can be changed due to collisions.
		!real(8) weight_asym   ! the weight factor due to asymptotic boundary condition 
		real(8) weight_real   ! the real weight used for calculation
		type(binary) byot, byot_ini, byot_bf

		type(track_type),allocatable:: track(:)
		contains
		procedure::track_init
		procedure::init=>particle_sample_init
		procedure::read_info=>read_sample_info
		procedure::write_info=>write_sample_info
		procedure::print=>print_particle_sample
		!procedure::reset_to_ini=>reset_to_ini_particle
	end type
    integer,parameter::track_length_expand_block=100
	integer,parameter::nint_particle=12, nreal_particle=19
	integer,parameter::flag_ini_ini=1
    integer,parameter::flag_ini_bd=2
    integer,parameter::flag_ini_or=5
    integer,parameter::record_track_nes=1, record_track_detail=2, record_track_all=3
	integer,parameter::source_bk=2
	integer,parameter::state_ae_evl=1
	integer,parameter::state_emax=19, state_plunge=199,state_emri=198
	integer,parameter::state_td=71
	
    integer,parameter::exit_normal=1
    integer,parameter::exit_other=100   
	integer,parameter::exit_stellar_evl_supnov=105     
	integer,parameter::exit_emri_single=180, exit_lc=181 
    
    integer,parameter::exit_tidal=2
    integer,parameter::exit_by_ignore=20
    integer,parameter::exit_merge_eject=200
    integer,parameter::exit_max_reach=3
	integer,parameter::exit_boundary_min=4, exit_boundary_max=5
	integer,parameter::exit_invtransit=7, exit_tidal_empty=8, exit_tidal_full=9 
	
	integer,parameter:: star_type_BH=1,star_type_NS=2,star_type_MS=3, star_type_WD=4
	integer,parameter:: star_type_BD=5  ! brown dwarf 
	integer,parameter:: star_type_RG=6
	integer,parameter:: star_type_dark_matter=7
	integer,parameter:: star_type_nakedHe=8
	integer,parameter:: star_type_other=99

	private::particle_sample_init !, write_sample_info, read_sample_info
	private::print_particle_sample!, particle_sample_get_weight_clone
	!private::reset_to_ini_particle

contains
subroutine print_particle_sample(this, str)
	implicit none
	class(particle_sample_type)::this
	character*(*) str
	character*(9) typestr
	print*, str
	write(*,fmt="(A20,4F15.6)") "en=", this%en
	write(*,fmt="(A20,4F15.6)") "x,j=", this%x,this%jm
	write(*,fmt="(A20,4I15)") "within_jt=", this%within_jt
	write(*,fmt="(A20,4I15)") "rid,idx,length,source=", & 
			this%rid, this%idx, this%length, this%source
	write(*,fmt="(A20,4I15)") "id, exit_flag,N_gene=", & 
			this%id,this%EXIT_FLAG, this%n_GENE
    write(*,fmt="(A20,4I15)") "state_last=", this%state_flag_last
	call get_star_type(this%obtype,typestr)
	write(*,fmt="(A20,A15)") "obtype=",this%obtype,typestr
	!write(*,fmt="(A20,4I15)") "state_last=", this%state_flag_last
						
	write(*,fmt="(A20,4F15.6)") "otby:a,e, I, Om=", & 
			this%byot%a_bin, this%byot%e_bin, this%byot%Inc, this%byot%Om
	write(*,fmt="(A30,6E15.6)") "ctime, w(real,  clone, n)=", & 
			this%create_time, this%weight_real, this%weight_clone, &
            this%weight_n
end subroutine
	subroutine track_init(sp,n)
		implicit none
		class(particle_sample_type)::sp
		integer n
		if(allocated(sp%track))then
		deallocate(sp%track)
		endif
		allocate(sp%track(n))
		sp%length=0
		sp%length_to_expand=n
		!print*, "track init", n
	end subroutine
	subroutine particle_sample_init(sp)
		implicit none
		class(particle_sample_type) ::sp
		sp%write_down_track=0
		sp%source=0;sp%within_jt=0
		sp%byot%a_bin=0; sp%byot%e_bin=0;
		sp%track_step=1; sp%period=0d0
		sp%exit_flag=0;sp%m=0; sp%weight_real=-1d99
		sp%weight_clone=-1d99;sp%weight_N=-1d99
        sp%exit_time=0d0
		!if(allocated(sp%track)) deallocate(sp%track)
		!print*,"init:", sizeof(sp%track), allocated(sp%track)
		call track_init(sp,0)
		sp%sh%n=0
	end subroutine

subroutine read_sample_info(sp,funit,version_number)
	integer funit
	integer mypos,version_number
	class(particle_sample_type)::sp
    logical,save::first=.true.

	select case(version_number) 
	case(107)
		call read_sample_info_107(sp,funit)
	end select
end subroutine
subroutine read_sample_info_107(sp,funit)
	integer funit
	integer mypos,version_number
	class(particle_sample_type)::sp
    logical,save::first=.true.
	
	call sp%init()
!	print*, sp%byot%a_bin
	read(unit=funit) sp%exit_time,sp%r_lc,sp%m ,sp%en0,sp%jm0, &
                   sp%period,sp%rp,sp%tgw, sp%state_flag_last, sp%n_gene,sp%state_emri_last,sp%state_emri_current,sp%Lvl_clone
	
	read(unit=funit) sp%obtype, sp%obidx, sp%rid, sp%idx, sp%id 
	read(unit=funit) sp%byot,sp%byot_ini, sp%byot_bf 
	read(unit=funit) sp%length_to_expand,sp%exit_flag,sp%create_time,sp%En, sp%x, &
				   sp%jph,sp%Jm, sp%rp, sp%ra, sp%period, sp%jc,sp%raq, &
	               sp%write_down_track,sp%track_step,sp%source , &
                    sp%within_jt, sp%length, &
					 sp%weight_real, sp%weight_N,sp%weight_clone
	if(ieee_is_nan(sp%weight_real))then
		print*, "error! sp%weight_real = NaN", &
			sp%weight_real, sp%weight_N
		stop
	end if 
 	if(sp%length>0)Then 
		if(allocated(sp%track))deallocate(sp%track)
		allocate(sp%track(sp%length))
		read(unit=funit) sp%track(1:sp%length) 
	end if

	call sp%sh%read_info(funit)
	
end subroutine
 
character*(8) function star_type_str(type_int)
	implicit none
	integer type_int

	call get_star_type(type_int,star_type_str)

end function
subroutine write_sample_info(sp,funit,version_number)
	implicit none
	integer funit,version_number
	class(particle_sample_type)::sp
    logical,save::first=.true.
	select case(version_number)
	case(107)
		call write_sample_info_107(sp,funit)
	end select
end subroutine
 

subroutine write_sample_info_107(sp,funit)
	implicit none
	integer funit
	class(particle_sample_type)::sp
    logical,save::first=.true.

	write(unit=funit) sp%exit_time,sp%r_lc,sp%m ,sp%en0,sp%jm0, &
                   sp%period,sp%rp,sp%tgw, sp%state_flag_last, sp%n_gene,sp%state_emri_last,sp%state_emri_current,sp%Lvl_clone
	
	write(unit=funit) sp%obtype, sp%obidx, sp%rid, sp%idx, sp%id

	write(unit=funit) sp%byot,sp%byot_ini, sp%byot_bf 
	write(unit=funit)  sp%length_to_expand,sp%exit_flag,sp%create_time,sp%En, sp%x, &
						sp%jph,sp%Jm, sp%rp, sp%ra, sp%period, sp%jc,sp%raq, &
						sp%write_down_track,sp%track_step,sp%source , &
						sp%within_jt, sp%length, &
						sp%weight_real, sp%weight_N,sp%weight_clone	
	if(ieee_is_nan(sp%weight_real))then
		print*, "error! sp%weight_real = NaN",  &
			sp%weight_real, sp%weight_N
		stop
	end if 
 	if(sp%length>0)Then 
		print*, "sp%length=",sp%length
		write(unit=funit) sp%track(1:sp%length)
	end if

	call sp%sh%write_info(funit) 
end subroutine



end module