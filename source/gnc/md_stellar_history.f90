module md_stellar_history
    use constant
    type type_stellar_history
        integer n
        integer cur_idx
        integer,allocatable:: ktype(:)
        integer,allocatable:: kwtype(:)
        real(8),allocatable::mass(:)
        real(8),allocatable::radius(:)
        real(8),allocatable::time(:)
        contains
            procedure::init=>init_stellar_history
            procedure::write_info=>write_info_stellar_history
            procedure::read_info=>read_info_stellar_history
            procedure::print=>print_stellar_history
            procedure::get_history=>get_history_stellar_history!,get_history_stellar_history_brown_dwarf
    end type
    integer,parameter::nrecord_max=2000
    real(8) brown_data_age(4)
    real(8) brown_data_xs(14,4), brown_data_ys(14,4), brown_data_y2(14,4)

    private::init_stellar_history,print_stellar_history,get_history_stellar_history
    private::get_history_stellar_history_brown_dwarf
    
contains
    subroutine init_stellar_history(this, n)
        implicit none
        class(type_stellar_history)::this
        integer n
        if(allocated(this%mass))then
            deallocate(this%mass,this%time,this%ktype,  this%radius)
        end if
        allocate(this%mass(n),this%time(n),this%ktype(n), this%radius(n))
        this%n=n

        if(allocated(this%kwtype))then
            deallocate(this%kwtype)
        end if
        allocate(this%kwtype(n))
        this%cur_idx=0
    end subroutine
    subroutine print_stellar_history(this)
        use md_mobse_stellar_single,only:get_kstar_type,get_kw_type
        class(type_stellar_history)::this
        integer i
        character*(9) kwtype,kstartype
        
        write(*,fmt="(8A15)") "time", "mass","radius(rsun)", "type", "event"
        do i=1, this%n
            kstartype=get_kstar_type(this%ktype(i))
            kwtype=get_kw_type(this%kwtype(i) )
            write(*,fmt="(3F15.3, A15, A15)") this%time(i), this%mass(i), this%radius(i)/rd_sun, &
                 kstartype,kwtype 
        end do
    end subroutine    
    subroutine write_info_stellar_history(this, funit)
        implicit none
        class(type_stellar_history)::this
        integer funit
        write(unit=funit) this%n
		if(this%n<1) return
        write(unit=funit) this%cur_idx
        write(unit=funit) this%ktype(1:this%n), &
            this%mass(1:this%n), this%time(1:this%n), &
             this%radius(1:this%n), this%kwtype(1:this%n)
    end subroutine
    subroutine read_info_stellar_history(this, funit)
        implicit none
        class(type_stellar_history)::this
        integer funit
        read(unit=funit) this%n
        !print*, "this%n,cur_idx=",this%n,this%cur_idx
        call this%init(this%n)
		if(this%n<1) return
        read(unit=funit) this%cur_idx
        read(unit=funit) this%ktype(1:this%n), &
        this%mass(1:this%n), this%time(1:this%n) ,&
        this%radius(1:this%n), this%kwtype(1:this%n)
    end subroutine
    subroutine get_history_stellar_history_ms(this,kstar0,m, z, time_tot,time_start,output_flag)
        implicit none
        class(type_stellar_history)::this
        real(8) m, z, time_tot,time_start
        integer nseq,output_flag
        real(8) time_seq(nrecord_max), mass_seq(nrecord_max),radius_seq(nrecord_max)
        integer kstar_seq(nrecord_max), kw_seq(nrecord_max),kstar0
        if(m>0.1d0)then
#ifdef  MOBSE
            call get_evl_single_mobse(m, time_tot,kstar0,z,100,time_seq,&
            kstar_seq, mass_seq, kw_seq, radius_seq, nseq, nrecord_max, .true., output_flag)
#endif
            call this%init(nseq)
        else
            print*, "error!, this subroutine works only for m>0.1"
            stop
        end if

        this%ktype(1:nseq)=kstar_seq(1:nseq)
        this%kwtype(1:nseq)=kw_seq(1:nseq)
        this%mass(1:nseq)=mass_seq(1:nseq)
        this%time(1:nseq)=time_seq(1:nseq)+time_start
        this%radius(1:nseq)=10**radius_seq(1:nseq)*rd_sun
    end subroutine
    subroutine get_history_stellar_history(this,kstar0,m, z, time_tot,time_start,output_flag)
        implicit none
        class(type_stellar_history)::this
        real(8) m, z, time_tot,time_start
        integer kstar0,output_flag
        if(m>0.1d0)then
            call get_history_stellar_history_ms(this,kstar0,m, z, time_tot,time_start,output_flag)
        else
            call get_history_stellar_history_brown_dwarf(this,kstar0,m, z, time_tot,time_start)
        end if
    end subroutine
    subroutine get_history_stellar_history_brown_dwarf(this,kstar0,m, z, time_tot,time_start)
        implicit none
        class(type_stellar_history)::this
        real(8) m, z, time_tot
        integer nseq,i
        real(8) time_seq(100), mass_seq(100),radius_seq(100),time_start
        integer kstar_seq(100), kw_seq(100),kstar0
        nseq=4
        if(m<=0.1d0)then
            do i=1, 4
                time_seq(i)=brown_data_age(i)
                call get_brown_dwarf_radius(m,brown_data_xs(1:14,i),brown_data_ys(1:14,i),&
                    brown_data_y2(1:14,i),radius_seq(i))
                mass_seq(i)=m
                kstar_seq(i)=kstar0
                kw_seq(i)=2
            end do
            kw_seq(1)=1
        else
            print*, "error!, this subroutine works only for m<=0.1"
            stop
        end if
        call this%init(nseq)
        this%ktype(1:nseq)=kstar_seq(1:nseq)
        this%kwtype(1:nseq)=kw_seq(1:nseq)
        this%mass(1:nseq)=mass_seq(1:nseq)
        this%time(1:nseq)=time_seq(1:nseq)+time_start
        this%radius(1:nseq)=radius_seq(1:nseq)*rd_sun
    end subroutine
    
    subroutine get_current_stellar_info(this,curtime, ktype,kwtype,mass, radius)
        implicit none
        class(type_stellar_history)::this
        real(8) curtime, mass, radius
        integer,intent(out):: ktype, kwtype
        integer i
        real(8),parameter::tiny=1d-9
        if(this%n.eq.0)then
            print*, "error stellar history not calculated"
            stop
        end if
        if(curtime<=this%time(1))then
            this%cur_idx=1
        elseif(curtime>this%time(this%n)+tiny)then
            !print*, "curtime not in rage", ctime, sh%time(1),sh%time(sh%n)
            !flag=-1
            this%cur_idx=this%n
            !return
        else

loop1:		do i=1, this%n-1
                if((this%time(i).le.curtime).and.(this%time(i+1)+tiny>=curtime))then
                    this%cur_idx=i
                    exit loop1
                end if
                this%cur_idx=this%n
            end do loop1
        end if

        ktype=this%ktype(this%cur_idx)
        kwtype=this%kwtype(this%cur_idx)
        mass=this%mass(this%cur_idx)
        radius=this%radius(this%cur_idx)
        ! if(mass<0.1d0)then
        !     print*, "mass, curtime,this%cur_idx=",mass, curtime, this%cur_idx, radius
        ! end if
    end subroutine
end module