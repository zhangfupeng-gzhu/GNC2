module md_fc_type
    use md_s1d_hst_type
    implicit none
	!Count the all kinds of frequency or numbers of a data in bins, including histogram
    type,extends(s1d_hst_basic_type):: sts_fc_type
        integer, allocatable:: cnb(:)
		real(8), allocatable:: cnbw(:), cfx(:), cfxw(:), pn(:), pnw(:)
        real(8), allocatable:: y2_y(:)

        ! cnb: unweighted cumulative number of samples
        ! cnbw:  weighted cumulative number of samples  
        ! fx:  unweighted normalized fraction of samples per pin
        ! fxw:   weighted normalized fraction of samples per pin
        ! cfx:  unweighted normalized cumulative fraction of samples 
        ! cfxw:   weighted normalized cumulative fraction of samples 
        ! pn: unweighted possion noise associated in each bin
		! pnw:  weighted possion noise associated in each bin
        integer flaglog  
        ! flaglog=fc_spacing_log :    Output the bins in log scale
        ! flaglog=fc_spacing_linear:  Output the bins in linear scale
        !integer nbin
        !nbin: the total number of bins
        integer type_int_size
		integer type_real_size
		integer type_log_size
		contains
		    procedure::write_fc
            procedure::read_fc
            procedure::set_range=>set_range_fc
            procedure::print=>print_fc
            procedure::normalize=>fc_normalize
            procedure::init=>init_fc
            procedure::init_intrp_y=>init_intrp_fc_y
            procedure::output_fc
            procedure::set_fc_xb
            procedure::deallocate=>deallocate_fc
            procedure::get_value_l_y=>get_value_fc_linear_y
            procedure::get_value_s_y=>get_value_fc_spline_y
			!procedure::get_value_d=>get_value_fc_direct
            procedure::sort=>sort_fc
#ifdef USE_HDF5                
            procedure::save_hdf5=>save_fc_hdf5
            procedure::read_hdf5=>read_fc_hdf5
#endif            
            !procedure::declare_fc_mpi_type
			generic :: write(unformatted) => write_fc
		  	generic :: read(unformatted) => read_fc
    end type
    integer,parameter::fc_spacing_log=1, fc_spacing_linear=0
    private::print_fc, write_fc, read_fc, fc_normalize,deallocate_fc
	!private::get_value_fc_direct
    private::get_value_fc_linear_y,get_value_fc_spline_y
    private::sort_fc, init_fc
#ifdef USE_HDF5    
    private::save_fc_hdf5,read_fc_hdf5
#endif    
contains
    subroutine print_fc(this,str_,print_y2_in)
        class(sts_fc_type)::this
        character*(15) str_fc_spacing
        character*(*) , optional::str_
        logical,optional::print_y2_in
        logical print_y2
        integer i

        if(present(str_))then
            write(*,*) "for fc=", trim(adjustl(str_))
        end if
        if(this%use_weight)then
            write(unit=*, fmt="(10A15)") "spacing", "nbin", "xmin", "xmax", "xstep", "nsw", "ns"
            select case (this%flaglog)
            case (fc_spacing_log)
                str_fc_spacing="LOG"
            case (fc_spacing_linear)
                str_fc_spacing="LINEAR"
            end select
            write(unit=*, fmt="(A15, I15, 4E15.5, I15)") str_fc_spacing, this%nbin, &
                this%xmin, this%xmax, this%xstep, this%nsw, this%ns
            write(unit=*, fmt="(20A20)") "X", "FX", "FXW", "NB", "NBW", "cfx", "cfxw"
            do i=1, this%nbin
                write(unit=*, fmt="(3E20.10, I20, 3E20.10)") this%xb(i), this%fx(i), this%fxw(i), this%nb(i),&
                 this%nbw(i), this%cfx(i), this%cfxw(i)
            end do
        else
            write(unit=*, fmt="(10A15)") "spacing", "nbin", "xmin", "xmax", "xstep",  "ns"
            select case (this%flaglog)
            case (fc_spacing_log)
                str_fc_spacing="LOG"
            case (fc_spacing_linear)
                str_fc_spacing="LINEAR"
            end select
            write(unit=*, fmt="(A15, I15, 3E15.5, I15)") str_fc_spacing, this%nbin, &
                this%xmin, this%xmax, this%xstep, this%ns
            write(unit=*, fmt="(20A20)") "X", "FX",  "NB", "cfx"
            do i=1, this%nbin
                write(unit=*, fmt="(2E20.10, I20, E20.10)") this%xb(i), this%fx(i), this%nb(i), this%cfx(i)
            end do
        end if
    end subroutine
    subroutine write_fc(fc, unit, iostat, iomsg)
       class(sts_fc_type), intent(in) :: fc
       integer, intent(in) :: unit
       integer, intent(out) :: iostat
       character(*), intent(inout) :: iomsg
       integer n

       ! Write a record giving sizes for the allocation
       write(unit, iostat=iostat, iomsg=iomsg) fc%nbin, fc%xmin, fc%xmax, &
					fc%nsw, fc%xstep, fc%ns, fc%flaglog, fc%use_weight, fc%is_spline_prepared
	   n=fc%nbin
	   if(n<1) return
       if(fc%use_weight)then
	        write(unit, iostat=iostat, iomsg=iomsg) fc%xb(1:n), fc%fx(1:n), &
		 fc%fxw(1:n), fc%nbw(1:n), fc%cnbw(1:n), fc%cfx(1:n), fc%cfxw(1:n), &
		 fc%pn(1:n), fc%pnw(1:n), fc%nb(1:n), fc%cnb(1:n), fc%y2_y(1:n)
        else
           write(unit, iostat=iostat, iomsg=iomsg) fc%xb(1:n), fc%fx(1:n), &
		 fc%cfx(1:n), fc%pn(1:n), fc%nb(1:n), fc%cnb(1:n), fc%y2_y(1:n)
       endif
       if(fc%is_spline_prepared)then
        write(unit, iostat=iostat, iomsg=iomsg) fc%y2(1:n)
        end if
     end subroutine write_fc

    subroutine read_fc(fc, unit, iostat, iomsg)
       class(sts_fc_type), intent(inout) :: fc
       integer, intent(in) :: unit
       integer, intent(out) :: iostat
       character(*), intent(inout) :: iomsg
       integer n, flaglog, ns
	   real(8) xmin, xmax, xstep, nsw
       logical use_weight, is_spline_prepared
       ! Write a record giving sizes for the allocation
       read(unit, iostat=iostat, iomsg=iomsg) n, xmin, xmax, &
            fc%nsw, fc%xstep, fc%ns, flaglog, use_weight, fc%is_spline_prepared
        ! print*, n, xmin, xmax, &
        ! nsw, xstep, ns, flaglog, use_weight, fc%is_spline_prepared
	   if(n<1) return
	   call init_fc(fc, xmin, xmax,n, flaglog, use_weight)

	   n=fc%nbin

       if(use_weight)then
        read(unit, iostat=iostat, iomsg=iomsg) fc%xb(1:n), fc%fx(1:n), &
            fc%fxw(1:n), fc%nbw(1:n), fc%cnbw(1:n), fc%cfx(1:n), fc%cfxw(1:n), &
            fc%pn(1:n), fc%pnw(1:n), fc%nb(1:n), fc%cnb(1:n), fc%y2_y(1:n)
       else
        read(unit, iostat=iostat, iomsg=iomsg) fc%xb(1:n), fc%fx(1:n), &
		 fc%cfx(1:n), fc%pn(1:n), fc%nb(1:n), fc%cnb(1:n), fc%y2_y(1:n)
       end if     
       if(fc%is_spline_prepared)then
        if(.not.allocated(fc%y2))allocate(fc%y2(n))
        read(unit, iostat=iostat, iomsg=iomsg) fc%y2(1:n)
     end if 
     end subroutine read_fc
     subroutine deallocate_fc(this)
        implicit none
        class(sts_fc_type)::this
        call this%s1d_basic_type%deallocate()
        if(allocated(this%cfx)) then
            deallocate(this%cfx,  this%nb, this%cnb, &
            this%pn, this%y2_y)
            if(this%use_weight)then
                deallocate(this%nbw, this%fxw, this%cfxw, this%cnbw, this%pnw)
            end if
        endif
     end subroutine
    subroutine init_fc(this, xmin,xmax, n, flag_log,use_weight)
        implicit none
		! Initialize the data type of fc
        class(sts_fc_type)::this
        logical,optional:: use_weight
        integer n, flag_log
        real(8) xmin,xmax
	
       !call this%deallocate()
        if(n<1)then
            print*, "error: n<1 in init_fc, n=", n
            stop
        end if
        if(present(use_weight))then
            call init_s1_hst_basic(this, xmin,xmax,n, use_weight)
        else
            call init_s1_hst_basic(this, xmin,xmax,n, .false.)
        end if
        if(allocated(this%cfx))then
            deallocate(this%cfx, this%cnb,this%y2_y, this%pn)
        end if

        allocate(this%cfx(n))
        allocate(this%cnb(n))
        allocate(this%y2_y(n), this%pn(n))
        this%cnb=0; 
        this%cfx=0;  this%pn=0 ; this%y2_y=0
        if(this%use_weight)then
            if(allocated(this%cfxw))then
                deallocate(this%cfxw,this%cnbw,this%pnw)
            end if
            allocate(this%cfxw(n), this%cnbw(n))
            allocate(this%pnw(n))
            this%cfxw=0; this%pnw=0; this%cnbw=0
        end if
        this%flaglog=flag_log
        this%type_log_size=2
        this%type_int_size=this%nbin*2+3
        this%type_real_size=this%nbin*10+4
    end subroutine
    
	subroutine set_fc_xb(fc)
		implicit none
        real(8) xbg, xed
        class(sts_fc_type)::fc
        select case (fc%flaglog)
        case (fc_spacing_log)
            xbg=log10(fc%xmin); xed=log10(fc%xmax)
        case (fc_spacing_linear)
            xbg=fc%xmin;   xed= fc%xmax
        end select
        call set_range(fc%xb, fc%nbin, xbg,xed, 0)
        select case (fc%flaglog)
        case (fc_spacing_log)
            fc%xb=10**fc%xb
        case (fc_spacing_linear)
        end select
        fc%xstep=(xed-xbg)/dble(fc%nbin)		
	end subroutine
    subroutine set_range_fc(this)
        implicit none
        class(sts_fc_type)::this
        real(8) xbg, xed

        select case (this%flaglog)
        case (fc_spacing_log)
            xbg=log10(this%xmin); xed=log10(this%xmax)
        case (fc_spacing_linear)
            xbg=this%xmin;   xed= this%xmax
        end select

		! set the range of the bins
        call set_range(this%xb, this%nbin, xbg,xed, 0)

		! set the logscale or linear scale for bins
        select case (this%flaglog)
        case (fc_spacing_log)
            this%xb=10**this%xb
        case (fc_spacing_linear)
        end select

    end subroutine

    subroutine get_fc(x, n, fc)
        implicit none
		! Get the statistics of data x(n) and save it to fc according to the weight array w(n)
        integer n,i, indx
        real(8) x(n), w(n)
        real(8),allocatable:: xtp(:)
        real(8) xbg, xed
        type(sts_fc_type)::fc
        logical original_use_weight

        allocate(xtp(n))
        select case (fc%flaglog)
        case (fc_spacing_log)
            xbg=log10(fc%xmin); xed=log10(fc%xmax)
            xtp=log10(x)
        case (fc_spacing_linear)
            xbg=fc%xmin;   xed= fc%xmax
            xtp=x
        end select

		! set the range of the bins
        call set_range(fc%xb, fc%nbin, xbg,xed, 0)

		! set the logscale or linear scale for bins
        select case (fc%flaglog)
        case (fc_spacing_log)
            fc%xb=10**fc%xb
        case (fc_spacing_linear)
        end select

		! set the step size of bins
        fc%xstep=(xed-xbg)/dble(fc%nbin)

		! get nbw, nb, ns and nsw
        !original_use_weight=fc%use_weight
        !fc%use_weight=.true.
        !print*, xtp
        !print*, xbg, xed
        call get_dstr_num_in_each_bin(xtp(1:n),n,xbg,fc%xstep, fc%nbin, &
			fc%nb, fc%ns)
        !fc%fx=dble(fc%nb)/fc%xstep            
        !    do i=1, fc%nbin
        !        indx=int((xtp(i)-xbg)/fc%xstep+1)
        !        if(indx>0.and.indx<=fc%nbin)then
        !            fc%nb(indx)=fc%nb(indx)+1
        !            fc%ns=fc%ns+1
        !        end if
        !    end do
        !else
        !    call fc%get_hst(xtp,n)
        !!fc%use_weight=original_use_weight
        !end if
        
	    do i=1,fc%nbin
		    fc%fx(i)=fc%nb(i)/dble(fc%ns)/dble(fc%xstep)
           ! fc%fxw(i)=fc%nbw(i)/fc%ns/dble(fc%xstep)
            if(fc%fx(i).ne.0)then
                fc%pn(i)=fc%fx(i)/sqrt(dble(fc%nb(i)))
            !    fc%pnw(i)=fc%fxw(i)/sqrt(dble(fc%nb(i)))
            else
                fc%pn(i)=0d0
			!	fc%pnw(i)=0d0
            end if
	    end do
	    fc%cnb(1)=fc%nb(1); !fc%cnbw(1)=fc%nbw(1)
	    fc%cfx(1)=fc%fx(1)*fc%xstep
	    !fc%cfxw(1)=fc%fxw(1)*fc%xstep
	    do i=1,fc%nbin-1
		    fc%cnb(i+1)=fc%cnb(i)+fc%nb(i+1)
            !fc%cnbw(i+1)=fc%cnbw(i)+fc%nbw(i+1)
            fc%cfx(i+1) =fc%cfx(i)+fc%fx(i+1)*fc%xstep
            !fc%cfxw(i+1)=fc%cfxw(i)+fc%fxw(i+1)*fc%xstep
	    end do
    end subroutine

    subroutine get_fc_weight(x,w, n, fc)
        implicit none
		! Get the statistics of data x(n) and save it to fc according to the weight array w(n)
        integer n,i, indx
        real(8) x(n), w(n)
        real(8),allocatable:: xtp(:)
        real(8) xbg, xed
        type(sts_fc_type)::fc
        logical original_use_weight

        allocate(xtp(n))
        select case (fc%flaglog)
        case (fc_spacing_log)
            xbg=log10(fc%xmin); xed=log10(fc%xmax)
            xtp=log10(x)
        case (fc_spacing_linear)
            xbg=fc%xmin;   xed= fc%xmax
            xtp=x
        end select
        if(.not.fc%use_weight)then
            print*, "error! fc%use_weight=FALSE"
            stop
        end if
		! set the range of the bins
        call set_range(fc%xb, fc%nbin, xbg,xed, 0)

		! set the logscale or linear scale for bins
        select case (fc%flaglog)
        case (fc_spacing_log)
            fc%xb=10**fc%xb
        case (fc_spacing_linear)
        end select

		! set the step size of bins
        fc%xstep=(xed-xbg)/dble(fc%nbin)

		! get nbw, nb, ns and nsw
        !if(fc%flaglog.eq.fc_spacing_log)then
        !    do i=1, fc%nbin
        !        indx=int((xtp(i)-xbg)/fc%xstep+1)
        !        if(indx>0.and.indx<=fc%nbin)then
        !            fc%nbw(indx)=fc%nbw(indx)+w(i)
        !            fc%nb(indx)=fc%nb(indx)+1
        !            fc%nsw=fc%nsam+w(i)
        !            fc%ns=fc%ns+1
        !        end if
        !    end do
        !else
        !if(fc%use_weight=.false.)then
        !    print*, "error! fc%use_weight=false"
        !    stop
        !endif
        !original_use_weight=fc%use_weight
        !fc%use_weight=.true.
        !call fc%get_hst(xtp,w,n)
        !fc%use_weight=original_use_weight
    
        call get_dstr_num_in_each_bin(xtp(1:n),n,xbg,fc%xstep, fc%nbin, &
        fc%nb, fc%ns)
        !fc%fx=dble(fc%nb)/fc%xstep  
        call get_dstr_num_in_each_bin_weight(xtp(1:n),w(1:n),n,xbg,fc%xstep, fc%nbin, &
        fc%nbw, fc%nsw)
        !fc%fxw=fc%nbw/fc%xstep

        do i=1,fc%nbin
            fc%fx(i)=fc%nb(i)/dble(fc%ns)/dble(fc%xstep)
            fc%fxw(i)=fc%nbw(i)/fc%nsw/dble(fc%xstep)
            if(fc%fx(i).ne.0)then
                fc%pn(i)=fc%fx(i)/sqrt(dble(fc%nb(i)))
                fc%pnw(i)=fc%fxw(i)/sqrt(dble(fc%nb(i)))
            else
                fc%pn(i)=0d0
                fc%pnw(i)=0d0
            end if
        end do
        fc%cnb(1)=fc%nb(1); fc%cnbw(1)=fc%nbw(1)
        fc%cfx(1)=fc%fx(1)*fc%xstep
        fc%cfxw(1)=fc%fxw(1)*fc%xstep
        do i=1,fc%nbin-1
            fc%cnb(i+1)=fc%cnb(i)+fc%nb(i+1)
            fc%cnbw(i+1)=fc%cnbw(i)+fc%nbw(i+1)
            fc%cfx(i+1) =fc%cfx(i)+fc%fx(i+1)*fc%xstep
            fc%cfxw(i+1)=fc%cfxw(i)+fc%fxw(i+1)*fc%xstep
        end do
        !end if
    end subroutine


    subroutine fc_normalize(this)
        implicit none
        class(sts_fc_type)::this
        ! this normalize the fc even if it is not initially normalized
        real(8) wtot
        integer i

        wtot=0
        do i=1, this%nbin
            wtot=wtot+this%xstep*this%fx(i)
        end do
        this%fx=this%fx/wtot
    end subroutine
    subroutine output_sts_data_weight(x,w,n, xmin,xmax,rn,flaglog, fn)
        implicit none
		!output the statistics of data into file
        type(sts_fc_type)::fc
        character*(*) fn
        integer flaglog, rn, n
        real(8) x(n), w(n), xmin,xmax
        call init_fc(fc, xmin, xmax, rn, flaglog)
        call get_fc_weight(x(1:n), w(1:n) , n, fc)
        call output_fc_weight(fc, trim(adjustl(fn)))
        call output_nfc_weight(fc, trim(adjustl(fn)))
    end subroutine

    subroutine output_sts_data_weight_auto(x,w,n,rn,flaglog, fn)
        implicit none
        type(sts_fc_type)::fc
        character*(*) fn
        integer flaglog, rn, n
        real(8) x(n), w(n), xmin,xmax
		xmin=minval(x(1:n)); xmax=maxval(x(1:n))
        call init_fc(fc, xmin, xmax, rn, flaglog)
        call get_fc_weight(x(1:n), w(1:n) , n, fc)
        call output_fc_weight(fc, trim(adjustl(fn)))
        call output_nfc_weight(fc, trim(adjustl(fn)))
    end subroutine

    subroutine output_sts_data(x,n, xmin,xmax,rn,flaglog, fn)
        implicit none
        type(sts_fc_type)::fc
        character*(*) fn
        integer flaglog, rn, n
        real(8) x(n), xmin,xmax
        call init_fc(fc, xmin, xmax, rn, flaglog)
        call get_fc(x(1:n) , n, fc)
        call output_fc(fc, trim(adjustl(fn)))
        call output_nfc(fc, trim(adjustl(fn)))
    end subroutine

    subroutine output_sts_data_auto(x,n,rn,flaglog, fn)
        implicit none
        type(sts_fc_type)::fc
        character*(*) fn
        integer flaglog, rn, n
        real(8) x(n), w(n), xmin,xmax

		xmin=minval(x(1:n)); xmax=maxval(x(1:n))
        call init_fc(fc, xmin, xmax, rn, flaglog)
        call get_fc(x(1:n), n, fc)
        call output_fc(fc, trim(adjustl(fn)))
        call output_nfc(fc, trim(adjustl(fn)))
    end subroutine

    subroutine get_sts_data_weight(x,w,n, xmin,xmax,rn,flaglog, fc)
        implicit none
        type(sts_fc_type)::fc
        integer flaglog, rn, n
        real(8) x(n), w(n), xmin,xmax
        call init_fc(fc, xmin, xmax, rn, flaglog,use_weight=.true.)
        call get_fc_weight(x(1:n), w(1:n) , n, fc)
    end subroutine

    subroutine get_sts_data_weight_auto(x,w,n, rn,flaglog, fc)
        implicit none
        type(sts_fc_type)::fc
        integer flaglog, rn, n
        real(8) x(n), w(n), xmin,xmax
        xmin=minval(x(1:n)); xmax=maxval(x(1:n))
        call init_fc(fc, xmin, xmax, rn, flaglog,use_weight=.true.)
        call get_fc_weight(x(1:n), w(1:n) , n, fc)
    end subroutine

    subroutine get_sts_data_auto(x,n,rn,flaglog, fc)
        implicit none
        type(sts_fc_type)::fc
        integer flaglog, rn, n
        real(8) x(n), xmin,xmax
		xmin=minval(x(1:n)); xmax=maxval(x(1:n))
        call init_fc(fc, xmin, xmax, rn, flaglog)
        call get_fc(x(1:n) , n, fc)
    end subroutine

    subroutine get_sts_data(x,n, xmin,xmax,rn,flaglog, fc)
        implicit none
        type(sts_fc_type)::fc
        integer flaglog, rn, n
        real(8) x(n), xmin,xmax
        call init_fc(fc, xmin, xmax, rn, flaglog)
        call get_fc(x(1:n) , n, fc)
    end subroutine

    subroutine output_fc_bin(fc, fn)
        implicit none
        type(sts_fc_type)::fc
        character*(*) fn
        integer i
        open(unit=1999,file=trim(adjustl(fn))//"_fc.bin",access='stream', form='unformatted')
        write(unit=1999) fc
        close(unit=1999)
    end subroutine

    subroutine input_fc_bin(fc, fn)
        implicit none
        type(sts_fc_type)::fc
        character*(*) fn
        integer i,flaglog,rn
		open(unit=1999,file=trim(adjustl(fn))//"_fc.bin",access='stream', form='unformatted',status='old')
        read(unit=1999) fc
        close(unit=1999)

    end subroutine
    subroutine output_fc(fc, fn)
        implicit none
        class(sts_fc_type)::fc
        character*(*) fn
        integer i
        open(unit=999,file=trim(adjustl(fn))//"_fc.txt")
        write(unit=999,fmt="(10A27)") "xb", "fx", "pn" , "cfx" 
        do i=1, fc%nbin
			!if(isnan(fc%fxw(i)).or.fc%ns.eq.0)then
			!	print*, "warnning: i, fx, fxw, ns=", i,fc%fx(i), fc%fxw(i), fc%ns
			!end if
            write(unit=999,fmt="(1P10E27.12E4)") fc%xb(i), fc%fx(i), fc%pn(i) , fc%cfx(i)
        end do
    end subroutine

    subroutine output_fc_weight(fc, fn)
        implicit none
        class(sts_fc_type)::fc
        character*(*) fn
        integer i
        open(unit=999,file=trim(adjustl(fn))//"_fc.txt")
        write(unit=999,fmt="(10A27)") "xb", "fx", "pn", "fxw", "pnw", "cfx", "cfxw"
        do i=1, fc%nbin
			!if(isnan(fc%fxw(i)).or.fc%ns.eq.0)then
			!	print*, "warnning: i, fx, fxw, ns=", i,fc%fx(i), fc%fxw(i), fc%ns
			!end if
            write(unit=999,fmt="(1P10E27.12E4)") fc%xb(i), fc%fx(i), fc%pn(i), fc%fxw(i), fc%pnw(i), fc%cfx(i), fc%cfxw(i)
        end do
    end subroutine

    subroutine input_fc(fc, fn,xmin,xmax,rn,flaglog)
        implicit none
        type(sts_fc_type)::fc
        character*(*) fn
        integer i,flaglog,rn
		real(8) xmin,xmax

        open(unit=999,file=trim(adjustl(fn))//"_fc.txt")
        read(unit=999,fmt=*) 
		call init_fc(fc, xmin, xmax, rn, flaglog)
        do i=1, fc%nbin
            read(unit=999,fmt="(1P10E27.12E4)") fc%xb(i), fc%fx(i), fc%pn(i), fc%fxw(i), fc%pnw(i), fc%cfx(i), fc%cfxw(i)
			if(isnan(fc%fxw(i)))then
				print*, "warnning: i, fx, fxw=", i,fc%fx(i), fc%fxw(i), fc%ns
			end if
        end do
    end subroutine

    subroutine output_nfc(fc, fn)
        implicit none
        type(sts_fc_type)::fc
        character*(*) fn
        integer i
        open(unit=999,file=trim(adjustl(fn))//"_nfc.txt")
        write(unit=999,fmt="(1A27,10A10)") "xb",  "nb", "cnb"
        do i=1, fc%nbin
            write(unit=999,fmt="(1P1E27.12E4, 10I10)") fc%xb(i),   fc%nb(i), fc%cnb(i)
        end do
    end subroutine
    subroutine output_nfc_weight(fc, fn)
        implicit none
        type(sts_fc_type)::fc
        character*(*) fn
        integer i
        open(unit=999,file=trim(adjustl(fn))//"_nfc.txt")
        write(unit=999,fmt="(3A27,10A10)") "xb", "nbw", "cnbw", "nb", "cnb"
        do i=1, fc%nbin
            write(unit=999,fmt="(1P3E27.12E4, 10I10)") fc%xb(i), fc%nbw(i), fc%cnbw(i), fc%nb(i), fc%cnb(i)
        end do
    end subroutine
	subroutine get_percentage_position(fc, percent, pos, n, method)
		implicit none
		type(sts_fc_type)::fc
		integer i, n
		real(8) percent(n), pos(4, n)
		integer method,idx
		integer,parameter::method_linear=1	
		do i=1, n
			call search_for_position(fc%cfx, fc%xb, fc%nbin, percent(i), pos(1, i))
			call search_for_position(fc%cfxw, fc%xb, fc%nbin, percent(i), pos(2, i))
		end do
	end subroutine
	subroutine search_for_position(y, x, n, per, pos)
		implicit none
		integer i, n
		real(8) x(n), y(n), per, pos
		do i=1, n-1
			if(per.eq.y(i))then
				pos=x(i)
				return
			else
				if(per>y(i).and.per<y(i+1))then
					pos=(x(i+1)-x(i))/(y(i+1)-y(i))*(per-y(i))+x(i)
				end if
				return
			end if
		end do
	end subroutine
    subroutine get_value_fc_linear_y(fc,va,xout)
		implicit none
		class(sts_fc_type)::fc
		real(8) va,xout
		integer n
		
		n=fc%nbin
        call linear_int_arb(fc%fx(1:n),fc%xb(1:n),fc%nbin, va,xout)
        
	end subroutine

    subroutine get_value_fc_spline_y(fc,va,yout)
		implicit none
		class(sts_fc_type)::fc
		real(8) va,yout
		integer n,ier

        if(fc%fx(1)>va)then
			yout=fc%xb(1)
		else if(fc%fx(fc%nbin)<va)then
			yout=fc%xb(fc%nbin)
		else
			n=fc%nbin
			call splint_mylib(fc%fx(1:n),fc%xb(1:n),fc%y2_y(1:n),fc%nbin, va, yout,ier)
            !print*, "fc%y2_y=",fc%y2_y
            !print*, "xb=",fc%xb
            !print*, "fc=",fc%fx
			!print*, "va,yout=",va, yout
		end if
		
	end subroutine
	!subroutine get_value_fc_direct(fc,va,yout)
	!	implicit none
	!	class(sts_fc_type)::fc
	!	real(8) va,yout
	!	integer idx
	!	
!
	!	if(fc%xb(1)>va)then
	!		yout=fc%fx(1)
	!	else if(fc%xb(fc%nbin)<va)then
	!		yout=fc%fx(fc%nbin)
	!	else
	!		call return_idx(va,fc%xmin,fc%xmax, fc%nbin, idx,0)
	!	end if
	!end subroutine
	
	subroutine get_value_at_x_fc(fc, va, y_out, method)
		implicit none
		type(sts_fc_type)::fc
		real(8) va, y_out(6)
		integer method,idx,n,ier
        integer,parameter::method_linear=1
        integer,parameter::method_intp=2
        integer,parameter::method_direct=3
        !print*, fc%xmax<va.or.fc%xmin>va
        if(fc%flaglog.eq.fc_spacing_log)then
            print*, "warnning!, fc spacing log could cause inaccuracy"
        end if
		if(fc%xb(1)>va.or.fc%xb(fc%nbin)<va)then
		!if(fc%xmin>va.or.fc%xmax<va)then
			y_out(1:2)=0d0
			!print*, "????"
			!stop
			!if(fc%xmax<va)then
			if(fc%xb(fc%nbin)<va)then
				y_out(3:4)=1d0
			end if
		else
			select case(method)
			case(method_intp)
				call splint_mylib(fc%xb,fc%fx,fc%y2,fc%nbin, va, y_out(1),ier)
			case(method_linear)
                n=fc%nbin
				call linear_int(fc%xb(1:n),fc%fx(1:n),fc%nbin, va,y_out(1))
				call linear_int(fc%xb(1:n),fc%fxw(1:n),fc%nbin, va,y_out(2))
				call linear_int(fc%xb(1:n),fc%cfx(1:n),fc%nbin, va,y_out(3))
				call linear_int(fc%xb(1:n),fc%cfxw(1:n),fc%nbin, va,y_out(4))
				call linear_int(fc%xb(1:n),dble(fc%cnb(1:n)),fc%nbin, va,y_out(5))
				call linear_int(fc%xb(1:n),fc%cnbw(1:n),fc%nbin, va,y_out(6))
			case(method_direct)
				call return_idx(va,fc%xmin,fc%xmax, fc%nbin, idx,0)
                !call fc%print()
				!print*, "va, idx=",va, idx
                !read(*,*)
                if(idx>0.and.idx<=fc%nbin)then
				    y_out(1)=fc%fx(idx)
                else
                    y_out(1)=0
                end if
			!	call return_idx(va,fc%xmin,fc%xmax, fc%nbin, idx,fc%flaglog)
				!print*, "idx=",idx
			!	y_out=fc%fx(idx)
			case default
				print*, "error, defined the method", method
				stop
			end select
		end if
	end subroutine

    subroutine init_intrp_fc_y(fc)
        implicit none
        class(sts_fc_type)::fc
        integer n
        real(8) yp1,ypn
        yp1=1d30;ypn=1d30
        n=fc%nbin
        call spline_mylib(fc%fx(1:n),fc%xb(1:n),fc%nbin,yp1,ypn,fc%y2_y(1:n))
    end subroutine
    subroutine sort_fc(this,ind_col)
        !use stdlib_sorting 
        implicit none
        class(sts_fc_type)::this
        integer ind_col
        integer,parameter:: ind_xb=0,ind_fx=1,ind_fxw=2
        integer indx(this%nbin)

        select case(ind_col)
        case(ind_xb)
            call indexx(this%nbin,this%xb, indx)
        case(ind_fx)
            call indexx(this%nbin,this%fx, indx)
        case(ind_fxw)
            call indexx(this%nbin,this%fxw, indx)
        end select
       ! print*, "indx=",indx
        call sort_arr_index(this%nbin,this%xb,indx)
        call sort_arr_index(this%nbin,this%fx,indx)
        call sort_arr_index(this%nbin,this%fxw,indx)
        
    end subroutine


	subroutine conv_fc_int_real_arrays(fc, intarr, realarr, logarr)
		!use com_main_gw
		implicit none
		type(sts_fc_type)::fc
		integer nint, nreal
        logical logarr(fc%type_log_size)
		integer intarr(fc%type_int_size)
		real(8) realarr(fc%type_real_size)

        logarr=(/fc%use_weight, fc%is_spline_prepared/)

		intarr(1:3)=(/fc%flaglog,fc%nbin, fc%ns/)
		intarr(4:3+fc%nbin)=fc%nb(1:fc%nbin)
		intarr(4+fc%nbin:3+fc%nbin*2)=fc%cnb(1:fc%nbin)

		realarr(1:fc%nbin)=fc%xb(1:fc%nbin)
		realarr(fc%nbin+1:2*fc%nbin)=fc%fx(1:fc%nbin)
		realarr(fc%nbin*5+1:6*fc%nbin)=fc%cfx(1:fc%nbin)
		realarr(fc%nbin*7+1:8*fc%nbin)=fc%pn(1:fc%nbin)
		if(fc%is_spline_prepared)then
		    realarr(fc%nbin*9+1:10*fc%nbin)=fc%y2(1:fc%nbin)
        end if
        if(fc%use_weight)then
            realarr(fc%nbin*2+1:3*fc%nbin)=fc%fxw(1:fc%nbin)
            realarr(fc%nbin*3+1:4*fc%nbin)=fc%nbw(1:fc%nbin)
            realarr(fc%nbin*4+1:5*fc%nbin)=fc%cnbw(1:fc%nbin)
            realarr(fc%nbin*8+1:9*fc%nbin)=fc%pnw(1:fc%nbin)
            realarr(fc%nbin*6+1:7*fc%nbin)=fc%cfxw(1:fc%nbin)
		end if        
		realarr(fc%nbin*10+1:10*fc%nbin+4)=&
				(/fc%xmin,fc%xmax,fc%nsw,fc%xstep/)

	end subroutine
	subroutine conv_int_real_arrays_fc(fc, intarr, realarr,logarr)
		!use com_main_gw
		implicit none
		type(sts_fc_type)::fc
		integer nint, nreal
        logical logarr(fc%type_log_size)
		integer intarr(fc%type_int_size)
		real(8) realarr(fc%type_real_size)

        fc%use_weight=logarr(1)
        fc%is_spline_prepared=logarr(2)
		fc%flaglog=intarr(1); fc%nbin=intarr(2)
        fc%ns=intarr(3)
		fc%nb(1:fc%nbin)=intarr(4:3+fc%nbin)
		fc%cnb(1:fc%nbin)=intarr(4+fc%nbin:3+fc%nbin*2)

		fc%xb(1:fc%nbin)=realarr(1:fc%nbin)
		fc%fx(1:fc%nbin)=realarr(1+fc%nbin:fc%nbin*2)
		fc%cfx(1:fc%nbin)=realarr(1+fc%nbin*5:fc%nbin*6)
		fc%pn(1:fc%nbin)=realarr(1+fc%nbin*7:fc%nbin*8)
		if(fc%is_spline_prepared)then
            if(.not.allocated(fc%y2))allocate(fc%y2(fc%nbin))
            fc%y2(1:fc%nbin)=realarr(1+fc%nbin*9:fc%nbin*10)
        end if
        if(fc%use_weight)then
            fc%fxw(1:fc%nbin)=realarr(1+fc%nbin*2:fc%nbin*3)
            fc%nbw(1:fc%nbin)=realarr(1+fc%nbin*3:fc%nbin*4)
            fc%cnbw(1:fc%nbin)=realarr(1+fc%nbin*4:fc%nbin*5)
            fc%pnw(1:fc%nbin)=realarr(1+fc%nbin*8:fc%nbin*9)
            fc%cfxw(1:fc%nbin)=realarr(1+fc%nbin*6:fc%nbin*7)
        end if
		fc%xmin=realarr(fc%nbin*10+1)
		fc%xmax=realarr(fc%nbin*10+2)
		fc%nsw=realarr(fc%nbin*10+3)
		fc%xstep=realarr(fc%nbin*10+4)
		!fc%ns=realarr(fc%nbin*10+5)
	end subroutine

#ifdef USE_HDF5	
	subroutine save_fc_hdf5(fc, group_id,  tablename)
		use md_hdf5_table
		implicit none
		class(sts_fc_type)::fc
		character*(*) tablename
		INTEGER(HID_T) :: group_id      ! Group identifier
        !INTEGER(HSIZE_T), PARAMETER :: nfields  = 11            ! nfields
        type(hdf5_table_type)::htable
		if(fc%nbin.eq.0) return
        

        if(fc%use_weight)then
            !print*, "fc%nbin,tablename=",fc%nbin,tablename,allocated(fc%xb)
            call htable%init_table(11, fc%nbin,  tablename)
            htable%field_names=(/"   X","  FX"," FXW"," cfx", "cfxw", &
                        "  PN"," PNW"," nbw","cnbw","  nb", " cnb"/)	
            !print*, "ini finished"
            htable%field_types(1:9)=H5T_NATIVE_DOUBLE
            htable%field_types(10:11)=H5T_NATIVE_INTEGER

            call htable%prepare_write_table(group_id)
            !print*, "prepare finished"
            CALL htable%write_column_real(fc%xb)
            CALL htable%write_column_real(fc%fx)
            CALL htable%write_column_real(fc%fxw)
            CALL htable%write_column_real(fc%cfx)
            CALL htable%write_column_real(fc%cfxw)
            CALL htable%write_column_real(fc%pn)
            CALL htable%write_column_real(fc%pnw)
            CALL htable%write_column_real(fc%nbw)
            CALL htable%write_column_real(fc%cnbw)
            CALL htable%write_column_int(fc%nb)
            CALL htable%write_column_int(fc%cnb)
        else
            call htable%init_table(6, fc%nbin,  tablename)
            htable%field_names=(/"   X","  FX", " cfx", "  PN","  nb", " cnb"/)	

            htable%field_types(1:4)=H5T_NATIVE_DOUBLE
            htable%field_types(5:6)=H5T_NATIVE_INTEGER

            call htable%prepare_write_table(group_id)

            CALL htable%write_column_real(fc%xb)
            CALL htable%write_column_real(fc%fx)
            CALL htable%write_column_real(fc%cfx)
            CALL htable%write_column_real(fc%pn)
            CALL htable%write_column_int(fc%nb)
            CALL htable%write_column_int(fc%cnb)
		endif
	end subroutine
	subroutine read_fc_hdf5(fc, group_id, tablename,use_weight)
		use md_hdf5_table
		implicit none
		class(sts_fc_type)::fc
		character*(*) tablename
		INTEGER(HID_T) :: group_id      ! Group identifier
        !INTEGER(HSIZE_T), PARAMETER :: nfields  = 11            ! nfields
        type(hdf5_table_type)::htable
        logical use_weight
		if(fc%nbin.eq.0) return
        if(use_weight)then
            call htable%init_table(11, fc%nbin,  tablename)
            
            !nrecords=fc%nbin;
            htable%field_names=(/"   X","  FX"," FXW"," cfx", "cfxw", & 
                        "  PN"," PNW"," nbw","cnbw", "  nb", " cnb"/)	

            htable%field_types(1:9)=H5T_NATIVE_DOUBLE
            htable%field_types(10:11)=H5T_NATIVE_INTEGER

            call htable%prepare_read_table(group_id)

            CALL htable%read_column_real(fc%xb)
            CALL htable%read_column_real(fc%fx)
            CALL htable%read_column_real(fc%fxw)
            CALL htable%read_column_real(fc%cfx)
            CALL htable%read_column_real(fc%cfxw)
            CALL htable%read_column_real(fc%pn)
            CALL htable%read_column_real(fc%pnw)
            CALL htable%read_column_real(fc%nbw)
            CALL htable%read_column_real(fc%cnbw)
            CALL htable%read_column_int(fc%nb)
            CALL htable%read_column_int(fc%cnb)
        else
            call htable%init_table(6, fc%nbin,  tablename)
            
            !nrecords=fc%nbin;
            htable%field_names=(/"   X","  FX", " cfx", "  PN","  nb", " cnb"/)		

            htable%field_types(1:9)=H5T_NATIVE_DOUBLE
            htable%field_types(10:11)=H5T_NATIVE_INTEGER

            call htable%prepare_read_table(group_id)

            CALL htable%read_column_real(fc%xb)
            CALL htable%read_column_real(fc%fx)
            CALL htable%read_column_real(fc%cfx)
            CALL htable%read_column_real(fc%pn)
            CALL htable%read_column_int(fc%nb)
            CALL htable%read_column_int(fc%cnb)
        end if
		
	end subroutine
#endif

end module
