
module md_s1d_hst_basic_type
	use my_intgl
	!use md_sts_fc
	use md_s1d_basic_type
	integer,parameter:: f_log=0, f_linear=1
	integer,parameter::dct_x=1,dct_y=2
	
	type,extends(s1d_basic_type):: s1d_hst_basic_type
		integer, allocatable:: nb(:)
		real(8), allocatable:: fxw(:), nbw(:)
		real(8) nsw
        integer ns
		logical use_weight
        ! xb: the central position of each bin
        ! nb: unweighted number of samples per bin
		!nbw:   weighted number of samples per bin
		! fx:  unweighted distribution function of samples per pin
        ! fxw:   weighted distribution function of samples per pin
        ! nsw:      weighted total number of samples in all bins
        ! ns:      unweighted total number of samples in all bins

        contains
        procedure::print=>print_s1_hst
        procedure::get_s1d_hst
        procedure::get_s1d_hst_weight        
        generic::get_hst=>get_s1d_hst, get_s1d_hst_weight
	end type
	!interface get_hst
	!	module procedure get_hst_no_weight
	!	module procedure get_hst_weight
	!end interface
    private::print_s1_hst,get_s1d_hst,get_s1d_hst_weight
    
contains
	subroutine init_s1_hst_basic(this, xmin,xmax, n, use_weight)
		implicit none
		class(s1d_hst_basic_type)::this
		integer n
		real(8) xmin,xmax
		logical,optional:: use_weight
		if(n==0)then
			print*, "init_s1d:warnning s1d%nbin=0"
		end if
        !print*, "1"
		call init_s1d_basic(this, xmin,xmax,n, sts_type_dstr)
        !print*, "2"
		if(present(use_weight))then
			this%use_weight=use_weight
		else
			this%use_weight=.false.
		end if
        !print*, "3"
		if(allocated(this%nb))then
			deallocate(this%nb)
			if(this%use_weight)then
				deallocate(this%nbw, this%fxw)
			end if
		end if
		!print*, "4"
		allocate(this%nb(n))
		this%nb=0;this%ns=0
		if(this%use_weight)then
			allocate(this%nbw(n), this%fxw(n))
			this%nbw=0; this%fxw=0; this%nsw=0
		end if
		
	end subroutine

    subroutine print_s1_hst(this,str_,print_y2_in)
        implicit none
        class(s1d_hst_basic_type)::this
        character*(*),optional ::str_
		logical,optional::print_y2_in
		logical::print_y2
        character*(15) str_bin_type		
        integer i
        !print*, "?",str_,trim(adjustl(str_))
        !print*, present(str_)
		if(present(str_))then
        !print*, "0"
            write(*,*) "for sts s1d=", trim(adjustl(str_))
        end if
        !print*, "???"
		select case (this%bin_type)
        case (sts_type_grid)
            str_bin_type="GRID"
        case (sts_type_dstr)
            str_bin_type="DSTR"
        end select
        !print*, "1"
        
        !print*, "3"
		if(this%use_weight)then
			write(*, fmt="(10A15)") "bin_type", "nbin", "xmin", "xmax", "xstep", "NS", "NSW"
			write(*, fmt="(A15, I15, 3E15.5, I15, E15.5 )") str_bin_type, this%nbin, &
            this%xmin, this%xmax, this%xstep, this%ns, this%nsw
			write(unit=*, fmt="(20A20)") "X", "FX", "FXW", "NB", "NBW"
			do i=1, this%nbin
				write(unit=*, fmt="(3E20.10, I20, E20.10)") this%xb(i), this%fx(i), this%fxw(i), this%nb(i), this%nbw(i)
			end do
		else
			write(*, fmt="(10A15)") "bin_type", "nbin", "xmin", "xmax", "xstep", "NS"
			!print*, "2"
			write(*, fmt="(A15, I15, 3E15.5, I15 )") str_bin_type, this%nbin, &
            this%xmin, this%xmax, this%xstep, this%ns
			write(unit=*, fmt="(20A20)") "X", "FX", "NB"
			do i=1, this%nbin
				write(unit=*, fmt="(2E20.10, I20)") this%xb(i), this%fx(i), this%nb(i)
			end do
		end if
        !print*, "4"
    end subroutine
	
!	subroutine get_value_s1d(s1d, va, y_out, method)
!		implicit none
!		type(s1d_hst_basic_type)::s1d
!		real(8) va, y_out
!		integer method,idx
!		if(s1d%xmax<va.or.s1d%xmin>va)then
!			y_out=0d0
!		else
!			select case(method)
!			case(method_intp_s1d)
!				call splint_mylib(s1d%xb,s1d%fx,s1d%y2,s1d%nbin, va, y_out)
!			case(method_linear_s1d)
!				call linear_int(s1d%xb,s1d%fx,s1d%nbin, va,y_out)
!			case(method_direct_s1d)
!				call return_idx(va,s1d%xmin,s1d%xmax, s1d%nbin, idx,s1d%bin_type)
!				!print*, "idx=",idx
!				y_out=s1d%fx(idx)
!			end select
!		end if
!	end subroutine

	

	!subroutine intp_int_s1d(s1d, sums)
	!	implicit none
	!	type(s1d_hst_basic_type)::s1d
	!	integer i,j
!	!	integer,optional,intent(in):: use_ylog
	!	real(8) sums
	!	sums=0
	!	call my_integral_interpolate(s1d%xmin,s1d%xmax, s1d%xb,s1d%fx,s1d%nbin, sums)
	!	
	!end subroutine
	!subroutine intp_int_s1d_acc(s1d, atol, rtol, sums)
	!	implicit none
	!	type(s1d_hst_basic_type)::s1d
	!	integer i,j
	!	real(8) atol(1),rtol(1)
!	!	integer,optional,intent(in):: use_ylog
	!	real(8) sums
	!	sums=0
	!	call my_integral_interpolate_acc(s1d%xmin,s1d%xmax, s1d%xb,s1d%fx,s1d%nbin, atol, rtol, sums)
	!	
	!end subroutine
	!subroutine intp_int_s1d_ylog(s1d, sums)
	!	implicit none
	!	type(s1d_hst_basic_type)::s1d
	!	integer i,j
!	!	integer,optional,intent(in):: use_ylog
	!	real(8) sums
	!	sums=0
	!	call my_integral_interpolate_ylog(s1d%xmin,s1d%xmax, s1d%xb,log10(s1d%fx),s1d%nbin, sums)
	!	
	!end subroutine
	!subroutine intp_int_s1d_par(s1d, ylog, atol, rtol, sums)
	!	implicit none
	!	type(s1d_hst_basic_type)::s1d
	!	integer i,j
	!	real(8) atol(1),rtol(1)
	!	logical ylog
	!	real(8) sums
	!	sums=0
	!	if(ylog)then
	!		call my_integral_interpolate_par(s1d%xmin,s1d%xmax, s1d%xb,log10(s1d%fx),s1d%nbin, atol, rtol,ylog,sums)
	!		
	!	else
	!		call my_integral_interpolate_par(s1d%xmin,s1d%xmax, s1d%xb,s1d%fx,s1d%nbin, atol,rtol,ylog,sums)
	!		
	!	end if
	!end subroutine
	
	subroutine get_s1d_hst_weight(s1_hst,x,w,n)
		implicit none
		class(s1d_hst_basic_type)::s1_hst
		!type(sts_fc_type)::s1d
		integer n
		real(8) x(n),w(n)
		integer i, indx
		!call init_fc(s1d,s1d%xmin,s1d%xmax,s1d%nbin,fc_spacing_linear)
		if(.not.s1_hst%use_weight)then
			print*, "error! s1_hst%use_weight=FALSE"
			stop
		end if
		call get_dstr_num_in_each_bin_weight(x,w,n,s1_hst%xmin,s1_hst%xstep, s1_hst%nbin, &
			s1_hst%nbw, s1_hst%nsw)
		s1_hst%fxw(1:s1_hst%nbin)=s1_hst%nbw(1:s1_hst%nbin)/s1_hst%xstep

		call get_dstr_num_in_each_bin(x(1:n),n,s1_hst%xmin,s1_hst%xstep, s1_hst%nbin, &
			s1_hst%nb, s1_hst%ns)
		s1_hst%fx(1:s1_hst%nbin)=dble(s1_hst%nb(1:s1_hst%nbin))/s1_hst%xstep
	end subroutine
	subroutine get_s1d_hst(s1_hst,x,n)
		implicit none
		class(s1d_hst_basic_type)::s1_hst
		!type(sts_fc_type)::s1d
		integer n
		real(8) x(n)!,w(n)
		integer i, indx
		!call init_fc(s1d,s1d%xmin,s1d%xmax,s1d%nbin,fc_spacing_linear)
		if(s1_hst%use_weight)then
			print*, "warnning: s1_hst%use_weight=True, but there is no input weighting data "
		endif
		!w=1d0
        !print*, "start"

		call get_dstr_num_in_each_bin(x(1:n),n,s1_hst%xmin,s1_hst%xstep, s1_hst%nbin, &
			s1_hst%nb, s1_hst%ns)
		s1_hst%fx(1:s1_hst%nbin)=dble(s1_hst%nb(1:s1_hst%nbin))/s1_hst%xstep
        !print*, "end"
	end subroutine	

	!subroutine output_sts_1d(s1d,fn)
	!	implicit none
	!	type(s1d_hst_basic_type)::s1d
	!	character*(*) fn
	!	integer i
	!	character*(8) tmp
!
	!	open(unit=999,file=trim(adjustl(fn)))
	!	do i=1, s1d%nbin
	!		write(unit=999,fmt="(3E20.10E4)") s1d%xb(i), s1d%fx(i), s1d%rnx(i)
	!	end do
	!	close(unit=999)
!
	!end subroutine
	subroutine output_s1d_hst(s1d_hst, fn)
        implicit none
        type(s1d_hst_basic_type)::s1d_hst
        character*(*) fn
        integer i
        open(unit=999,file=trim(adjustl(fn))//"_s1d_hst.txt")
		if(s1d_hst%use_weight)then
			write(unit=999,fmt="(10A27)") "xb", "fx", "fxw", "nb", "nbw"
			do i=1, s1d_hst%nbin
				!if(isnan(s1d%fxw(i)).or.s1d%ns.eq.0)then
				!	print*, "warnning: i, fx, fxw, ns=", i,s1d%fx(i), s1d%fxw(i), s1d%ns
				!end if
				write(unit=999,fmt="(1P3E27.12E4, I27, 1PE27.12E4)") s1d_hst%xb(i), s1d_hst%fx(i), &
					s1d_hst%fxw(i), s1d_hst%nb(i), s1d_hst%nbw(i)
			end do
		else
			write(unit=999,fmt="(10A27)") "xb", "fx", "nb"
			do i=1, s1d_hst%nbin
				write(unit=999,fmt="(1P2E27.12E4, I27)") s1d_hst%xb(i), s1d_hst%fx(i), s1d_hst%nb(i)
			end do
		end if
		close(999)
    end subroutine
	
end module

module md_s1d_hst_type
	use md_s1d_hst_basic_type
	type,extends(s1d_hst_basic_type):: s1d_hst_type
		integer type_int_size
		integer type_real_size
		integer type_log_size
	contains
		procedure::init=>init_s1d_hst

		procedure::read_s1d_hst
		procedure::write_s1d_hst
#ifdef USE_HDF5		
		procedure::save_hdf5=>save_s1d_hst_hdf5
		procedure::read_hdf5=>read_s1d_hst_hdf5
#endif
		generic :: read(unformatted) => read_s1d_hst
		generic :: write(unformatted) => write_s1d_hst
		
	end type
	private init_s1d_hst
	private::read_s1d_hst
	private::write_s1d_hst
#ifdef USE_HDF5	
	private::save_s1d_hst_hdf5, read_s1d_hst_hdf5
#endif
contains 
	subroutine init_s1d_hst(this, xmin,xmax, n, use_weight)
		implicit none
		class(s1d_hst_type)::this
		integer n
		real(8) xmin,xmax
		logical,optional:: use_weight
		if(present(use_weight))then
			call init_s1_hst_basic(this,  xmin,xmax,n, use_weight)
		else
			call init_s1_hst_basic(this,  xmin,xmax,n, .false.)
		end if
		this%type_log_size=2
        this%type_int_size=this%nbin+3
        this%type_real_size=this%nbin*6+4
	end subroutine

	subroutine read_s1d_hst(s1d,file_unit, iostat, iomsg)
		implicit none
		class(s1d_hst_type), intent(inout)::s1d
		integer,intent(in):: file_unit
		integer,intent(out)::iostat
		character(*), intent(inout) :: iomsg
		integer n
		read(unit=file_unit) s1d%nbin, s1d%xmin,s1d%xmax
		read(unit=file_unit) s1d%use_weight
		n=s1d%nbin
		if(n.eq.0)return
		call s1d%init(s1d%xmin,s1d%xmax,n,s1d%use_weight)
		read(unit=file_unit) s1d%xb(1:n), s1d%fx(1:n), s1d%nb(1:n), s1d%ns
		if(s1d%use_weight)then
			read(unit=file_unit)  s1d%nbw(1:n), s1d%fxw(1:n), s1d%nsw
		endif
	end subroutine
	subroutine write_s1d_hst(s1d,file_unit, iostat, iomsg)
		implicit none
		class(s1d_hst_type), intent(in)::s1d
		integer,intent(in):: file_unit
		integer,intent(out)::iostat
		character(*), intent(inout) :: iomsg
		integer n

		write(unit=file_unit) s1d%nbin, s1d%xmin,s1d%xmax
		write(unit=file_unit) s1d%use_weight
		n=s1d%nbin
		if(n.eq.0)return
		write(unit=file_unit) s1d%xb(1:n), s1d%fx(1:n), s1d%nb(1:n), s1d%ns
		if(s1d%use_weight)then
			write(unit=file_unit)  s1d%nbw(1:n), s1d%fxw(1:n), s1d%nsw
		end if
	end subroutine
	subroutine conv_s1d_hst_int_real_arrays(fc, intarr, realarr, logarr)
		!use com_main_gw
		implicit none
		type(s1d_hst_type)::fc
		integer nint, nreal
        logical logarr(fc%type_log_size)
		integer intarr(fc%type_int_size)
		real(8) realarr(fc%type_real_size)

        logarr=(/fc%use_weight, fc%is_spline_prepared/)

		intarr(1:2)=(/fc%nbin, fc%ns/)
		intarr(3:2+fc%nbin)=fc%nb(1:fc%nbin)

		realarr(1:fc%nbin)=fc%xb(1:fc%nbin)
		realarr(fc%nbin+1:2*fc%nbin)=fc%fx(1:fc%nbin)
		if(fc%is_spline_prepared)then
		    realarr(fc%nbin*5+1:6*fc%nbin)=fc%y2(1:fc%nbin)
        end if
        if(fc%use_weight)then
            realarr(fc%nbin*2+1:3*fc%nbin)=fc%fxw(1:fc%nbin)
            realarr(fc%nbin*3+1:4*fc%nbin)=fc%nbw(1:fc%nbin)
		end if        
		realarr(fc%nbin*6+1:6*fc%nbin+4)=&
				(/fc%xmin,fc%xmax,fc%nsw,fc%xstep/)

	end subroutine
	subroutine conv_int_real_arrays_s1d_hst(fc, intarr, realarr,logarr)
		!use com_main_gw
		implicit none
		type(s1d_hst_type)::fc
		integer nint, nreal
        logical logarr(fc%type_log_size)
		integer intarr(fc%type_int_size)
		real(8) realarr(fc%type_real_size)

        fc%use_weight=logarr(1)
        fc%is_spline_prepared=logarr(2)
		fc%nbin=intarr(1)
        fc%ns=intarr(2)
		fc%nb(1:fc%nbin)=intarr(3:2+fc%nbin)

		fc%xb(1:fc%nbin)=realarr(1:fc%nbin)
		fc%fx(1:fc%nbin)=realarr(1+fc%nbin:fc%nbin*2)
		if(fc%is_spline_prepared)then
            if(.not.allocated(fc%y2))allocate(fc%y2(fc%nbin))
            fc%y2(1:fc%nbin)=realarr(1+fc%nbin*5:fc%nbin*6)
        end if
        if(fc%use_weight)then
            fc%fxw(1:fc%nbin)=realarr(1+fc%nbin*2:fc%nbin*3)
            fc%nbw(1:fc%nbin)=realarr(1+fc%nbin*3:fc%nbin*4)
        end if
		fc%xmin=realarr(fc%nbin*6+1)
		fc%xmax=realarr(fc%nbin*6+2)
		fc%nsw=realarr(fc%nbin*6+3)
		fc%xstep=realarr(fc%nbin*6+4)
		!fc%ns=realarr(fc%nbin*10+5)
	end subroutine

#ifdef USE_HDF5	
	subroutine save_s1d_hst_hdf5(s1d, group_id,  tablename)
		use md_hdf5_table
		implicit none
		class(s1d_hst_type)::s1d
		character*(*) tablename
		INTEGER(HID_T) :: group_id      ! Group identifier
        !INTEGER(HSIZE_T), PARAMETER :: nfields  = 11            ! nfields
        type(hdf5_table_type)::htable
		if(s1d%nbin.le.0.or.(.not.allocated(s1d%xb))) return
        

        if(s1d%use_weight)then
            !nrecords=s1d%nbin;
            call htable%init_table(5, s1d%nbin,  tablename)
            htable%field_names(1:5)=(/"   X","  FX"," FXW", "  nb"," nbw"/)	

            htable%field_types(1:3)=H5T_NATIVE_DOUBLE
			htable%field_types(5)=H5T_NATIVE_DOUBLE
            htable%field_types(4)=H5T_NATIVE_INTEGER

            call htable%prepare_write_table(group_id)
			if(.not.allocated(s1d%xb))then
				print*, "s1d not allocated ", trim(adjustl(tablename)), s1d%nbin
				stop
			end if
            CALL htable%write_column_real(s1d%xb)
            CALL htable%write_column_real(s1d%fx)
            CALL htable%write_column_real(s1d%fxw)
			CALL htable%write_column_int(s1d%nb)
            CALL htable%write_column_real(s1d%nbw)
            
        else
            call htable%init_table(3, s1d%nbin,  tablename)
            htable%field_names(1:3)=(/"   X","  FX","  nb"/)	

            htable%field_types(1:2)=H5T_NATIVE_DOUBLE
            htable%field_types(3)=H5T_NATIVE_INTEGER

            call htable%prepare_write_table(group_id)

            CALL htable%write_column_real(s1d%xb)
            CALL htable%write_column_real(s1d%fx)
            CALL htable%write_column_int(s1d%nb)
		endif
	end subroutine
	subroutine read_s1d_hst_hdf5(s1d, group_id, tablename,use_weight)
		use md_hdf5_table
		implicit none
		class(s1d_hst_type)::s1d
		character*(*) tablename
		INTEGER(HID_T) :: group_id      ! Group identifier
        !INTEGER(HSIZE_T), PARAMETER :: nfields  = 11            ! nfields
        type(hdf5_table_type)::htable
        logical use_weight
		if(s1d%nbin.eq.0) return
        if(use_weight)then
            call htable%init_table(5, s1d%nbin,  tablename)
            htable%field_names(1:5)=(/"   X","  FX"," FXW", "  nb"," nbw"/)	

            htable%field_types(1:3)=H5T_NATIVE_DOUBLE
			htable%field_types(5)=H5T_NATIVE_DOUBLE
            htable%field_types(4)=H5T_NATIVE_INTEGER

            call htable%prepare_read_table(group_id)

            CALL htable%read_column_real(s1d%xb)
            CALL htable%read_column_real(s1d%fx)
            CALL htable%read_column_real(s1d%fxw)
            CALL htable%read_column_int(s1d%nb)
            CALL htable%read_column_real(s1d%nbw)
        else
            
			call htable%init_table(3, s1d%nbin,  tablename)
            htable%field_names(1:3)=(/"   X","  FX","  nb"/)	

            htable%field_types(1:2)=H5T_NATIVE_DOUBLE
            htable%field_types(3)=H5T_NATIVE_INTEGER

            call htable%prepare_read_table(group_id)

            CALL htable%read_column_real(s1d%xb)
            CALL htable%read_column_real(s1d%fx)
            CALL htable%read_column_int(s1d%nb)
        end if
		
	end subroutine
#endif
end module
