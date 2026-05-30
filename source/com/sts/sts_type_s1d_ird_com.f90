module md_s1d_ird_basic_type
    use bin_constants
    use md_bf_type
	type,extends(bin_function_type)::type_s1d_ird_basic
	 		real(8) xmin,xmax
            real(8),allocatable::xsteps(:)
			integer:: bin_type
			!logical,private::is_range_set
	contains
			procedure::print=>print_s1d
!			procedure,public::get_value_d=>get_value_d_s1d_basic
	end type
	private::set_range_s1d, print_s1d!,get_value_d_s1d_basic
contains
    subroutine init_s1d_ird_basic(this,  xmin,xmax,n, bin_type)
        implicit none
        class(type_s1d_ird_basic)::this
        integer n,bin_type
        real(8) xmin,xmax
        !character*(*),optional::s1d_name
        
        if(n==0)then
            print*, "init_s1d_ird_basic:warnning s1d%nbin=0"
        end if
        call init_bin_function(this, n)
        this%xmin=xmin;this%xmax=xmax
        this%bin_type=bin_type
        if(allocated(this%xsteps)) deallocate(this%xsteps)
        allocate(this%xsteps(n))
        if(this%xmin>this%xmax)then
            print*, "init_s1d_ird_basic:warnning xmin>xmax", xmin, xmax
        end if
        !print*, "xmin,xmax=",xmin, xmax
    end subroutine
	subroutine print_s1d(this,str_,print_y2_in)
		implicit none
        class(type_s1d_ird_basic)::this
        character*(15) str_bin_type
		character*(*) , optional::str_
        logical,optional::print_y2_in
        logical print_y2
        integer i
		if(present(str_))then
            write(*,*) "for s1d=", trim(adjustl(str_))
        end if
		select case (this%bin_type)
        case (sts_type_grid)
            str_bin_type="GRID"
        case (sts_type_dstr)
            str_bin_type="DSTR"
        end select

        write(unit=*, fmt="(10A15)") "bin_type", "nbin", "xmin", "xmax"
        write(unit=*, fmt="(A15, I15, 5E15.5)") str_bin_type, this%nbin, &
            this%xmin, this%xmax 

        if(present(print_y2_in))then
            print_y2=print_y2_in
        else
            print_y2=.false.
        end if
        if(.not.print_y2)then
            write(unit=*, fmt="(20A20)") "X", "FX", "XSTEP"
            do i=1, this%nbin
                write(unit=*, fmt="(20E20.10)") this%xb(i), this%fx(i), this%xsteps(i)
            end do
        else
            write(unit=*, fmt="(20A20)") "X", "FX", "XSTEP", "Y2"
            do i=1, this%nbin
                write(unit=*, fmt="(20E20.10)") this%xb(i), this%fx(i), this%xsteps(i), this%y2(i)
            end do
        end if
    end subroutine
end module


module md_s1d_ird_type
    use md_s1d_ird_basic_type
    implicit none
	type,extends(type_s1d_ird_basic)::s1d_ird_type
		integer type_int_size
		integer type_real_size
		integer type_log_size
	contains
		procedure::init=>init_s1d_ird
		procedure,private::read_s1d
		procedure,private::write_s1d
		!procedure::conv_to_array=>conv_s1d_int_real_arrays
		!procedure::conv_from_array=>conv_int_real_arrays_s1d
#ifdef USE_HDF5		
		procedure::save_hdf5=>save_s1d_hdf5
		procedure::read_hdf5=>read_s1d_hdf5
#endif		
		generic::read(unformatted)=>read_s1d
		generic::write(unformatted)=>write_s1d
	end type
	private init_s1d_ird
	private::read_s1d
	private::write_s1d
	!private::conv_s1d_int_real_arrays
	!private::conv_int_real_arrays_s1d
#ifdef USE_HDF5
	private::save_s1d_hdf5
	private::read_s1d_hdf5
#endif

contains
subroutine init_s1d_ird(this,  xmin,xmax,n, bin_type)
    implicit none
    class(s1d_ird_type)::this
    integer n,bin_type
    real(8) xmin,xmax
    call init_s1d_ird_basic(this,  xmin,xmax,n, bin_type)
    this%type_int_size=2
    this%type_real_size=this%nbin*4+2
    this%type_log_size=1
end subroutine
subroutine read_s1d(s1d,file_unit, iostat, iomsg)
    implicit none
    class(s1d_ird_type), intent(inout)::s1d
    integer,intent(in):: file_unit
    integer,intent(out)::iostat
    character(*), intent(inout) :: iomsg
    integer n, bin_type
    read(unit=file_unit) s1d%nbin, s1d%xmin,s1d%xmax, bin_type
    n=s1d%nbin
    call s1d%init(s1d%xmin,s1d%xmax,n,bin_type)
    read(unit=file_unit) s1d%xb(1:n), s1d%fx(1:n),s1d%xsteps(1:n)
end subroutine
subroutine write_s1d(s1d,file_unit, iostat, iomsg)
    implicit none
    class(s1d_ird_type), intent(in)::s1d
    integer,intent(in):: file_unit
    integer,intent(out)::iostat
    character(*), intent(inout) :: iomsg
    integer n

    write(unit=file_unit) s1d%nbin, s1d%xmin,s1d%xmax,s1d%bin_type
    n=s1d%nbin
    write(unit=file_unit) s1d%xb(1:n), s1d%fx(1:n),s1d%xsteps(1:n)
end subroutine
subroutine conv_s1d_ird_int_real_arrays(s1d, intarr, realarr,logarr)
    !use com_main_gw
    implicit none
    class(s1d_ird_type)::s1d
    integer intarr(s1d%type_int_size)
    logical logarr(s1d%type_log_size)
    real(8) realarr(s1d%type_real_size)

    intarr(1:2)=(/s1d%nbin,s1d%bin_type/)
    logarr=(/s1d%is_spline_prepared/)
    realarr(1:s1d%nbin)=s1d%xb(1:s1d%nbin)
    realarr(s1d%nbin+1:2*s1d%nbin)=s1d%fx(1:s1d%nbin)
    if(s1d%is_spline_prepared)then
        realarr(s1d%nbin*2+1:3*s1d%nbin)=s1d%y2(1:s1d%nbin)
    end if
    realarr(s1d%nbin*3+1:s1d%nbin*4)=s1d%xsteps(1:s1d%nbin)
    realarr(s1d%nbin*4+1:4*s1d%nbin+2)=&
            (/s1d%xmin,s1d%xmax/)

end subroutine
subroutine conv_int_real_arrays_s1d_ird(s1d, intarr,  realarr,logarr)
    !use com_main_gw
    implicit none
    class(s1d_ird_type)::s1d
    integer nint, nreal
    integer intarr(s1d%type_int_size)
    real(8) realarr(s1d%type_real_size)
    logical logarr(s1d%type_log_size)
    
    s1d%nbin=intarr(1); s1d%bin_type=intarr(2)
    s1d%is_spline_prepared=logarr(1)
    s1d%xb(1:s1d%nbin)=realarr(1:s1d%nbin)
    s1d%fx(1:s1d%nbin)=realarr(1+s1d%nbin:s1d%nbin*2)
    if(s1d%is_spline_prepared)then
        if(.not.allocated(s1d%y2)) allocate(s1d%y2(s1d%nbin))
        s1d%y2(1:s1d%nbin)=realarr(1+s1d%nbin*2:s1d%nbin*3)
    end if
    s1d%xsteps(1:s1d%nbin)=realarr(s1d%nbin*3+1:s1d%nbin*4)
    s1d%xmin=realarr(s1d%nbin*4+1)
    s1d%xmax=realarr(s1d%nbin*4+2)
end subroutine
#ifdef USE_HDF5	
subroutine save_s1d_hdf5(s1d, group_id,  tablename)
    use md_hdf5_table
    implicit none
    class(s1d_ird_type)::s1d
    character*(*) tablename
    INTEGER(HID_T) :: group_id      ! Group identifier
    !INTEGER(HSIZE_T), PARAMETER :: nfields  = 11            ! nfields
    type(hdf5_table_type)::htable

    if(s1d%nbin.eq.0) return

    call htable%init_table(3, s1d%nbin,  tablename)        
    htable%field_names=(/"   X","  FX", "  DX"/)	

    htable%field_types(1:3)=H5T_NATIVE_DOUBLE

    call htable%prepare_write_table(group_id)
    CALL htable%write_column_real(s1d%xb)
    CALL htable%write_column_real(s1d%fx)
	CALL htable%write_column_real(s1d%xsteps)

end subroutine

subroutine read_s1d_hdf5(s1d, group_id, tablename)
    use md_hdf5_table
    implicit none
    class(s1d_ird_type)::s1d
    character*(*) tablename
    INTEGER(HID_T) :: group_id      ! Group identifier
    !INTEGER(HSIZE_T), PARAMETER :: nfields  = 11            ! nfields
    type(hdf5_table_type)::htable

    if(s1d%nbin.eq.0) return
    call htable%init_table(2, s1d%nbin,  tablename)
    
    !nrecords=s1d%nbin;
    htable%field_names=(/"   X","  FX", "  DX" /)	

    htable%field_types(1:2)=H5T_NATIVE_DOUBLE 

    call htable%prepare_read_table(group_id)

    CALL htable%read_column_real(s1d%xb)
    CALL htable%read_column_real(s1d%fx) 
	CALL htable%read_column_real(s1d%xsteps)
end subroutine
#endif
end module


