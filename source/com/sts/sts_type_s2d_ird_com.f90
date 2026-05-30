!module md_sts
!	type(stats_1d_type)
!		real(8),allocatable::bin(:),fprb(:),fcmp(:),fncp(:)
!		integer,allocatable::ffra(:),fcfr(:),fncf(:)
!		real(8) mean,std
!	end type
!end module

module md_s2d_ird_basic_type
    !use md_bf_type
	use bin_constants
	type ::s2d_ird_basic_type
			real(8) xmin,xmax,ymin,ymax
			real(8),allocatable::xsteps(:),ysteps(:)
			integer:: bin_type
			logical:: is_spline_prepared
			real(8),allocatable:: xcenter(:),ycenter(:),fxy(:,:)
			integer nx,ny
			real(8),allocatable,private::y2(:,:)
	contains
			!procedure::set_range=>set_range_s2d
			procedure::print=>print_s2d_ird
			!procedure::get_bin_type
			procedure::prepare_spline=>prepare_spline_s2d_ird
			procedure::get_value_d=>get_value_d_s2d_ird
			procedure::get_value_s=>get_value_s_s2d_ird
			procedure::get_value_l=>get_value_l_s2d_ird
	end type
	
	private::print_s2d_ird, prepare_spline_s2d_ird
	private::get_value_d_s2d_ird, get_value_s_s2d_ird
contains
	subroutine init_s2d_ird_basic(s2d,nx,ny,xmin,xmax,ymin,ymax,bin_type)
		implicit none
		class(s2d_ird_basic_type)::s2d
		integer nx,ny,bin_type
		integer i
		real(8) xmin,xmax,ymin,ymax
		if(nx.eq.0.or.ny.eq.0)then
			print*, "init_s2d: warnning s2d%nx=0 or s2d%ny=0"
		end if
		if (allocated(s2d%xcenter))then
			!print*, "deallocating..."
			deallocate(s2d%xcenter,s2d%ycenter, s2d%fxy,s2d%xsteps,s2d%ysteps)
		end if
		s2d%nx=nx;s2d%ny=ny
		allocate(s2d%xcenter(s2d%nx))
		allocate(s2d%ycenter(s2d%ny))
		allocate(s2d%fxy(s2d%nx,s2d%ny))
		allocate(s2d%xsteps(s2d%nx),s2d%ysteps(s2d%ny))
		s2d%xmin=xmin;s2d%xmax=xmax;s2d%ymin=ymin;s2d%ymax=ymax
		if(s2d%xmax<s2d%xmin)then
			print*, "error! s2d%xmax<s2d%xmin", s2d%xmax,s2d%xmin
			stop
		end if
		s2d%bin_type=bin_type
		s2d%fxy=0; s2d%xcenter=0; s2d%ycenter=0
		s2d%is_spline_prepared=.false.
	end subroutine
	subroutine prepare_spline_s2d_ird(s2d)
		implicit none
		class(s2d_ird_basic_type)::s2d
		if(allocated(s2d%y2))deallocate(s2d%y2)
		allocate(s2d%y2(s2d%nx,s2d%ny))

		call splie2(s2d%xcenter,s2d%ycenter,s2d%fxy,s2d%nx,s2d%ny,s2d%y2)	
		s2d%is_spline_prepared=.true.
	end subroutine
	
	integer function get_bin_type(this)
		implicit none
		class(s2d_ird_basic_type)::this
		get_bin_type=this%bin_type
	end function

	subroutine print_s2d_ird(this,str_)
		implicit none
        class(s2d_ird_basic_type)::this
        character*(15) str_bin_type
		character*(*) , optional::str_
        integer i
		if(present(str_))then
            write(*,*) "s2d=", trim(adjustl(str_))
        end if
		select case (this%bin_type)
        case (sts_type_grid)
            str_bin_type="GRID"
        case (sts_type_dstr)
            str_bin_type="DSTR"
        end select

        write(unit=*, fmt="(10A15)") "bin_type", "nx", "xmin", "xmax", &
			 "ny", "ymin", "ymax"
        write(unit=*, fmt="(A15, 2(I15, 2E15.5))") str_bin_type, this%nx,&
            this%xmin, this%xmax, this%ny, this%ymin,this%ymax
        
		do i=1, this%nx
			write(*, fmt="(100E12.3)") this%fxy(i,:)
		end do
    end subroutine
	subroutine get_value_s_s2d_ird(s2d,vx,vy,yout)
		implicit none
		class(s2d_ird_basic_type)::s2d
		real(8) vx,vy, yout
		integer idx,idy
		if(.not.s2d%is_spline_prepared)then
			print*, "s2d error: spline table is not prepared"
			stop
		end if
		if(s2d%xmax<vx.or.s2d%xmin>vx.or. s2d%ymax<vy.or.s2d%ymin>vy)then
			yout=0d0
		else
!				select case(fscale)
!				case(f_linear)
				call splin2(s2d%xcenter,s2d%ycenter,s2d%fxy,s2d%y2,s2d%nx,s2d%ny,vx,vy, yout)
!				case(f_log)
!					call splin2(s2d%xcenter,s2d%ycenter,log10(s2d%fxy),s2d%y2,s2d%nx,s2d%ny,vx,vy, y_out)
!				end select
		end if
	end subroutine
	subroutine get_value_d_s2d_ird(s2d, vx, vy, yout)
		implicit none
		class(s2d_ird_basic_type)::s2d
		real(8) vx,vy, yout
		integer idx,idy

		if(s2d%xmax<vx.or.s2d%xmin>vx.or. s2d%ymax<vy.or.s2d%ymin>vy)then
			yout=0d0
		else
			call return_idxy(vx,vy,s2d%xmin,s2d%xmax,s2d%ymin,s2d%ymax,s2d%nx,s2d%ny,idx,idy,s2d%bin_type)
			if(idx>0.and.idx<=s2d%nx.and.idy>0.and.idy<=s2d%ny)then
				yout=s2d%fxy(idx,idy)
			else
				yout=0d0
			end if
		end if
	end subroutine
	subroutine get_value_l_s2d_ird(s2d, vx, vy, yout)
		implicit none
		class(s2d_ird_basic_type)::s2d
		real(8) vx,vy, yout, rdx,rdy
		integer idx,idy

		if(s2d%xmax<vx.or.s2d%xmin>vx.or. s2d%ymax<vy.or.s2d%ymin>vy)then
			yout=0d0
		else
			!print*, "finished the code: get_value_l_s2d in ird"
			!stop
			call linear_int_2d_xy(idx,idy,rdx,rdy, s2d%fxy, s2d%nx,s2d%ny, yout)
		end if
	end subroutine

end module



module md_s2d_ird_type 
    use md_s2d_ird_basic_type
	use bin_constants
	type,extends(s2d_ird_basic_type) ::s2d_ird_type
		integer type_int_size
		integer type_real_size
		integer type_log_size
	contains
			procedure::init=>init_s2d_ir
			procedure::read_s2d_ir
			procedure::write_s2d_ir
#ifdef USE_HDF5		
			procedure::save_hdf5=>save_s2d_ir_hdf5
			procedure::read_hdf5=>read_s2d_ir_hdf5
			procedure::output_hdf5=>output_s2d_ir_hdf5
#endif								
			generic ::read(unformatted)=>read_s2d_ir
			generic ::write(unformatted)=>write_s2d_ir
	end type
	private::init_s2d
#ifdef USE_HDF5
	private::save_s2d_irhdf5, write_s2d_ir
	private::read_s2d_irhdf5, read_s2d_ir
#endif	
contains
	subroutine init_s2d_ir(s2d,nx,ny,xmin,xmax,ymin,ymax,bin_type)
		implicit none
		class(s2d_ird_type)::s2d
		integer nx,ny,bin_type
		integer i
		real(8) xmin,xmax,ymin,ymax

		call init_s2d_ird_basic(s2d,nx,ny,xmin,xmax,ymin,ymax,bin_type)
		s2d%type_int_size=3
		s2d%type_real_size=s2d%nx+s2d%ny+6
		s2d%type_log_size=1
	end subroutine

	subroutine read_s2d_ir(s2d,file_unit, iostat, iomsg)
		implicit none
		class(s2d_ird_type),intent(inout)::s2d
		integer nx,ny,i, bin_type
		integer,intent(in):: file_unit
		integer, intent(out) :: iostat
		character(*), intent(inout) :: iomsg

		read(unit=file_unit) nx, ny,s2d%xmin,s2d%xmax,s2d%ymin,s2d%ymax, bin_type
		call s2d%init(nx,ny,s2d%xmin,s2d%xmax,s2d%ymin,s2d%ymax, bin_type)
		read(unit=file_unit) s2d%xcenter(1:nx), s2d%ycenter(1:ny), s2d%fxy(1:nx,1:ny)
		read(unit=file_unit) s2d%xsteps, s2d%ysteps
	end subroutine

	subroutine write_s2d_ir(s2d,file_unit, iostat, iomsg)
		implicit none
		class(s2d_ird_type),intent(in)::s2d
		integer nx,ny,i
		integer,intent(in):: file_unit
		integer, intent(out) :: iostat
		character(*), intent(inout) :: iomsg

		write(unit=file_unit) s2d%nx,s2d%ny,s2d%xmin,s2d%xmax,s2d%ymin,s2d%ymax, s2d%bin_type
		nx=s2d%nx;	ny=s2d%ny;
		write(unit=file_unit) s2d%xcenter(1:nx), s2d%ycenter(1:ny), s2d%fxy(1:nx,1:ny)
		write(unit=file_unit) s2d%xsteps, s2d%ysteps
		
	end subroutine
	subroutine conv_s2d_int_real_arrays_ird(s2d, intarr, realarr,s2darr,logarr)
		!use com_main_gw
		implicit none
		class(s2d_ird_type)::s2d
		integer intarr(s2d%type_int_size)
		logical logarr(s2d%type_log_size)
		real(8) s2darr(s2d%nx,s2d%ny)
		real(8) realarr(s2d%type_real_size)
		integer::counter

		intarr(1:s2d%type_int_size)=(/s2d%nx,s2d%ny, s2d%bin_type/)
		logarr=(/s2d%is_spline_prepared/)
		s2darr=s2d%fxy
		counter=0
		!print*, s2d%type_real_size
		realarr(counter+1:counter+s2d%nx)=s2d%xcenter(1:s2d%nx)
		counter=counter+s2d%nx
		realarr(counter+1:counter+s2d%ny)=s2d%ycenter(1:s2d%ny)
		counter=counter+s2d%ny
		realarr(counter+1:counter+6)=(/s2d%xmin,s2d%xmax,s2d%ymin,&
			s2d%ymax,s2d%xsteps,s2d%ysteps/)

	end subroutine

	subroutine conv_int_real_arrays_s2d_ird(s2d, intarr, realarr,s2darr, logarr)
		!use com_main_gw
		implicit none
		class(s2d_ird_type)::s2d
		integer intarr(s2d%type_int_size)
		logical logarr(s2d%type_log_size)
		real(8) realarr(s2d%type_real_size)
		real(8) s2darr(s2d%nx,s2d%ny)
		integer::counter
		s2d%nx=intarr(1)
		s2d%ny=intarr(2)
		s2d%bin_type=intarr(3)

		s2d%is_spline_prepared=logarr(1)
		s2d%fxy=s2darr
		counter=0
		s2d%xcenter(1:s2d%nx)=realarr(counter+1:counter+s2d%nx)
		counter=counter+s2d%nx
		s2d%ycenter(1:s2d%ny)=realarr(counter+1:counter+s2d%ny)
		counter=counter+s2d%ny
		s2d%xmin=realarr(counter+1)
		s2d%xmax=realarr(counter+2)
		s2d%ymin=realarr(counter+3)
		s2d%ymax=realarr(counter+4)
		s2d%xsteps=realarr(counter+5)
		s2d%ysteps=realarr(counter+6)
	end subroutine

#ifdef USE_HDF5
	subroutine save_s2d_ir_hdf5(s2d,group_id, s2dname)
		use hdf5
		use h5lt	
		implicit none
		character*(*) s2dname
		INTEGER(HID_T) :: group_id,sub_group_id      ! Group identifier
		INTEGER(HSIZE_T), DIMENSION(2) :: dims2  ! Dataset dimensions
		INTEGER(HSIZE_T), DIMENSION(1) :: dimsx,dimsy  ! Dataset dimensions
		integer error
		class(s2d_ird_type)::s2d
		!print*,allocated(s2d%fxy)
		if(allocated(s2d%fxy))then
			call h5gcreate_f(group_id, trim(adjustl(s2dname)), sub_group_id, error)

			dims2=(/s2d%nx,s2d%ny/)
			dimsx=s2d%nx;dimsy=s2d%ny
			!print*, "1"
			!print*, "allocated?", allocated(s2d%xcenter), allocated(s2d%ycenter), allocated(s2d%fxy)
			call h5ltmake_dataset_double_f(sub_group_id, "X",  1, dimsx,   s2d%xcenter, error)
			call h5ltmake_dataset_double_f(sub_group_id, "Y",  1, dimsy,   s2d%ycenter, error)
			!print*, "2"
			call h5ltmake_dataset_double_f(sub_group_id, "FXY", 2, dims2,  s2d%fxy(:,:), error)
			call h5gclose_f(sub_group_id, error) 
		end if
	end subroutine
	
	subroutine read_s2d_ir_hdf5(s2d,group_id, s2dname)
		use hdf5
		use h5lt	
		implicit none
		character*(*) s2dname
		INTEGER(HID_T) :: group_id,sub_group_id      ! Group identifier
		INTEGER(HSIZE_T), DIMENSION(2) :: dims2  ! Dataset dimensions
		INTEGER(HSIZE_T), DIMENSION(1) :: dimsx,dimsy  ! Dataset dimensions
		integer error
		class(s2d_ird_type)::s2d	

		call h5gcreate_f(group_id, trim(adjustl(s2dname)), sub_group_id, error)

		dims2=(/s2d%nx,s2d%ny/)
		dimsx=s2d%nx;dimsy=s2d%ny
		
		call h5ltread_dataset_f(sub_group_id, "X",H5T_NATIVE_DOUBLE,  s2d%xcenter, dimsx,  error)
		call H5LTread_dataset_f(sub_group_id, "Y",H5T_NATIVE_DOUBLE,   s2d%ycenter, dimsy,  error)
		call H5LTread_dataset_f(sub_group_id, "FXY",H5T_NATIVE_DOUBLE,  s2d%fxy(:,:), dims2, error)
		call h5gclose_f(sub_group_id, error) 
	end subroutine
	subroutine output_s2d_ir_hdf5(s2d,s2dname, fl)
		use md_hdf5
		implicit none
		type(hdf5_file_type)::hf
		class(s2d_ird_type)::s2d
		character*(*) s2dname, fl
		call hf%open(trim(adjustl(fl))//".hdf5")
		call s2d%save_hdf5(hf%file_id,s2dname)
		call hf%close()
	end subroutine
#endif

end module
