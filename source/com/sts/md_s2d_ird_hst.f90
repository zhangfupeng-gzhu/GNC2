
module md_s2d_hst_ird_basic_type
	use my_intgl
	!use md_sts_fc
	use md_s2d_ird_basic_type 
	integer,parameter:: f_log=0, f_linear=1
	integer,parameter::dct_x=1,dct_y=2
	
	type,extends(s2d_ird_basic_type):: s2d_hst_ird_basic_type
		integer, allocatable:: nxy(:,:)
		real(8), allocatable:: fxyw(:,:), nxyw(:,:)
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
        procedure::print=>print_s2_hst_ird
        procedure::get_s2d_ird_hst
        procedure::get_s2d_ird_hst_weight        
		procedure::get_s2d_ird_hst_weight_kernel
        generic::get_hst=>get_s2d_ird_hst, get_s2d_ird_hst_weight,get_s2d_ird_hst_weight_kernel
	end type
	!interface get_hst
	!	module procedure get_hst_no_weight
	!	module procedure get_hst_weight
	!end interface
    private::print_s2_hst_ird, get_s2d_ird_hst,get_s2d_ird_hst_weight
    
contains
	subroutine init_s2_hst_ird_basic(this, nx, ny, xmin,xmax, ymin,ymax, use_weight)
		implicit none
		class(s2d_hst_ird_basic_type)::this
		integer nx,ny
		real(8) xmin,xmax, ymin,ymax
		logical,optional:: use_weight
		if(nx==0.or.ny==0)then
			print*, "init_s2d_ird:warnning s2d%nx or ny=0", nx, ny
		end if
        !print*, "1"
		call init_s2d_ird_basic(this,nx,ny, xmin,xmax, ymin, ymax, sts_type_dstr)

        if(present(use_weight))then
			this%use_weight=use_weight
		else
			this%use_weight=.false.
		end if

        if(allocated(this%nxy))then
			deallocate(this%nxy)
			if(this%use_weight)then
				deallocate(this%nxyw, this%fxyw)
			end if
		end if
		
		allocate(this%nxy(nx,ny))
		this%nxy=0
		if(this%use_weight)then
			allocate(this%nxyw(nx,ny), this%fxyw(nx,ny))
			this%nxyw=0; this%fxyw=0;
		end if
        this%ns=0; this%nsw=0
	end subroutine
	
	subroutine get_s2d_ird_hst_weight(s2_hst,x,y,w,n)
		implicit none
		class(s2d_hst_ird_basic_type)::s2_hst
		!type(sts_fc_type)::s2d
		integer n
		real(8) x(n),y(n),w(n)
		integer i,j
		!call init_fc(s2d,s2d%xmin,s2d%xmax,s2d%nbin,fc_spacing_linear)
		if(.not.s2_hst%use_weight)then
			print*, "error! s2_hst%use_weight=FALSE"
			stop
		end if
		!print*, "x=",x(1:10), n, s2_hst%nbin
		call get_dstr_num_in_each_ird_bin_weight_s2d(x(1:n),y(1:n),w(1:n),n,s2_hst%xcenter, &
            s2_hst%xsteps, s2_hst%nx, s2_hst%ycenter, &
            s2_hst%ysteps, s2_hst%ny, s2_hst%nxyw, s2_hst%nsw)
        do i=1, s2_hst%nx
            do j=1, s2_hst%ny
                s2_hst%fxyw(i,j)=dble(s2_hst%nxyw(i,j))/s2_hst%xsteps(i)/s2_hst%ysteps(j)
            end do
        end do

		call get_dstr_num_in_each_ird_bin_s2d(x(1:n),y(1:n),n,s2_hst%xcenter, &
        s2_hst%xsteps, s2_hst%nx, s2_hst%ycenter, &
        s2_hst%ysteps, s2_hst%ny, s2_hst%nxy, s2_hst%ns)
        do i=1, s2_hst%nx
            do j=1, s2_hst%ny
		        s2_hst%fxy(i,j)=dble(s2_hst%nxy(i,j))/s2_hst%xsteps(i)/s2_hst%ysteps(j)
            end do
        end do

	end subroutine
	subroutine get_s2d_ird_hst_weight_kernel(s2_hst,x,y,w,n,kernel_flag)
		implicit none
		class(s2d_hst_ird_basic_type)::s2_hst
		!type(sts_fc_type)::s2d
		integer n
		real(8) x(n),y(n),w(n),xmin,xmax
		integer i,j,kernel_flag
		integer,parameter::kernel_normal=2,kernel_delta=1
		!call init_fc(s2d,s2d%xmin,s2d%xmax,s2d%nbin,fc_spacing_linear)
		if(.not.s2_hst%use_weight)then
			print*, "error! s2_hst%use_weight=FALSE"
			stop
		end if
		select case(kernel_flag)
		case(kernel_delta)
			call get_dstr_num_in_each_ird_bin_weight_s2d_delta_kernel(x(1:n),y(1:n),w(1:n),n,s2_hst%xmin,s2_hst%xmax,s2_hst%xcenter, &
            s2_hst%xsteps, s2_hst%nx, s2_hst%ycenter, &
            s2_hst%ysteps, s2_hst%ny, s2_hst%nxyw, s2_hst%nsw)
			
			do i=1, s2_hst%nx
				do j=1, s2_hst%ny
					s2_hst%fxyw(i,j)=dble(s2_hst%nxyw(i,j))/s2_hst%xsteps(i)/s2_hst%ysteps(j)
				end do
			end do

			call get_dstr_num_in_each_ird_bin_s2d_delta_kernel(x(1:n),y(1:n),n,s2_hst%xmin,s2_hst%xmax,s2_hst%xcenter, &
			s2_hst%xsteps, s2_hst%nx, s2_hst%ycenter, &
			s2_hst%ysteps, s2_hst%ny, s2_hst%nxy, s2_hst%ns)
			do i=1, s2_hst%nx
				do j=1, s2_hst%ny
					s2_hst%fxy(i,j)=dble(s2_hst%nxy(i,j))/s2_hst%xsteps(i)/s2_hst%ysteps(j)
				end do
			end do

		case(kernel_normal)
			call get_dstr_num_in_each_ird_bin_weight_s2d_normal_kernel(x(1:n),y(1:n),w(1:n),n,s2_hst%xmin,s2_hst%xmax,s2_hst%xcenter, &
            s2_hst%xsteps, s2_hst%nx, s2_hst%ycenter, &
            s2_hst%ysteps, s2_hst%ny, s2_hst%nxyw, s2_hst%nsw)
			
			do i=1, s2_hst%nx
				do j=1, s2_hst%ny
					s2_hst%fxyw(i,j)=dble(s2_hst%nxyw(i,j))/s2_hst%xsteps(i)/s2_hst%ysteps(j)
				end do
			end do

			call get_dstr_num_in_each_ird_bin_s2d_normal_kernel(x(1:n),y(1:n),n,s2_hst%xmin,s2_hst%xmax,s2_hst%xcenter, &
			s2_hst%xsteps, s2_hst%nx, s2_hst%ycenter, &
			s2_hst%ysteps, s2_hst%ny, s2_hst%nxy, s2_hst%ns)
			do i=1, s2_hst%nx
				do j=1, s2_hst%ny
					s2_hst%fxy(i,j)=dble(s2_hst%nxy(i,j))/s2_hst%xsteps(i)/s2_hst%ysteps(j)
				end do
			end do
		case default
			print*, "error! define kernel type"
			stop
		end select
	end subroutine
	subroutine get_s2d_ird_hst(s2_hst,x,y,n)
		implicit none
		class(s2d_hst_ird_basic_type)::s2_hst
		!type(sts_fc_type)::s2d
		integer n
		real(8) x(n),y(n)
		integer i,j
		!call init_fc(s2d,s2d%xmin,s2d%xmax,s2d%nbin,fc_spacing_linear)
		if(s2_hst%use_weight)then
			print*, "warnning: s2_hst%use_weight=True, but there is no input weighting data "
		endif
		!w=1d0
        !print*, "start"

		call get_dstr_num_in_each_ird_bin_s2d(x(1:n),y(1:n),n,s2_hst%xcenter, &
        s2_hst%xsteps, s2_hst%nx, s2_hst%ycenter, &
        s2_hst%ysteps, s2_hst%ny, s2_hst%nxy, s2_hst%ns)

		do i=1, s2_hst%nx
            do j=1, s2_hst%ny
		        s2_hst%fxy(i,j)=dble(s2_hst%nxy(i,j))/s2_hst%xsteps(i)/s2_hst%ysteps(j)
            end do
        end do
        !print*, "end"
	end subroutine	

    subroutine print_s2_hst_ird(this,str_)
		implicit none
        class(s2d_hst_ird_basic_type)::this
        character*(15) str_bin_type
		character*(*) , optional::str_
        !integer,optional::pflag_
        integer i,pflag
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
        !if(present(pflag_))then
        !    pflag=pflag_
        !else
            pflag=3
        !end if
		do i=1, this%nx
            select case(pflag)
            case(1)
			    write(*, fmt="(100E12.3)") this%fxy(i,:)
            case(2)
                write(*, fmt="(100I12.3)") this%nxy(i,:)
            case(3)
                write(*, fmt="(100E12.3)") this%nxyw(i,:)
            case default
                print*, "error, define print_s2_hst_ird: pflag=",pflag
                return
            end select
		end do
    end subroutine

end module

module md_s2d_hst_ird_type
	use md_s2d_hst_ird_basic_type
	type,extends(s2d_hst_ird_basic_type):: s2d_hst_ird_type
		integer type_int_size
		integer type_real_size
		integer type_log_size
        !integer type_s2d_size
	contains
		procedure::init=>init_s2d_hst

		procedure::read_s2d_hst
		procedure::write_s2d_hst
#ifdef USE_HDF5		
		procedure::save_hdf5=>save_s2d_hst_hdf5
		!procedure::read_hdf5=>read_s2d_hst_hdf5
#endif
		generic :: read(unformatted) => read_s2d_hst
		generic :: write(unformatted) => write_s2d_hst
		
	end type
	private init_s2d_hst
	private::read_s2d_hst
	private::write_s2d_hst
#ifdef USE_HDF5	
	private::save_s2d_hst_hdf5!, read_s2d_hst_hdf5
#endif
contains 
	subroutine init_s2d_hst(this,nx,ny, xmin,xmax, ymin, ymax, use_weight)
		implicit none
		class(s2d_hst_ird_type)::this
		integer nx,ny
		real(8) xmin,xmax, ymin, ymax
		logical,optional:: use_weight
		if(present(use_weight))then
			call init_s2_hst_ird_basic(this,nx, ny,  xmin,xmax,ymin, ymax, use_weight)
		else
			call init_s2_hst_ird_basic(this,nx, ny,  xmin,xmax,ymin, ymax, .false.)
		end if
		this%type_log_size=2
        this%type_int_size=this%nx*this%ny+3
        this%type_real_size=this%nx*this%ny*3+this%nx+this%ny+5
	end subroutine

	subroutine read_s2d_hst(s2d,file_unit, iostat, iomsg)
		implicit none
		class(s2d_hst_ird_type), intent(inout)::s2d
		integer,intent(in):: file_unit
		integer,intent(out)::iostat
		character(*), intent(inout) :: iomsg
		integer n
		read(unit=file_unit, iostat=iostat, iomsg=iomsg) s2d%nx, s2d%ny, s2d%xmin,s2d%xmax, s2d%ymin,s2d%ymax
		read(unit=file_unit, iostat=iostat, iomsg=iomsg) s2d%use_weight
		if(s2d%nx.eq.0.or.s2d%ny.eq.0)return
		call s2d%init(s2d%nx, s2d%ny, s2d%xmin,s2d%xmax, s2d%ymin,s2d%ymax,s2d%use_weight)
		read(unit=file_unit, iostat=iostat, iomsg=iomsg) s2d%xcenter(1:s2d%nx), s2d%ycenter(1:s2d%ny), &
        s2d%fxy(1:s2d%nx,1:s2d%ny), s2d%nxy(1:s2d%nx, 1:s2d%ny), &
        s2d%xsteps(1:s2d%nx),s2d%ysteps(1:s2d%ny), s2d%ns
		if(s2d%use_weight)then
			read(unit=file_unit, iostat=iostat, iomsg=iomsg)  s2d%nxyw(1:s2d%nx,1:s2d%ny), s2d%fxyw(1:s2d%nx,1:s2d%ny), s2d%nsw
		endif
	end subroutine
	subroutine write_s2d_hst(s2d,file_unit, iostat, iomsg)
		implicit none
		class(s2d_hst_ird_type), intent(in)::s2d
		integer,intent(in):: file_unit
		integer,intent(out)::iostat
		character(*), intent(inout) :: iomsg
		integer n

		write(unit=file_unit, iostat=iostat, iomsg=iomsg) s2d%nx, s2d%ny, s2d%xmin,s2d%xmax, s2d%ymin,s2d%ymax
		write(unit=file_unit, iostat=iostat, iomsg=iomsg) s2d%use_weight
		if(s2d%nx.eq.0.or.s2d%ny.eq.0)return
		write(unit=file_unit, iostat=iostat, iomsg=iomsg) s2d%xcenter(1:s2d%nx), s2d%ycenter(1:s2d%ny), &
        s2d%fxy(1:s2d%nx,1:s2d%ny), s2d%nxy(1:s2d%nx, 1:s2d%ny), &
        s2d%xsteps(1:s2d%nx),s2d%ysteps(1:s2d%ny), s2d%ns
		if(s2d%use_weight)then
			write(unit=file_unit, iostat=iostat, iomsg=iomsg)  s2d%nxyw(1:s2d%nx,1:s2d%ny), s2d%fxyw(1:s2d%nx,1:s2d%ny), s2d%nsw
		end if
	end subroutine
	subroutine conv_s2d_ird_hst_int_real_arrays(fc, intarr, realarr, logarr)
		!use com_main_gw
		implicit none
		type(s2d_hst_ird_type)::fc
		integer nint, nreal
        logical logarr(fc%type_log_size)
		integer intarr(fc%type_int_size)
		real(8) realarr(fc%type_real_size)
        real(8) s2darr(fc%nx, fc%ny)
        integer ibg
        integer i

        logarr=(/fc%use_weight, fc%is_spline_prepared/)

		intarr(1:3)=(/fc%nx,fc%ny, fc%ns/)
        ibg=3
        do i=1, fc%ny
		    intarr(ibg+1:ibg+fc%nx)=fc%nxy(1:fc%nx,i)
            ibg=ibg+fc%nx
        end do
        
        
        realarr(1:fc%nx)=fc%xcenter(1:fc%nx)
        realarr(fc%nx+1:fc%nx+fc%ny)=fc%ycenter(1:fc%ny)
        ibg=fc%nx+fc%ny
        do i=1, fc%ny
            realarr(ibg+1:ibg+fc%nx)=fc%fxy(1:fc%nx,i)
            ibg=ibg+fc%nx
        end do

        if(fc%use_weight)then
            do i=1, fc%ny
                realarr(ibg+1:ibg+fc%nx)=fc%fxyw(1:fc%nx,i)
                ibg=ibg+fc%nx
            end do
            do i=1, fc%ny
                realarr(ibg+1:ibg+fc%nx)=fc%nxyw(1:fc%nx,i)
                ibg=ibg+fc%nx
            end do
		end if        

		realarr(ibg+1:ibg+5)=&
				(/fc%xmin,fc%xmax,fc%ymin, fc%ymin, fc%nsw/)

	end subroutine
	subroutine conv_int_real_arrays_s2d_ird_hst(fc, intarr, realarr,logarr)
		!use com_main_gw
		implicit none
		type(s2d_hst_ird_type)::fc
		integer nint, nreal,ibg
        logical logarr(fc%type_log_size)
		integer intarr(fc%type_int_size)
		real(8) realarr(fc%type_real_size)
        integer i

        fc%use_weight=logarr(1)
        fc%is_spline_prepared=logarr(2)
		fc%nx=intarr(1)
        fc%ny=intarr(2)
        fc%ns=intarr(3)
        ibg=3
        do i=1, fc%ny
		    fc%nxy(1:fc%nx,i)=intarr(ibg+1:ibg+fc%nx)
            ibg=ibg+fc%nx
        end do

		fc%xcenter(1:fc%nx)=realarr(1:fc%nx)
        fc%ycenter(1:fc%ny)=realarr(fc%nx+1:fc%nx+fc%ny)
        ibg=fc%nx+fc%ny
        do i=1, fc%ny
		    fc%nxy(1:fc%nx, i)=realarr(ibg+1:ibg+fc%nx)
            ibg=ibg+fc%nx
        end do
		
        if(fc%use_weight)then
            do i=1, fc%ny
                fc%fxyw(1:fc%nx,i)= realarr(ibg+1:ibg+fc%nx)
                ibg=ibg+fc%nx
            end do
            do i=1, fc%ny
                fc%nxyw(1:fc%nx,i)=realarr(ibg+1:ibg+fc%nx)
                ibg=ibg+fc%nx
            end do
        end if
        fc%xmin=realarr(ibg+1)
		fc%xmax=realarr(ibg+2)
        fc%ymin=realarr(ibg+3)
		fc%ymax=realarr(ibg+4)
		fc%nsw=realarr(ibg+5)
		!fc%xsteps=realarr(fc%nbin*6+4)
		!fc%ns=realarr(fc%nbin*10+5)
	end subroutine

#ifdef USE_HDF5	
subroutine save_s2d_hst_hdf5(s2d, group_id, s2dname)
    use md_hdf5
    implicit none
    class(s2d_hst_ird_type)::s2d
    integer(HID_T):: sub_group_id, group_id
    character*(*) s2dname
    integer error
    INTEGER(HSIZE_T), DIMENSION(2) :: dims2,dimsx,dimsy  ! Dataset dimensions
    real(8) xg(2,s2d%nx), yg(2,s2d%ny)
	if(s2d%nx.le.0.or. s2d%ny.le.0.or.(.not.allocated(s2d%xcenter).or.(.not.allocated(s2d%ycenter))))then
		return
	end if
    call h5gcreate_f(group_id, trim(adjustl(s2dname)), sub_group_id, error)
    if(s2d%nx.ne.size(s2d%xcenter))then
		print*, "save_s2d_hst_hdf5:s2d%nx.ne.size(s2d%xcenter)", s2d%nx, size(s2d%xcenter), "  "//trim(adjustl(s2dname))
		return
	end if
	if(s2d%ny.ne.size(s2d%ycenter))then
		print*, "save_s2d_hst_hdf5:s2d%ny.ne.size(s2d%ycenter)", s2d%nx, size(s2d%ycenter), "  "//trim(adjustl(s2dname))
		return
	end if
	if(s2d%nx.gt.0.and.s2d%ny.gt.0)then
		dims2=(/s2d%nx,s2d%ny/)
		dimsx=(/2,s2d%nx/);dimsy=(/2,s2d%ny/)

		xg(1,:)=s2d%xcenter
		xg(2,:)=s2d%xsteps
		yg(1,:)=s2d%ycenter
		yg(2,:)=s2d%ysteps
		
		call h5ltmake_dataset_double_f(sub_group_id, "X",  2, dimsx,   xg, error)
		call h5ltmake_dataset_double_f(sub_group_id, "Y",  2, dimsy,   yg, error)
		call h5ltmake_dataset_double_f(sub_group_id, "FXY", 2, dims2,  s2d%fxy(:,:), error)
		call h5ltmake_dataset_int_f(sub_group_id, "NXY", 2, dims2,  s2d%nxy(:,:), error)
		if(s2d%use_weight)then
			call h5ltmake_dataset_double_f(sub_group_id, "FXYW", 2, dims2, s2d%fxyw(:,:), error)
			call h5ltmake_dataset_double_f(sub_group_id, "NXYW", 2, dims2, s2d%nxyw(:,:), error)
		end if
	end if
    call h5gclose_f(sub_group_id, error)   
end subroutine

#endif
end module
