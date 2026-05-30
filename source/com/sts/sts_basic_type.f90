
module bin_constants
	implicit none
	integer,parameter::sts_type_grid=1, sts_type_dstr=0, sts_type_dstr_log=2
	integer:: sts_type_default=sts_type_grid
	integer,parameter::method_intp_s1d=1,method_direct_s1d=2, method_linear_s1d=3
    integer,parameter::method_linear_log_s1d=4
end module
module md_bf_type
	implicit none
    type bin_function_type
        real(8),allocatable:: xb(:), fx(:)
		real(8),allocatable::y2(:)
        integer:: nbin
		logical::is_spline_prepared
        !character*(10) nname
        contains
        	!procedure,public::init_bin_function_type
			procedure,public::prepare_spline=>prepare_spline_d2
			procedure,public::get_value_l=>get_y_at_x_linear_d2
			procedure,public::get_value_s=>get_y_at_x_spline_d2
			!procedure,public::get_value_d=>get_y_at_x_direct_d2
        	procedure,public::form_d2
			procedure,public::print=>print_d2
			procedure,public::deallocate=>deallocate_d2
			procedure,public::sort=>sort_bf
    end type
	!interface bin_function_type
	!	module procedure initialize
	!end interface
	private prepare_spline_d2, get_y_at_x_linear_d2,get_y_at_x_spline_d2
	private print_d2, deallocate_d2!,get_y_at_x_direct_d2
	!init_bin_function_type
contains
    subroutine init_bin_function(this, n)
        implicit none
        class(bin_function_type)::this
        integer i,n        
        call deallocate_d2(this)
        this%nbin=n;
        allocate(this%xb(n), this%fx(n))
        this%xb=0; this%fx=0;
		!this%y2=0
		!bf=this
		this%is_spline_prepared=.false.
	end subroutine

	subroutine sort_bf(this,ind_col)
        !use stdlib_sorting 
        implicit none
        class(bin_function_type)::this
        integer ind_col
        integer,parameter:: ind_xb=0,ind_fx=1
        integer indx(this%nbin)

        select case(ind_col)
        case(ind_xb)
            call indexx(this%nbin,this%xb, indx)
        case(ind_fx)
            call indexx(this%nbin,this%fx, indx)
        end select
       ! print*, "indx=",indx
        call sort_arr_index(this%nbin,this%xb,indx)
        call sort_arr_index(this%nbin,this%fx,indx)

    end subroutine

	subroutine deallocate_d2(this)
		implicit none
		class(bin_function_type)::this
		!integer n
		if(allocated(this%fx))then
            deallocate(this%xb, this%fx)
        end if
	end subroutine
    subroutine print_d2(this,str_,print_y2_in)
		implicit none 
		class(bin_function_type)::this
		character*(*) , optional::str_
		logical,optional::print_y2_in
		logical print_y2
		integer i

		if(present(str_))then
            write(*,*) "Name=", trim(adjustl(str_))
        end if

		if(present(print_y2_in))then
			print_y2=print_y2_in
		else
			print_y2=.false.
		end if
		if(.not.print_y2)then
			write(unit=*, fmt="(20A20)") "X", "Y"
			do i=1, this%nbin
				write(unit=*, fmt="(20E20.10)") this%xb(i), this%fx(i)
			end do
		else
			write(unit=*, fmt="(20A20)") "X", "Y", "Y2"
			do i=1, this%nbin
				write(unit=*, fmt="(20E20.10)") this%xb(i), this%fx(i),this%y2(i)
			end do
		endif

	end subroutine	

    subroutine form_d2(this, x, y, n, reverse)
		implicit none
		class(bin_function_type)::this
		integer n
		real(8),intent(in):: x(n), y(n)
		integer, optional::reverse
		integer reverse_in
		if(present(reverse))then
			reverse_in=1
		else 
			reverse_in=0
		end if

		call init_bin_function(this,n)
		!print*, y(1:n)
		if(reverse_in.eq.0)then
			this%xb(1:n)=x(1:n); this%fx(1:n)=y(1:n)
		else
			this%xb(1:n)=x(n:1:-1); this%fx(1:n)=y(n:1:-1)
		end if
		!call this%print_s2()
		!read(*,*)
	end subroutine   
    subroutine prepare_spline_d2(this,yp1_in,ypn_in)
		implicit none
		class(bin_function_type)::this
		real(8)::yp1=1d30,ypn=1d30
		real(8),optional:: yp1_in, ypn_in
		if(present(yp1_in))then
			yp1=yp1_in
		end if
		if(present(ypn_in))then
			ypn=ypn_in
		end if

		if(allocated(this%y2))deallocate(this%y2)
		allocate(this%y2(this%nbin))
			if(this%nbin.eq.1)then
				print*, "prepare_spline_d2:???, this%nbin=1"
				stop
			end if
			if(.not.allocated(this%xb))then
				print*, "??? this%xb not allocated"
				stop
			end if
		  if(this%xb(2)<this%xb(1))then
				print*, "d2:in prepare_spline"
			   	print*,"x is ascending,x(2)<x(1), x(:)=", this%fx
			   	print*,"y(:)=", this%fx
				pause
		  end if
		call spline_mylib(this%xb,this%fx,this%nbin,yp1,ypn,this%y2)
		this%is_spline_prepared=.true.
	end subroutine

	

	subroutine get_y_at_x_linear_d2(this, x, y, expolate)
		implicit none
		class(bin_function_type)::this
		integer n
		real(8) x, y
		logical,optional::expolate
		logical::expolate_in	

		if(present(expolate))then
			expolate_in=expolate
		else
			expolate_in=.false.
		end if
		if(.not.expolate_in)then
			if(this%xb(1)>x)then
				y=this%fx(1)
				return
			else if(this%xb(this%nbin)<x)then
				y=this%fx(this%nbin)
				return	
			end if
		end if		
		n=this%nbin
		call linear_int(this%xb(1:n),this%fx(1:n),this%nbin, x,y)
	
	end subroutine

!	subroutine get_y_at_x_direct_d2(s1d,x,y)
!		class(bin_function_type)::s1d
!		real(8) x,y
!		integer i
!		do i=1, s1d%nbin-1
!			if(x(i)<x.and.x(i+1)>x)then
!			end if
!		end do
!	end subroutine

	subroutine get_y_at_x_spline_d2(this, x, y)
		use ,intrinsic::ieee_arithmetic
		implicit none
		class(bin_function_type)::this
		real(8) x, y
		real(8) xmin,xmax
		integer n,ier
		if(.not.this%is_spline_prepared)then
			call prepare_spline_d2(this)
			!print*, "spline is prepared"
		endif
		!if(this%xmin>x)then
		!	y=this%fx(1)
		!else if(this%xmax<x)then
		!	y=this%fx(this%nbin)
		!else
			n=this%nbin
			call splint_mylib(this%xb(1:n),this%fx(1:n),this%y2(1:n),this%nbin, x, y,ier)
			!print*, "get_y_at_x_spline_d2:x,y=", x,y
			if(ier.ne.0)then
				print*, "stop"
				stop
			end if
			!call linear_int(s1d%xb(1:n),s1d%fx(1:n),s1d%nbin, va,yout)
		!end if
	end subroutine 
	subroutine get_dfdx_d2(s1d, dfdy)
		implicit none
		type(bin_function_type)::s1d,dfdy
		dfdy=s1d
		call get_dfdx(s1d%xb, s1d%fx, s1d%nbin, dfdy%fx)
	end subroutine

end module

module md_s1d_basic_type
	use bin_constants
    use md_bf_type
	type,extends(bin_function_type)::s1d_basic_type
	 		real(8) xmin,xmax,xstep
			integer:: bin_type
			!logical,private::is_range_set
	contains
			procedure::set_range=>set_range_s1d
			!procedure::init=>init_s1d
			procedure::print=>print_s1d
			procedure,public::get_value_d=>get_value_d_s1d_basic
			procedure::get_dfdx_spline_fx
			procedure::get_dfdx_spline_s1d
			generic::get_dfdx_spline=>get_dfdx_spline_fx,get_dfdx_spline_s1d
			!procedure::get_bin_type
			!procedure::get_if_spline_prepared
	end type
	!interface s1d_basic_type
	!	module procedure init_s1d
	!end interface
	private::get_dfdx_spline_s1d,get_dfdx_spline_fx
	
	private::set_range_s1d, print_s1d,get_value_d_s1d_basic
contains
	subroutine init_s1d_basic(this,  xmin,xmax,n, bin_type)
		implicit none
		class(s1d_basic_type)::this
		integer n,bin_type
		real(8) xmin,xmax
		!character*(*),optional::s1d_name
		
		if(n==0)then
			print*, "init_s1d:warnning s1d%nbin=0"
		end if
		call init_bin_function(this, n)
		this%xmin=xmin;this%xmax=xmax
		this%bin_type=bin_type
		if(this%xmin>this%xmax)then
			print*, "warnning: xmin>xmax", xmin, xmax
		end if
		!print*, "xmin,xmax=",xmin, xmax
	end subroutine
	subroutine set_range_s1d(this)
		implicit none
		class(s1d_basic_type)::this
		
		call set_range(this%xb,this%nbin,this%xmin,this%xmax,this%bin_type)
		!print*, this%xb(1:2),this%nbin,this%xmin,this%xmax,this%bin_type
		this%xstep=this%xb(2)-this%xb(1)
	end subroutine
	subroutine print_s1d(this,str_,print_y2_in)
		implicit none
        class(s1d_basic_type)::this
        character*(15) str_bin_type
		character*(*) , optional::str_
		logical,optional::print_y2_in
		logical::print_y2
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

        write(unit=*, fmt="(10A15)") "bin_type", "nbin", "xmin", "xmax", "xstep"
        write(unit=*, fmt="(A15, I15, 5E15.5)") str_bin_type, this%nbin, &
            this%xmin, this%xmax, this%xstep
		if(present(print_y2_in))then
			print_y2=print_y2_in
		else
			print_y2=.false.
		end if
		if(.not.print_y2)then
			write(unit=*, fmt="(20A20)") "X", "FX"
			do i=1, this%nbin
				write(unit=*, fmt="(20E20.10)") this%xb(i), this%fx(i)
			end do
		else
			write(unit=*, fmt="(20A20)") "X", "FX", "Y2"
			do i=1, this%nbin
				write(unit=*, fmt="(20E20.10)") this%xb(i), this%fx(i),this%y2(i)
			end do
		endif

    end subroutine
	subroutine get_value_d_s1d_basic(s1d,x,y)
		implicit none
		class(s1d_basic_type)::s1d
		real(8) x, y
		integer idx
		call return_idx(x,s1d%xmin,s1d%xmax,s1d%nbin,idx,s1d%bin_type)
		if(idx.lt.1.or.idx.gt.s1d%nbin)then
			y=0
		else
			y=s1d%fx(idx)
		end if
	end subroutine
	subroutine get_dir_sum_s1d(s1d, sums)
		implicit none
		class(s1d_basic_type)::s1d
		integer i
		real(8) sums
		sums=0
		select case (s1d%bin_type)
		case(sts_type_dstr)
			do i=1, s1d%nbin
				sums=sums+s1d%xb(i)*s1d%xstep
			end do		
		case(sts_type_grid)
		!--------trapezoidal rule, see Numerical recipe
			sums=sum(s1d%xb(2:))*s1d%xstep
			!print*, sums
			!sums=0
			!do i=2,n
			!	sums=sums+x(i)*dx
			!end do
			!print*, sums
			!pause
			sums=sums+(s1d%xb(1)+s1d%xb(s1d%nbin))/2d0*s1d%xstep
		end select

	end subroutine
	!integer function get_bin_type(this)
	!	implicit none
	!	class(s1d_basic_type)::this
	!	get_bin_type=this%bin_type
	!end function
	!logical function get_if_spline_prepared(this)
	!	implicit none
	!	class(s1d_basic_type)::this
	!	get_if_spline_prepared=this%is_spline_prepared
	!end function
	
	subroutine get_dfdx_spline_s1d(this,s1d)
		implicit none
		class(s1d_basic_type)::this
		class(s1d_basic_type)::s1d		
		integer i
		call this%prepare_spline()
		do i=1,s1d%nbin
			call get_dfdx_spline_atx(this%xb,this%fx,this%y2,this%nbin,s1d%xb(i),s1d%fx(i))
		end do
	end subroutine
	subroutine get_dfdx_spline_fx(this,fx)
		implicit none
		class(s1d_basic_type)::this
		real(8) fx(this%nbin)	
		integer i
		external::get_dfdx_spline
		call this%prepare_spline()
		call get_dfdx_spline(this%xb,this%fx,this%nbin,fx)
	end subroutine

end module


module md_s2d_basic_type
    !use md_bf_type
	use bin_constants
	type ::s2d_basic_type
			real(8) xmin,xmax,ymin,ymax,xstep,ystep
			integer:: bin_type
			logical:: is_spline_prepared
			real(8),allocatable:: xcenter(:),ycenter(:),fxy(:,:)
			integer nx,ny
			real(8),allocatable,private::y2(:,:)
	contains
			procedure::set_range=>set_range_s2d
			procedure::print=>print_s2d
			!procedure::get_bin_type
			procedure::prepare_spline=>prepare_spline_s2d
			procedure::get_value_d=>get_value_d_s2d
			procedure::get_value_s=>get_value_s_s2d
			procedure::get_value_l=>get_value_l_s2d
	end type
	
	private::set_range_s2d, print_s2d, prepare_spline_s2d
	private::get_value_d_s2d, get_value_s_s2d
contains
	subroutine init_s2d_basic(s2d,nx,ny,xmin,xmax,ymin,ymax,bin_type)
		implicit none
		class(s2d_basic_type)::s2d
		integer nx,ny,bin_type
		integer i
		real(8) xmin,xmax,ymin,ymax
		if(nx.eq.0.or.ny.eq.0)then
			print*, "init_s2d: warnning s2d%nx=0 or s2d%ny=0"
		end if
		if (allocated(s2d%xcenter))then
			!print*, "deallocating..."
			deallocate(s2d%xcenter,s2d%ycenter, s2d%fxy)
		end if
		s2d%nx=nx;s2d%ny=ny
		allocate(s2d%xcenter(s2d%nx))
		allocate(s2d%ycenter(s2d%ny))
		allocate(s2d%fxy(s2d%nx,s2d%ny))
		s2d%xmin=xmin;s2d%xmax=xmax;s2d%ymin=ymin;s2d%ymax=ymax
		if(s2d%xmax<s2d%xmin)then
			print*, "error! s2d%xmax<s2d%xmin", s2d%xmax,s2d%xmin
			stop
		end if
		s2d%bin_type=bin_type
		s2d%fxy=0; s2d%xcenter=0; s2d%ycenter=0
		s2d%is_spline_prepared=.false.
	end subroutine
	subroutine prepare_spline_s2d(s2d)
		implicit none
		class(s2d_basic_type)::s2d
		if(allocated(s2d%y2))deallocate(s2d%y2)
		allocate(s2d%y2(s2d%nx,s2d%ny))

		call splie2(s2d%xcenter,s2d%ycenter,s2d%fxy,s2d%nx,s2d%ny,s2d%y2)	
		s2d%is_spline_prepared=.true.
	end subroutine
	
	integer function get_bin_type(this)
		implicit none
		class(s2d_basic_type)::this
		get_bin_type=this%bin_type
	end function
	subroutine set_range_s2d(this)
		implicit none
		class(s2d_basic_type)::this
		call set_range(this%xcenter,this%nx,this%xmin,this%xmax,this%bin_type)
		call set_range(this%ycenter,this%ny,this%ymin,this%ymax,this%bin_type)
		this%xstep=this%xcenter(2)-this%xcenter(1)
		this%ystep=this%ycenter(2)-this%ycenter(1)
	end subroutine
	subroutine print_s2d(this,str_)
		implicit none
        class(s2d_basic_type)::this
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

        write(unit=*, fmt="(10A15)") "bin_type", "nx", "xmin", "xmax", "xstep", &
			 "ny", "ymin", "ymax", "ystep"
        write(unit=*, fmt="(A15, 2(I15, 3E15.5))") str_bin_type, this%nx,&
            this%xmin, this%xmax, this%xstep,  this%ny, this%ymin,this%ymax, this%ystep
        
		do i=1, this%nx
			write(*, fmt="(100E12.3)") this%fxy(i,:)
		end do
    end subroutine
	subroutine get_value_s_s2d(s2d,vx,vy,yout)
		implicit none
		class(s2d_basic_type)::s2d
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
	subroutine get_value_d_s2d(s2d, vx, vy, yout)
		implicit none
		class(s2d_basic_type)::s2d
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
	subroutine get_value_l_s2d(s2d, vx, vy, yout)
		implicit none
		class(s2d_basic_type)::s2d
		real(8) vx,vy, yout
		integer idx,idy

		if(s2d%xmax<vx.or.s2d%xmin>vx.or. s2d%ymax<vy.or.s2d%ymin>vy)then
			yout=0d0
		else
			call linear_int_2d(s2d%xmin,s2d%ymin,s2d%nx,s2d%ny,s2d%xstep,s2d%ystep,s2d%fxy,vx,vy,yout)
		end if
	end subroutine

end module
