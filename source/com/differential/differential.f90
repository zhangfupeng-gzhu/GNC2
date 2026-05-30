subroutine get_dfdx(x, y, n, dydx)
	implicit none
!	external::spline_mylib, splint_mylib
	integer n,i,ier
	real(8) x(n), y(n), dydx(n)
	real(8) y2(n), yp1, ypn,y_1,y_2
	real(8) h, bg, bt
	yp1=1d60;ypn=1d60

	call spline_mylib(x,y, n,yp1,ypn,y2)

	do i=1, n
		if(i/=n .and. i/=1)then
			h=(x(i+1)-x(i))/50d0
			bg=h/2d0; bt=h/2d0
		else if(i==1)then
			h=(x(i+1)-x(i))/50d0
			bg=h; bt=0

		else if(i==n)then
			h=(x(i)-x(i-1))/50d0
			bg=0; bt=h
		end if

		call splint_mylib(x,y,y2, n, x(i)+bg, y_1,ier)
		call splint_mylib(x,y,y2, n, x(i)-bt, y_2,ier)
		dydx(i)=(y_1-y_2)/(bg+bt)
	end do
end subroutine
  
subroutine get_dfdx_direct(x, y, n, dydx)
	implicit none
!	external::spline_mylib, splint_mylib
	integer n,i
	real(8) x(n), y(n), dydx(n)
	real(8) h
	integer isb, ise

	do i=1, n
		isb=1; ise=1
        if(i.eq.1)then
            isb=1;ise=0
        endif
		if(i.eq.n)then
			isb=0; ise=1
		end if
		h=(x(i+isb)-x(i-ise))		
		dydx(i)=(y(i+isb)-y(i-ise))/h
	end do
end subroutine

subroutine get_dfdx_spline(x, y, n, dydx)
	implicit none
!	external::spline_mylib, splint_mylib
	integer n,i
	real(8) x(n), y(n), y2a(n), dydx(n)
	real(8) h
	integer isb, ise
	call spline_mylib(x,y,n,1d30,1d30,y2a)
	do i=1, n
		isb=1; ise=0
		if(i==n)then
			isb=0; ise=1
		end if
		h=(x(i+isb)-x(i-ise))		
		dydx(i)=(y(i+isb)-y(i-ise))/h-0.5d0*(x(i+isb)-x(i-ise))*y2a(i-ise)-1d0/6d0*(x(i+isb)-x(i-ise))*y2a(i+isb)
	end do
end subroutine

! subroutine get_dfdx_spline(x, y, n, dydx)
! 	implicit none
! !	external::spline_mylib, splint_mylib
! 	integer n,idx
! 	real(8) x(n), y(n), ypp(n), dydx(n)
! 	real(8) h
! 	real(8) b,c,d,dx
	
! 	call spline_mylib(x,y,n,1d30,1d30,ypp)
! 	do idx=1, n
! 		h = x(idx+1) - x(idx)
! 		dx = x - x(idx)
! 		! a = y(i)
! 		b = (y(idx+1)-y(idx))/h - h*(2*ypp(idx)+ypp(idx+1))/6
! 		c = ypp(idx)/2
! 		d = (ypp(idx+1)-ypp(idx))/(6*h)
		
! 		! 三次样条导数公式
! 		dydx= b + 2*c*dx + 3*d*dx**2
! 	end do
! end subroutine
subroutine find_interval_bisec(n, x, v, i)
	! 二分查找区间（与原始代码相同）
	integer, intent(in) :: n
	real(8), intent(in) :: x(n)
	real(8), intent(in) :: v
	integer, intent(out) :: i
	integer :: low, high, mid

	if (v < x(1) .or. v > x(n)) then
		i = -1
		return
	end if

	low = 1
	high = n
	do while (low <= high)
		mid = (low + high) / 2
		if (x(mid) < v) then
			low = mid + 1
		else
			high = mid - 1
		end if
	end do

	i = min(max(high, 1), n-1)
end subroutine
subroutine get_dfdx_spline_atx(x, y, ypp, n, x_prob, dydx)
	implicit none
!	external::spline_mylib, splint_mylib
	integer n,idx
	real(8) x(n), y(n), ypp(n), dydx,x_prob
	real(8) h
	real(8) b,c,d,dx
	
	if (x_prob < x(1) .or. x_prob > x(n)) then
		! print*, "x_prob=",x_prob, " not in the range of ", x(1), x(n)
		dydx=-99d0
		return
	end if
        
	! 查找所在区间
	call find_interval_bisec(n, x, x_prob, idx)

	! 计算样条导数
	h = x(idx+1) - x(idx)
	dx = x_prob - x(idx)
	! a = y(i)
	b = (y(idx+1)-y(idx))/h - h*(2*ypp(idx)+ypp(idx+1))/6
	c = ypp(idx)/2
	d = (ypp(idx+1)-ypp(idx))/(6*h)
	
	! 三次样条导数公式
	dydx= b + 2*c*dx + 3*d*dx**2
end subroutine

subroutine get_dfdx_fc(f1,df)
	use com_sts_type
	implicit none
	type(s1d_type)::f1,df
	call get_dfdx_direct(f1%xb,f1%fx,f1%nbin,df%fx)
end subroutine
subroutine get_dfdx_s2d(s2d, flag_dir, s2ddfdx)
    use md_s2_hst_type
	implicit none
    type(s2d_basic_type)::s2d,s2ddfdx
	integer i, flag_dir
	select case(flag_dir)
    case(0)
        do i=1, s2d%ny
            call get_dfdx(s2d%xcenter, s2d%fxy(:,i), s2d%nx, s2ddfdx%fxy(:,i))
        end do
    case(1)
        do i=1, s2d%nx
            call get_dfdx(s2d%ycenter, s2d%fxy(i,:), s2d%ny, s2ddfdx%fxy(i,:))
        end do
    end select
	
end subroutine
subroutine get_dfdx_s2d_direct(s2d, flag_dir, s2ddfdx)
    use md_s2_hst_type
	implicit none
    type(s2d_basic_type)::s2d,s2ddfdx
	integer i, flag_dir
	select case(flag_dir)
    case(0)
        do i=1, s2d%ny
            call get_dfdx_direct(s2d%xcenter, s2d%fxy(:,i), s2d%nx, s2ddfdx%fxy(:,i))
        end do
    case(1)
        do i=1, s2d%nx
            call get_dfdx_direct(s2d%ycenter, s2d%fxy(i,:), s2d%ny, s2ddfdx%fxy(i,:))
        end do
    end select
	
end subroutine