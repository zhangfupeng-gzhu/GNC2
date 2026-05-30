!     IDID        REPORTS ON SUCCESSFULNESS UPON RETURN:
!                   IDID= 1  COMPUTATION SUCCESSFUL,
!                   IDID= 2  COMPUT. SUCCESSFUL (INTERRUPTED BY SOLOUT)
!                   IDID=-1  INPUT IS NOT CONSISTENT,
!                   IDID=-2  LARGER NMAX IS NEEDED,
!                   IDID=-3  STEP SIZE BECOMES TOO SMALL.
!                   IDID=-4  PROBLEM IS PROBABLY STIFF (INTERRUPTED).

module my_intgl
	implicit none
    !interface my_integral
    !    module procedure my_integral_multi
    !    module procedure my_integral_multi2
    !    module procedure my_integral_none
    !    module procedure my_integral_pars
    !    module procedure my_integral_acc
    !    module procedure my_integral_interpolate
    !    module procedure my_integral_interpolate_par
    !end interface
contains
    subroutine my_integral_all(xs,xe,y,n_y,FCN,rtol_in, atol_in, rpar_in, ipar_in,order_in)
        implicit none
        integer LWORK,ITOL,LIWORK,IDID
		integer n_y,ipar(100)
		real(8) y(n_y)
        integer,optional:: ipar_in(100),order_in
		real(8),optional:: RTOL_in(n_y),ATOL_in(n_y),rpar_in(100)
        real(8):: RTOL(n_y),ATOL(n_y),rpar(100)
		real(8)	WORK(n_y*100)
		real(8) x_start,x_end,xs,xe
		integer:: IOUT,order=5
		integer IWORK(n_y*20)
		external::dopri5,FCN,my_solout_empty

        if(present(rtol_in)) rtol=rtol_in
        if(present(atol_in)) atol=atol_in
        if(present(rpar_in)) rpar=rpar_in
        if(present(ipar_in)) ipar=ipar_in
        if(present(order_in)) order=order_in
        ITOL=1;LWORK=n_y*100;WORK=0;IWORK=0;LIWORK=n_y*20
		IWORK(1)=10000000;
		WORK(7)=0.
		IOUT=0

		x_start=xs
		x_end=xe
		select case(order)
		case(5)
			call dopri5(n_y,FCN,x_start,y,x_end,RTOL,ATOL,ITOL,my_solout_empty,&
				IOUT,WORK,LWORK,IWORK,LIWORK,RPAR,ipar,IDID)
		case(8)
			call DOP853(n_y,FCN,x_start,y,x_end,RTOL,ATOL,ITOL,my_solout_empty,&
				IOUT,WORK,LWORK,IWORK,LIWORK,RPAR,ipar,IDID)
		endselect
    end subroutine
	subroutine my_integral_multi(xs, xe, y, n_y,order,FCN, rpar, ipar)
	!!  subroutine FCN(N,X,Y,F,RPAR,IPAR)
		IMPLICIT NONE
		integer LWORK,ITOL,LIWORK,IDID
		integer n_y,ipar(100)
		real(8) y(n_y)
		real(8) RTOL(n_y),ATOL(n_y),rpar(100)
		real(8)	WORK(n_y*100)
		real(8) x_start,x_end,xs,xe
		integer IOUT,order
		integer IWORK(n_y*20)
		external::dopri5,FCN,my_solout_empty

		if(rpar(100).eq.0d0)then
			RTOL=1d-12
			ATOL=1d-12
		else
			RTOL=rpar(100)
			ATOL=rpar(99)
            if(rpar(100)<1d-20.or.rpar(99)<1d-20)then
                print*, "warnning: rpar(99:100) not initilized?"
            end if
		end if
		ITOL=1;LWORK=n_y*100;WORK=0;IWORK=0;LIWORK=n_y*20
		IWORK(1)=10000000;
		WORK(7)=0.
		IOUT=0

		x_start=xs
		x_end=xe
		select case(order)
		case(5)
			call dopri5(n_y,FCN,x_start,y,x_end,RTOL,ATOL,ITOL,my_solout_empty,&
				IOUT,WORK,LWORK,IWORK,LIWORK,RPAR,ipar,IDID)
		case(8)
			call DOP853(n_y,FCN,x_start,y,x_end,RTOL,ATOL,ITOL,my_solout_empty,&
				IOUT,WORK,LWORK,IWORK,LIWORK,RPAR,ipar,IDID)
		endselect
	end subroutine

	subroutine my_integral_multi2(xs, xe, y, n_y,FCN, rpar_in, ipar_in, hrlow_in,hrhigh_in)
	!!  subroutine FCN(N,X,Y,F,RPAR,IPAR)
		IMPLICIT NONE
		integer LWORK,ITOL,LIWORK,IDID
		integer n_y,ipar(100)
		real(8) y(n_y)
		real(8) RTOL(n_y),ATOL(n_y),rpar(100)
		real(8)	WORK(400)
		real(8) x_start,x_end,xs,xe
		integer IOUT
		integer IWORK(150)
		real(8),optional,intent(in):: hrlow_in, hrhigh_in
		real(8),optional,dimension(:):: rpar_in
		integer,optional,dimension(:):: ipar_in
		real(8) hrlow,hrhigh
		external::dopri5,FCN,my_solout_empty
		if(rpar(100).eq.0d0)then
			RTOL=1d-12
			ATOL=1d-12
		else
			RTOL=rpar(100)
			ATOL=rpar(99)
		end if
		ITOL=1;LWORK=400;WORK=0;IWORK=0;LIWORK=400
		WORK(1)=1d-30;  
		IWORK(3)=6;
		if(present(hrlow_in))then
			hrlow=hrlow_in
			hrhigh=hrhigh_in
			WORK(3)=hrlow;WORK(4)=hrhigh
		end if
		!print*,present(rpar_in)
		if(present(rpar_in))then
			rpar=rpar_in
			!print*,rpar_in
		end if
		if(present(ipar_in))then
			ipar=ipar_in
		end if
		IWORK(1)=10000000;
		WORK(7)=0.
		IOUT=0

		x_start=xs
		x_end=xe
		call dopri5(n_y,FCN,x_start,y,x_end,RTOL,ATOL,ITOL,my_solout_empty,IOUT,&
			WORK,LWORK,IWORK,LIWORK,RPAR,ipar,IDID)
	
	end subroutine

	subroutine my_integral_none(xs, xe, y, FCN,idid)
	!!  subroutine FCN(N,X,Y,F,RPAR,IPAR)
		IMPLICIT NONE
		integer ITOL,IDID
		integer n_y,ipar(100)
		real(8) y,yout(1)
		real(8) RTOL(1),ATOL(1),rpar(100)
		integer,parameter::LWORK=400,LIWORK=400
		real(8)	WORK(LWORK)
		integer IWORK(LIWORK)
		real(8) x_start,x_end,xs,xe
		integer IOUT
		real(8) hrlow,hrhigh
		external::dopri5,FCN,my_solout_empty

		RTOL=1d-12
		ATOL=1d-12
		ITOL=0;;WORK=0;IWORK=0
		!WORK(1)=1d-30;  
		!IWORK(3)=6;   ! print messages
		IWORK(1)=10000000;
		WORK(7)=0.
		IOUT=0
		x_start=xs
		x_end=xe
        idid=0
		if(x_start.eq.x_end) then
			y=0
			return
		end if
		yout(1)=y
		call dopri5(1,FCN,x_start,yout,x_end,RTOL,ATOL,ITOL,my_solout_empty,IOUT,&
			WORK,LWORK,IWORK,LIWORK,RPAR,ipar,IDID)
		y=yout(1)
	end subroutine
	subroutine my_integral_pars(xs, xe, y,ATOL, RTOL, RPAR, IPAR, FCN,idid)
	!!  subroutine FCN(N,X,Y,F,RPAR,IPAR)
		IMPLICIT NONE
		integer ITOL,IDID
		integer n_y,ipar(100)
		real(8) y,yout(1)
		real(8) RTOL(1),ATOL(1),rpar(100)
		integer,parameter::LWORK=400,LIWORK=400
		real(8)	WORK(LWORK)
		integer IWORK(LIWORK)
		real(8) x_start,x_end,xs,xe
		integer IOUT
		real(8) hrlow,hrhigh
		external::dopri5,FCN,my_solout_empty

		ITOL=0;WORK=0;IWORK=0
		!WORK(1)=1d-30;  
		!IWORK(3)=6;   ! print messages
		IWORK(1)=10000000;
		WORK(7)=0.
		IOUT=0

		x_start=xs
		x_end=xe
		yout(1)=y
		call dopri5(1,FCN,x_start,yout,x_end,RTOL,ATOL,ITOL,my_solout_empty,IOUT,WORK,LWORK,IWORK,LIWORK,RPAR,ipar,IDID)
		y=yout(1)
	end subroutine
	subroutine my_integral_acc(xs, xe, y,ATOL_in, RTOL_in, FCN,IDID)
	!!  subroutine FCN(N,X,Y,F,RPAR,IPAR)
		IMPLICIT NONE
		integer LWORK,ITOL,LIWORK,IDID
		integer n_y,ipar(100)
		real(8) y,yout(1)
		real(8) RTOL(1),ATOL(1),rpar(100),atol_in,rtol_in
		real(8)	WORK(400)
		real(8) x_start,x_end,xs,xe
		integer IOUT
		integer IWORK(150)
		real(8) hrlow,hrhigh
		external::dopri5,FCN,my_solout_empty

		ITOL=0;LWORK=400;WORK=0;IWORK=0;LIWORK=400
		!WORK(1)=1d-30;  
		!IWORK(3)=6;   ! print messages
		IWORK(1)=10000000;
		WORK(7)=0.
		IOUT=0

		x_start=xs
		x_end=xe
		rtol=rtol_in
		atol=atol_in
		yout(1)=y
		idid=0
		if(x_start.eq.x_end) then
			y=0
			return
		end if
		call dopri5(1,FCN,x_start,yout,x_end,RTOL,ATOL,ITOL,my_solout_empty,IOUT,&
			WORK,LWORK,IWORK,LIWORK,RPAR,ipar,IDID)
		y=yout(1)
	end subroutine
	subroutine my_integral_acc_solout(xs, xe, y,ATOL_in, RTOL_in, FCN,MY_SO_OUT,IDID)
		!!  subroutine FCN(N,X,Y,F,RPAR,IPAR)
			IMPLICIT NONE
			integer LWORK,ITOL,LIWORK,IDID
			integer n_y,ipar(100)
			real(8) y,yout(1)
			real(8) RTOL(1),ATOL(1),rpar(100),atol_in,rtol_in
			real(8)	WORK(400)
			real(8) x_start,x_end,xs,xe
			integer IOUT
			integer IWORK(150)
			real(8) hrlow,hrhigh
			external::dopri5,FCN,MY_SO_OUT
	
			ITOL=0;LWORK=400;WORK=0;IWORK=0;LIWORK=400
			!WORK(1)=1d-30;  
			!IWORK(3)=6;   ! print messages
			IWORK(1)=10000000;
			WORK(7)=0.
			IOUT=1
	
			x_start=xs
			x_end=xe
			rtol=rtol_in
			atol=atol_in
			yout(1)=y
			idid=0
			if(x_start.eq.x_end) then
				y=0
				return
			end if
			call dopri5(1,FCN,x_start,yout,x_end,RTOL,ATOL,ITOL,MY_SO_OUT,IOUT,&
				WORK,LWORK,IWORK,LIWORK,RPAR,ipar,IDID)
			y=yout(1)
		end subroutine
		subroutine my_integral_acc_solout_multi(xs, xe, y,n_y,ATOL_in, RTOL_in, FCN,MY_SO_OUT,IDID)
			!!  subroutine FCN(N,X,Y,F,RPAR,IPAR)
				IMPLICIT NONE
				integer LWORK,ITOL,LIWORK,IDID
				integer n_y,ipar(100)
				real(8) y(n_y),yout(n_y)
				real(8) RTOL(n_y),ATOL(n_y),rpar(100),atol_in,rtol_in
				real(8)	WORK(400)
				real(8) x_start,x_end,xs,xe
				integer IOUT
				integer IWORK(150)
				real(8) hrlow,hrhigh
				external::dopri5,FCN,MY_SO_OUT
		
				ITOL=0;LWORK=400;WORK=0;IWORK=0;LIWORK=400
				!WORK(1)=1d-30;  
				!IWORK(3)=6;   ! print messages
				IWORK(1)=10000000;
				WORK(7)=0.
				IOUT=1
		
				x_start=xs
				x_end=xe
				rtol=rtol_in
				atol=atol_in
				yout=y
				idid=0
				if(x_start.eq.x_end) then
					y=0
					return
				end if
				call dopri5(n_y,FCN,x_start,yout,x_end,RTOL,ATOL,ITOL,MY_SO_OUT,IOUT,&
					WORK,LWORK,IWORK,LIWORK,RPAR,ipar,IDID)
				y=yout
			end subroutine
			subroutine my_integral_acc_multi(xs, xe, y,n_y,ATOL_in, RTOL_in, FCN,IDID,PAR)
				!!  subroutine FCN(N,X,Y,F,RPAR,IPAR)
					IMPLICIT NONE
					integer LWORK,ITOL,LIWORK,IDID
					integer n_y,ipar(100)
					real(8) y(n_y),yout(n_y)
					real(8) RTOL(n_y),ATOL(n_y),rpar(100),atol_in,rtol_in
					real(8)	WORK(400)
					real(8),optional::par
					real(8) x_start,x_end,xs,xe
					integer IOUT
					integer IWORK(150)
					real(8) hrlow,hrhigh
					external::dopri5,FCN,my_solout_empty
			
					ITOL=0;LWORK=400;WORK=0;IWORK=0;LIWORK=400
					!WORK(1)=1d-30;  
					!IWORK(3)=6;   ! print messages
					IWORK(1)=10000000;
					WORK(7)=0.
					IOUT=0
			
					x_start=xs
					x_end=xe
					rtol=rtol_in
					atol=atol_in
					yout=y
					idid=0
					if(x_start.eq.x_end) then
						y=0
						return
					end if
					call dopri5(n_y,FCN,x_start,yout,x_end,RTOL,ATOL,ITOL,my_solout_empty,IOUT,&
						WORK,LWORK,IWORK,LIWORK,RPAR,ipar,IDID)
					y=yout
					if(present(par))then
						if(idid<0)then
							par=x_start
						end if
					endif
				end subroutine
	subroutine my_integral2(xs, xe, y, FCN, rpar_in, ipar_in, hrlow_in,hrhigh_in)
	!!  subroutine FCN(N,X,Y,F,RPAR,IPAR)
		IMPLICIT NONE
		integer LWORK,ITOL,LIWORK,IDID
		integer n_y,ipar(100)
		real(8) y,yout(1)
		real(8) RTOL(1),ATOL(1),rpar(100)
		real(8)	WORK(400)
		real(8) x_start,x_end,xs,xe
		integer IOUT
		integer IWORK(150)
		real(8),optional,intent(in):: hrlow_in, hrhigh_in
		real(8),optional,dimension(:):: rpar_in
		integer,optional,dimension(:):: ipar_in
		real(8) hrlow,hrhigh
		external::dopri5,FCN,my_solout_empty
		RTOL=1d-12
		ATOL=1d-12
		ITOL=1;LWORK=400;WORK=0;IWORK=0;LIWORK=400
		!WORK(1)=1d-30;  
		IWORK(3)=6;
		if(present(hrlow_in))then
			hrlow=hrlow_in
			hrhigh=hrhigh_in
			WORK(3)=hrlow;WORK(4)=hrhigh
		end if
		!print*,present(rpar_in)
		if(present(rpar_in))then
			rpar=rpar_in
			!print*,rpar_in
		end if
		if(present(ipar_in))then
			ipar=ipar_in
		end if
		IWORK(1)=10000000;
		WORK(7)=0.
		IOUT=0

		x_start=xs
		x_end=xe
		yout(1)=y
		call dopri5(1,FCN,x_start,yout,x_end,RTOL,ATOL,ITOL,my_solout_empty,IOUT,WORK,LWORK,IWORK,LIWORK,RPAR,ipar,IDID)
		y=yout(1)
	end subroutine
	
	subroutine my_integral_interpolate(xs, xe, xin, yin,ar_n, y)
		IMPLICIT NONE
		integer LWORK,ITOL,LIWORK,IDID
		integer ipar(100),ar_n
		real(8) yin(ar_n),xin(ar_n),yout(1),y
		real(8) RTOL(1),ATOL(1),rpar(100)
		real(8)	WORK(400)
		real(8) x_start,x_end,xs,xe
		integer IOUT,nout
		integer IWORK(150)
		real(8) hrlow,hrhigh
		real(8) yp1,ypn
		real(8),allocatable::y2(:)		
		external::dopri5,my_solout_empty,spline_mylib,splint_mylib
		RTOL=1d-12
		ATOL=1d-12
		ITOL=0;LWORK=400;WORK=0;IWORK=0;LIWORK=400
		!WORK(1)=1d-30;  
		!IWORK(3)=6;
		IWORK(1)=10000000;
		WORK(7)=0.
		IOUT=0
	
		x_start=xs
		x_end=xe

		yp1=1d33;ypn=1d33
		allocate(y2(ar_n))
		call spline_mylib(xin(1:ar_n),yin(1:ar_n),ar_n,yp1,ypn,y2(1:ar_n))	
		yout(1)=y
		call dopri5(1,my_FCN_interpolate,x_start,yout,x_end,RTOL,ATOL,ITOL,my_solout_empty,IOUT,WORK,LWORK,IWORK,LIWORK,RPAR,ipar,IDID)
		y=yout(1)
		contains
			subroutine my_FCN_interpolate(N,X,Y,F,RPAR,IPAR)
				implicit none
				integer N,ier
				real(8) X,Y(N),F(N)
				integer IPAR(100)
				real(8) RPAR(100)	
				real(8) ytmp
				call splint_mylib(xin(1:ar_n),yin(1:ar_n),y2(1:ar_n),ar_n,X,ytmp,ier)
				F(1)=ytmp
			end subroutine
	end subroutine

	subroutine my_integral_interpolate_par(xs, xe, xin, yin,ar_n,ATOL,RTOL,ylog, y)
		IMPLICIT NONE
		integer LWORK,ITOL,LIWORK,IDID
		integer ipar(100),ar_n
		real(8) yin(ar_n),xin(ar_n),y,yout(1)
		real(8) RTOL(1),ATOL(1),rpar(100)
		real(8)	WORK(400)
		real(8) x_start,x_end,xs,xe
		logical ylog
		integer IOUT,nout
		integer IWORK(150)
		real(8) hrlow,hrhigh
		real(8) yp1,ypn
		real(8),allocatable::y2(:)		
		external::dopri5,my_solout_empty,spline_mylib,splint_mylib
		ITOL=0;LWORK=400;WORK=0;IWORK=0;LIWORK=400
		!WORK(1)=1d-30;  
		!IWORK(3)=6;
		IWORK(1)=10000000;
		WORK(7)=0.
		IOUT=0
	
		x_start=xs
		x_end=xe

		yp1=1d33;ypn=1d33
		allocate(y2(ar_n))
		call spline_mylib(xin(1:ar_n),yin(1:ar_n),ar_n,yp1,ypn,y2(1:ar_n))	
		yout(1)=y
		call dopri5(1,my_FCN_interpolate,x_start,yout,x_end,RTOL,ATOL,ITOL,my_solout_empty,IOUT,WORK,LWORK,IWORK,LIWORK,RPAR,ipar,IDID)
		y=yout(1)
		contains
			subroutine my_FCN_interpolate(N,X,Y,F,RPAR,IPAR)
				implicit none
				integer N,ier
				real(8) X,Y(N),F(N)
				integer IPAR(100)
				real(8) RPAR(100)	
				real(8) ytmp
				call splint_mylib(xin(1:ar_n),yin(1:ar_n),y2(1:ar_n),ar_n,X,ytmp,ier)
				if(ylog)then
					F(1)=10**ytmp
				else
					F(1)=ytmp
				end if
			end subroutine
	end subroutine

	subroutine my_integral_interpolate_acc(xs, xe, xin, yin,ar_n,ATOL,RTOL, y)
		IMPLICIT NONE
		integer LWORK,ITOL,LIWORK,IDID
		integer ipar(100),ar_n
		real(8) yin(ar_n),xin(ar_n),yout(1),y
		real(8) RTOL(1),ATOL(1),rpar(100)
		real(8)	WORK(400)
		real(8) x_start,x_end,xs,xe
		integer IOUT,nout
		integer IWORK(150)
		real(8) hrlow,hrhigh
		real(8) yp1,ypn
		real(8),allocatable::y2(:)		
		external::dopri5,my_solout_empty,spline_mylib,splint_mylib
		ITOL=0;LWORK=400;WORK=0;IWORK=0;LIWORK=400
		!WORK(1)=1d-30;  
		!IWORK(3)=6;
		IWORK(1)=10000000;
		WORK(7)=0.
		IOUT=0
	
		x_start=xs
		x_end=xe

		yp1=1d33;ypn=1d33
		allocate(y2(ar_n))
		call spline_mylib(xin(1:ar_n),yin(1:ar_n),ar_n,yp1,ypn,y2(1:ar_n))	
		yout(1)=y
		call dopri5(1,my_FCN_interpolate,x_start,yout,x_end,RTOL,ATOL,ITOL,my_solout_empty,IOUT,WORK,LWORK,IWORK,LIWORK,RPAR,ipar,IDID)
		y=yout(1)
		contains
			subroutine my_FCN_interpolate(N,X,Y,F,RPAR,IPAR)
				implicit none
				integer N,ier
				real(8) X,Y(N),F(N)
				integer IPAR(100)
				real(8) RPAR(100)	
				real(8) ytmp
				call splint_mylib(xin(1:ar_n),yin(1:ar_n),y2(1:ar_n),ar_n,X,ytmp,ier)
				F(1)=ytmp
			end subroutine
	end subroutine

	subroutine my_integral_interpolate_ylog(xs, xe, xin, yin,ar_n, y)
		IMPLICIT NONE
		integer LWORK,ITOL,LIWORK,IDID
		integer ipar(100),ar_n
		real(8) yin(ar_n),xin(ar_n),yout(1),y
		real(8) RTOL(1),ATOL(1),rpar(100)
		real(8)	WORK(400)
		real(8) x_start,x_end,xs,xe
		integer IOUT,nout
		integer IWORK(150)
		real(8) hrlow,hrhigh
		real(8) yp1,ypn
		real(8),allocatable::y2(:)		
		external::dopri5,my_solout_empty,spline_mylib,splint_mylib
		RTOL=1d-12
		ATOL=1d-12
		ITOL=0;LWORK=400;WORK=0;IWORK=0;LIWORK=400
		!WORK(1)=1d-30;  
		!IWORK(3)=6;
		IWORK(1)=10000000;
		WORK(7)=0.
		IOUT=0
	
		x_start=xs
		x_end=xe

		yp1=1d33;ypn=1d33
		allocate(y2(ar_n))
		call spline_mylib(xin(1:ar_n),yin(1:ar_n),ar_n,yp1,ypn,y2(1:ar_n))	
		yout(1)=y
		call dopri5(1,my_FCN_interpolate,x_start,yout,x_end,RTOL,ATOL,ITOL,my_solout_empty,IOUT,WORK,LWORK,IWORK,LIWORK,RPAR,ipar,IDID)
		y=yout(1)
		contains
			subroutine my_FCN_interpolate(N,X,Y,F,RPAR,IPAR)
				implicit none
				integer N
				real(8) X,Y(N),F(N)
				integer IPAR(100),ier
				real(8) RPAR(100)	
				real(8) ytmp
				call splint_mylib(xin(1:ar_n),yin(1:ar_n),y2(1:ar_n),ar_n,X,ytmp,ier)
				F(1)=10**ytmp
			end subroutine
	end subroutine
	subroutine remove_zeros(x,y,n, xmd,ymd, nout)
		implicit none
		integer n,i,nout,j
		real(8) x(n),y(n),xmd(n),ymd(n)
		j=0
		do i=1, n
			if(y(i).ne.0d0)then
				j=j+1
				xmd(j)=x(i); ymd(j)=y(i)
			end if
		end do
		nout=j
	end subroutine
	subroutine my_integral_linear(xs, xe, x_in, y_in,ar_n, y,xlog_in,ylog_in)
		IMPLICIT NONE
		integer LWORK,ITOL,LIWORK,IDID
		integer ipar(100),ar_n
		real(8) y_in(ar_n),x_in(ar_n),yout(1),y
		real(8) RTOL(1),ATOL(1),rpar(100)
		real(8)	WORK(400)
		real(8) x_start,x_end,xs,xe
		integer IOUT,nout
		integer IWORK(150)
		real(8) hrlow,hrhigh
		real(8) yp1,ypn
		real(8),allocatable::xmd(:),ymd(:),y2md(:)		
		integer,optional,intent(in):: xlog_in, ylog_in		
		integer::xlog,ylog
		external::dopri5,my_solout_empty,spline_mylib,splint_mylib
		RTOL=1d-12
		ATOL=1d-12
		ITOL=0;LWORK=400;WORK=0;IWORK=0;LIWORK=400
		!WORK(1)=1d-30;  
		!IWORK(3)=6;
		IWORK(1)=10000000;
		WORK(7)=0.
		IOUT=0
	
		x_start=xs
		x_end=xe
		xlog=0;ylog=0
		if(present(xlog_in))then
			xlog=xlog_in
		end if
		if(present(ylog_in))then
			ylog=ylog_in
		end if	
	
		allocate(xmd(ar_n));allocate(ymd(ar_n));allocate(y2md(ar_n))

		yp1=1d30;ypn=1d30

!		call remove_zero(x,y, ar_n,xmd,ymd,nout)
!		if(x_start<xmd(1)) x_start=xmd(1)
!		if(x_end>xmd(nout)) x_end=xmd(nout)

		xmd(1:ar_n)=x_in(1:ar_n); ymd(1:ar_n)=y_in(1:ar_n); nout=ar_n
		yout(1)=y
		call dopri5(1,my_FCN_linear,x_start,yout,x_end,RTOL,ATOL,ITOL,my_solout_empty,IOUT,WORK,LWORK,IWORK,LIWORK,RPAR,ipar,IDID)
		y=yout(1)
		deallocate(xmd, ymd,y2md)
		contains
			subroutine my_FCN_linear(N,X,Y,F,RPAR,IPAR)
				implicit none
				integer N
				real(8) X,Y(N),F(N)
				real(8) RPAR(100),IPAR(100)	
				real(8) ytmp
				if(xlog.eq.0)then
					call linear_int(xmd(1:nout),ymd(1:nout),nout,X,ytmp)
				else
					call linear_int(xmd(1:nout),ymd(1:nout),nout,log10(X),ytmp)
				end if
				if(ylog.eq.0)then
					F(1)=ytmp
				else
					F(1)=10**ytmp
				end if				
			!	print*,X, F(1), y
			end subroutine
	end subroutine
	
	
end module
subroutine my_gs_integral(xbg,xed, res, alph,beta,epsabs,epsrel,integr,fx,ier)
	IMPLICIT NONE
!   ier is now used as an input flag, if ier=-99, then the error screen output is suppressed
!            INTEGR - Integer
!                     Indicates which WEIGHT function is to be used
!                     = 1  (X-A)**ALFA*(B-X)**BETA
!                     = 2  (X-A)**ALFA*(B-X)**BETA*LOG(X-A)
!                     = 3  (X-A)**ALFA*(B-X)**BETA*LOG(B-X)
!                     = 4  (X-A)**ALFA*(B-X)**BETA*LOG(X-A)*LOG(B-X)
!					  A, B are the lower and the upper limit of integration.
!                     If INTEGR.LT.1 or INTEGR.GT.4, the routine
!                     will end with IER = 6.
!		  			  Alph>-1, beta>-1
!                     fx in forms of fx(x)
	real(8),external::fx
	real(8) xbg,xed,alph,beta
	real(8) epsabs,epsrel ! absolute and relative accuracy
	integer INTEGR,NEVAL,IER,LAST
	integer,parameter::LIMIT=150, LENW=LIMIT*4
	integer IWORK(LIMIT)
	real(8) WORK(LENW),abserr 
	real(8) res      !the result of the integration
!	epsabs=1d-12;epsrel=1d-12;
	logical suppress_screen_error
	if(ier.eq.-99) then
		suppress_screen_error=.true.
	else
		suppress_screen_error=.false.
	end if
	call dqaws(fx,xbg,xed,alph,beta,INTEGR,epsabs,epsrel,res,abserr,neval,ier,limit,lenw,last,iwork,work,suppress_screen_error)

end subroutine

subroutine my_solout_empty(NR,XOLD,X,Y,N,CON,ICOMP,ND,RPAR,IPAR,IRTRN)
!     SOLOUT      NAME (EXTERNAL) OF SUBROUTINE PROVIDING THE
!                 NUMERICAL SOLUTION DURING INTEGRATION. 
!                 IF IOUT.GE.1, IT IS CALLED AFTER EVERY SUCCESSFUL STEP.
!                 SUPPLY A DUMMY SUBROUTINE IF IOUT=0. 
!                 IT MUST HAVE THE FORM
!                    SUBROUTINE SOLOUT (NR,XOLD,X,Y,N,CON,ICOMP,ND,
!                                       RPAR,IPAR,IRTRN)
!                    DIMENSION Y(N),CON(5*ND),ICOMP(ND)
!                   ....  
!                 SOLOUT FURNISHES THE SOLUTION "Y" AT THE NR-TH
!                    GRID-POINT "X" (THEREBY THE INITIAL VALUE IS
!                    THE FIRST GRID-POINT).
!                 "XOLD" IS THE PRECEEDING GRID-POINT.
!                 "IRTRN" SERVES TO INTERRUPT THE INTEGRATION. IF IRTRN
!                    IS SET <0, DOPRI5 WILL RETURN TO THE CALLING PROGRAM.
!                    IF THE NUMERICAL SOLUTION IS ALTERED IN SOLOUT,
!                    SET  IRTRN = 2
!           
!          -----  CONTINUOUS OUTPUT: -----
!                 DURING CALLS TO "SOLOUT", A CONTINUOUS SOLUTION
!                FOR THE INTERVAL [XOLD,X] IS AVAILABLE THROUGH
!                 THE FUNCTION
!                        >>>   CONTD5(I,S,CON,ICOMP,ND)   <<<
!                 WHICH PROVIDES AN APPROXIMATION TO THE I-TH
!                 COMPONENT OF THE SOLUTION AT THE POINT S. THE VALUE
!                 S SHOULD LIE IN THE INTERVAL [XOLD,X].
	implicit none
	INTEGER N, ND, NR
	real(8) X, XOLD
    real(8) Y(N),CON(5*ND),ICOMP(ND)
	integer IRTRN,IPAR(100)
	real(8) rpar(100)
	
end subroutine
!subroutine FCN(N,X,Y,F,RPAR,IPAR)
!integer N
!integer IPAR(100)
!real(8) rpar(100)
!real(8) X,Y(N),F(N)
!end subroutine
