real(8) function rnd(min,max)
		IMPLICIT NONE
		real(8) min,max,temp
		CALL RANDOM_NUMBER (temp)
		if(min<max)then
			rnd=temp*(max-min)+min
		end if
		if(min>max) then
			rnd=temp*(min-max)+max
		end if
end function

real(8) function rnd_external(xmin,xmax, ymin, ymax, frn,par)
    implicit none
    external:: frn
    real(8) xmin,xmax,ymin,ymax
	real(8) par(99),rnd
    real(8) xtmp,ytmp,y

10	xtmp=rnd(xmin,xmax) 
	ytmp=rnd(ymin,ymax)
    call frn(xtmp,y, par)
    if (ytmp>abs(y)) goto 10
 	rnd_external=xtmp
end function
SUBROUTINE set_system_random_seed()
	implicit none
	INTEGER :: i, n, clock
	INTEGER, DIMENSION(:), ALLOCATABLE :: seed

	CALL RANDOM_SEED(size = n)
	ALLOCATE(seed(n))

	CALL SYSTEM_CLOCK(COUNT=clock)

	seed = clock + 37 * (/ (i - 1, i = 1, n) /)
	CALL RANDOM_SEED(PUT = seed)

	DEALLOCATE(seed)
END SUBROUTINE

subroutine same_random_seed(seed_value)
    implicit none
    integer seed_int, seed_value
    integer,allocatable::seed(:)
 
   	call random_seed(size=seed_int)
	allocate(seed(seed_int))
	seed(1:seed_int)=seed_value
    call random_seed(put=seed)
end subroutine
function rndI(min,max)
		IMPLICIT NONE
		Integer min,max, rndI
		real(8) temp
		CALL RANDOM_NUMBER (temp)
		if(min<max)then
			rndI=int(temp*(max-min+1)+min)
		end if
		if(min>max) then
			rndI=int(temp*(min-max+1)+max)
		end if
		if(min.eq.max) rndI=min
end function

real(8) function fPowerLaw_rnd(yta,xmin,xmax)
	implicit none
	real(8) yta,xmin,xmax
	real(8),external::rnd
	if (yta==-1d0) then
		fPowerLaw_rnd=exp(rnd(log(xmin),log(xmax)))
	else
		fPowerLaw_rnd=(rnd(xmin**(1+yta),xmax**(1+yta)))**(1d0/(1d0+yta))
	end if
end function
real(8) function general_rnd(Inverse_F, para)
	implicit none
	real(8),external::Inverse_F
	real(8) tmp
	real(8) para(100)
	real(8),external::rnd
	tmp=rnd(0d0,1d0)
	general_rnd=Inverse_F(tmp, para)	
end function
real(8) function fBPowerLawN_part_rnd(alpha,xb,n,t,c,q,y1,y2)
	implicit none
	integer n, i
	real(8) alpha(n), c(n), xb(n+1), a(n), b(n), s(n+1),t(n+1)
	real(8) ymin, ymax, rnd, Q, tmpx
	real(8) x1,x2,y1,y2

	tmpx=rnd(y1,y2)
	!print*, "tmpx=",tmpx
	do i=1, n
		if(tmpx<t(i+1).and.tmpx>t(i))then
			if(alpha(i).ne.-1d0)then
				!print*, "tmpx, t(i), alpha(i), c(i), xb(i)"
				!print*, tmpx, t(i), alpha(i), c(i), xb(i)
				fBPowerLawN_part_rnd=(Q*(tmpx-t(i))*(1+alpha(i))/c(i)/xb(i)+1)**(1/(alpha(i)+1))*xb(i)
				!print*, "fpowerlawN=", fBpowerlawn_rnd
				!stop
			else
				fBPowerLawN_part_rnd=exp(q*(tmpx-t(i))/c(i)/xb(i))*xb(i)
			end if
		end if
	end do
end function
real(8) function fBPowerLawN_rnd(alpha,xb,n,t,c,q)
	implicit none
	integer n, i
	real(8) alpha(n), c(n), xb(n+1), a(n), b(n), s(n+1),t(n+1)
	real(8) ymin, ymax, rnd, Q, tmpx

	tmpx=rnd(0d0,1d0)
	!print*, "tmpx=",tmpx
	do i=1, n
		if(tmpx<t(i+1).and.tmpx>t(i))then
			if(alpha(i).ne.-1d0)then
				!print*, "tmpx, t(i), alpha(i), c(i), xb(i)"
				!print*, tmpx, t(i), alpha(i), c(i), xb(i)
				fBPowerLawN_rnd=(Q*(tmpx-t(i))*(1+alpha(i))/c(i)/xb(i)+1)**(1/(alpha(i)+1))*xb(i)
				!print*, "fpowerlawN=", fBpowerlawn_rnd
				!stop
			else
				fBPowerLawN_rnd=exp(q*(tmpx-t(i))/c(i)/xb(i))*xb(i)
			end if
		end if
	end do
end function
real(8) function fBPowerLaw_rnd(yta1, yta2, xmid,xmin,xmax)
	implicit none
! breaking power law
	real(8) yta1, yta2, xmid,xmin,xmax
	real(8) c1, c2, c
	real(8),external::rnd
	real(8) tmp,xx
	if(yta1/=-1d0)then
		c1=(xmid**(yta1+1)-xmin**(yta1+1))/(yta1+1)/xmid**yta1
	else
		c1=(log(xmid)-log(xmin))*xmid
	end if
	if(yta2/=-1d0)then
		c2=(xmax**(yta2+1)-xmid**(yta2+1))/(yta2+1)/xmid**yta2
	else
		c2=(log(xmax)-log(xmid))*xmid
	end if
	
	c=c1+c2
	
	tmp=rnd(0d0, 1d0)
	!print*, yta1,yta2
	!stop
	if (tmp<=c1/c) then
		if(yta1/=-1d0)then
			fBPowerLaw_rnd=(tmp*c*(yta1+1)*xmid**yta1+xmin**(yta1+1))**(1d0/(yta1+1))
		else
			fBPowerLaw_rnd=exp(tmp*c/xmid+log(xmin))
		end if
	else
		if(yta2/=-1d0)then
			fBPowerLaw_rnd=((tmp-c1/c)*(yta2+1)*xmid**yta2*c+xmid**(yta2+1))**(1d0/(yta2+1))
		else
			fBPowerLaw_rnd=exp((tmp-c1/c)*c/xmid+log(xmid))
		end if
	end if
	
end function
real(8) function gen_ran_from_data(x,y,n, xmin,xmax,ymax, int_type)
	use constant
	implicit none
	real(8) ranout
	integer n,int_type,ier
	real(8) x(n),y(n), rnd, ymax
	real(8) xmin,xmax, xtmp,ytmp, fy,yp1,ypn
	real(8),allocatable::y2(:)
	allocate(y2(n))
	if(xmin>xmax)pause "error in gen_ran_from_data"
	!print*, "xmin,xmax=", xmin, xmax
	!print*, "ymax=", ymax
100	xtmp=rnd(xmin,xmax)
	ytmp=rnd(0d0,ymax)
	select case(int_type)
	case(1)
		yp1=1d30;ypn=1d30
		call spline_mylib(x,y,n,yp1,ypn,y2)
		call splint_mylib(x,y,y2,n,xtmp,fy,ier)
	case(2)
		call linear_int(x,y,n,xtmp,fy)
	    !print*, "xtmp, fy=",xtmp,fy
	case default
		pause 'error, define int_type'
	end select
	if(ytmp>fy) then
		goto 100
	end if
	gen_ran_from_data=xtmp
    !print*, "gen_ran_from_data=",gen_ran_from_data
!	pause	
end function


integer function discrite_random_selection(plist,n)
	implicit none
	integer n,Isel
	real(8) plist(n)
    integer tmpI
    real(8) tmpp
    real(8),external::rnd
    integer,external::rndI 
100 tmpI=rndI(1,n)
    tmpp=rnd(0d0,1d0)
    if(plist(tmpI)>tmpp)then
        discrite_random_selection=tmpI
		return
    else
        goto 100
    end if
end function

subroutine ran_rotation_sph(vecin,vecout)
   use constant
   implicit none
   real(8) vecin(3),vecout(3)
   real(8) rnd,phi, theta
   phi=rnd(0d0,2*pi)
   theta=acos(rnd(-1d0,1d0))
   call rotation_sph(vecin,vecout, theta,phi)
end subroutine
real(8) function gen_cauchy(a,b)
	use constant
	real(8) a,b
	real(8),external::rnd
	gen_cauchy=b*(tan(pi*(rnd(0d0,1d0)-0.5d0)))+a
end function

real(8) function Rd_Sin(ri)
real(8),external:: rnd
real(8) tem,R,y,ri
integer flag
	tem=rnd(0d0,1d0) 
	y=rnd(0d0,ri)
 do while(tem>abs(sin(y)))
 	tem=rnd(0d0,1d0) 
	y=rnd(0d0,ri)
 end do
Rd_Sin=y
end function

