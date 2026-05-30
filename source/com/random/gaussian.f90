module md_gaussian
	implicit none
	integer::flag=0
	real(8)::v1,v2,s
!!$OMP threadprivate(flag, v1,v2,s)	
end module
real(8) function gen_gaussian(sigma)
	use md_gaussian
! generate an Gaussian variable with variation=sigma, centered at zero
! test show that it take 4.57s to complete 1d8 loop
	IMPLICIT NONE
	real(8) u1,u2,x,sigma
!	data flag/0/,s/0./
	real(8),external::rnd

	if(flag==0)then
		s=0.
		do while(s>=1.or.s==0.)
            call random_number(u1)
            call random_number(u2)
			!u1=rnd(0d0,1d0)
			!u2=rnd(0d0,1d0)
			v1=2*u1-1
			v2=2*u2-1
			s=v1*v1+v2*v2
		end do
			x=v1*sqrt(-2*log(s)/s)
	!		print*, "x1=",x
        flag=1
	else
		x=v2*sqrt(-2*log(s)/s)
        flag=0
	!	print*, "x2=",x
	end if
	!flag=1-flag
	x=x*sigma
	gen_gaussian=x
	!print*, "x3=",x, sigma
end function

real(8) function gen_normal(sigma)
! generate an Gaussian variable with variation=sigma, centered at zero
! test show that it take 4.88s to complete 1d8 loop 
IMPLICIT NONE
    double precision, parameter :: pi = 3.141592653589793239
	integer,save::flag=0
	real(8),save::v1,v2,s
	real(8) u1,u2,x,sigma

    call random_number(u1)
    call random_number(u2)
	if(flag.eq.0)then
        s=sqrt(-2d0*log(u1))*cos(2d0*pi*u2)
        gen_normal=sigma*s
        flag=1
	else
        s=sqrt(-2d0*log(u1))*sin(2d0*pi*u2)
        gen_normal=sigma*s
        flag=0
	end if

end function

subroutine gen_gaussian_correlate(y1,y2, coeff)
	!Generate two random variables y1, and y2, with mean value <y1>=<y2>=0, 
	!dispersion <y1^2>=<y2^2>=1 and cross-correlation <y1y2>=coeff,where coeff<1
	implicit none
	real(8) y1p,y2p,y1,y2,coeff, gen_gaussian
	y1p=gen_gaussian(1d0); y2p=gen_gaussian(1d0)
	y1=y1p
	if(abs(coeff).le.1)then
		y2=y1p*coeff+y2p*sqrt(1-coeff**2)	
	else
		y2=coeff/abs(coeff)*y1
	end if
end subroutine

!subroutine gen_binormal(y1,y2,coeff)
!    implicit none
!    real(8) y1p,y2p, gen_gaussian, coeff
!    real(8) y1, y2, um, rnd
!
!    y1p=gen_gaussian(1d0)
!    y2p=gen_gaussian(1d0)
!    um=rnd(0d0,1d0)
!    print*, um, coeff
!    if(um>coeff)then
!        y1=y1p; y2=y2p
!        return
!    end if
!    if(um<=coeff.and.y1p+y2p>0)then
!        y1=1d0;y2=1d0
!       ! y1=y1p;
!       ! y2=y1
!        return
!    end if
!    if(um<=coeff.and.y1p+y2p<=0)then
!        y1=-1d0;y2=-1d0
!    !    y1=y1p; 
!    !    y2=y1
!    end if
!    
!end subroutine
