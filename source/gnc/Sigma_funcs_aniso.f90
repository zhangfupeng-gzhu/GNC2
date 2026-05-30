 
 
subroutine get_gxjcr(gxj, gxjcr, energyx,jcr)
	use constant
	use com_sts_type
	!use 
	implicit none
	integer nbin
	type(s2d_type):: gxj
	real(8) gxjcr
	real(8) energyx, jcr
	real(8) int_out
	integer idid
	!call gxj%print("gxj")
	nbin=gxj%nx
	!print*, "gxj%xmin,xmax=",gxj%xmin,gxj%xmax
	!print*, "xc=",gxj%xcenter
	!print*, "nbin,type=",nbin, gxj%sts_type
	!print*, "energyx,jcr=",energyx,jcr
	int_out=0d0
	call my_integral_none(0d0,pi/2d0,int_out, fcn_y, idid)
	!print*, "out=",int_out
	!read(*,*)
	gxjcr=int_out
contains
    subroutine fcn_y(n, x, y, f, par, ipar)
		use, intrinsic :: ieee_arithmetic
		implicit none
		integer n, ipar(100)
		real(8) x, y(n), f(n), par(100)
		integer idx, idy
		real(8) evl
		evl=max(jcr+log10(sin(x)),gxj%ymin)
		!print*, "fcn:energyx,evl=", energyx,evl
		!print*, "gxj%xmin, xmax, ymin,ymax=", gxj%xmin,gxj%xmax,gxj%ymin,gxj%ymax
		!read(*,*)
		call return_idxy(energyx,evl,gxj%xmin,gxj%xmax,gxj%ymin,gxj%ymax,&
		nbin,nbin,idx,idy,sts_type_dstr)
		if(idy>gxj%nx.or.idy<1)then
			print*, "energyx, evl,jcr=",energyx, evl, jcr
			print*, "idx,idy=",idx,idy
			stop
		end if
        f(1)=sin(x)*gxj%fxy(idx,idy)
		!print*, "f(1)=",f(1)
	end subroutine
end subroutine
 