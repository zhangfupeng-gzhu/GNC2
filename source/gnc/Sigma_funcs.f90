  
subroutine get_sigma0(energyx,  funcs, sigma)
	use com_sts_type
	use constant
	use md_coeff
    use model_basic,only: fgx_g0,sample_emin,sample_emax,fc_ir_share
	implicit none
	integer i,idid
	real(8) energyx, sigma, ecc,vel,rad,semi, int_out,xmin,xmax
    !real(8) intne, intpo
	external::funcs
!	integer, parameter::nintn=25
	int_out=0
	!print*, "inf=",inf
    idid=0
	!if(energyx.ne.0d0)then
        !intne=0d0; intpo=0d0
    !    print*, emin_factor, energyx
    xmin=10**fc_ir_share%xmin
    xmax=10**fc_ir_share%xmax
    if(energyx/emin_factor>1d0+1d-7)then
		call my_integral_none(emin_factor/energyx, 1d0, int_out, fcn,idid)
        if(idid.lt.0)then
            print*, "error! idid=",idid, energyx, emin_factor
            stop
        end if
    end if
    !read(*,*)
        !if(emin_factor<energyx)then
        !    call my_integral(emin_factor/energyx, 1d0, intpo, fcn,idid)
        !end if
        !int_out=intne+intpo
	!else
!		call my_integral_none(inf, 1d0, int_out, fcn,idid)
!        stop
!	end if
	sigma=int_out!+fgx_g0/energyx
    !if(energyx>10)then
    !    call fc_ir_share%print("fc_share")
	    !print*, "sigma=", energyx, sigma,sample_emin,sample_emax
	    !read(*,*)
    !end if
	if(idid.eq.-2)print*, "sigma0=",sigma
contains
	subroutine fcn(n, x, y, f, par, ipar)
		implicit none
		integer n, ipar(100)
		real(8) x, y(n), f(n), par(100), Fx
		call funcs(x*energyx, fx,xmin,xmax)
		f(1)=fx
        !if(energyx>10)then
        !    print*, "x,energxy,E,fx=",x,energyx,x*energyx,fx
        !end if
		!if(coeff_chattery.ge.3) 
        !print*, x, energyx,x*energyx, fx
		!if(ieee_is_nan(fx))then
		!	print*, alpha, x*energyx, fx
		!	stop
		!end if
	end subroutine
end subroutine
    