 
real(8) function GET_T_GW(m1,m2,a, e)
! T_GW in unit of Myr
	use constant
	implicit none
	real(8) m1,m2,a,e ! m in unit of solar mass, a in unit of AU
	real(8) fe
	fe=(1-e**2)**(-3.5d0)*(1+73d0/24d0*e**2+37d0/96d0*e**4)
	GET_T_GW=5d0/64d0*my_unit_vel_c**5*a**4/(m1*m2*(m1+m2))/fe/4d0  !!! the factor of 4 should be divided!!
	GET_T_GW=GET_T_GW/1d6/(2*pi)
end function 
 
real(8) function mdot_edd_msun_yr(Rad_eff,mbh)
	! input mbh in unit of sun
	! output mass accretion rate msun per yr
	implicit none
	real(8) mbh
	real(8) Rad_eff ! radiative efficiency
	mdot_edd_msun_yr=2.22*0.1d0/Rad_eff*mbh/1e8
end function 
 
subroutine get_c0(ain,ein,m1,m2,c0)
	implicit none
	real(8) ain,ein,m1,m2,c0
	c0=ain*(1-ein**2)/((1+121d0/304d0*ein**2)**(870d0/2299d0))/ein**(12d0/19d0)
end subroutine

subroutine get_fgw_orbfreq( m1, m2, ain, ein,fgworb)
	use constant
	implicit none
    real(8) ain, ein, fgworb, m1,m2
    fgworb=((m1+m2)/(ain**3d0))**0.5d0/2d0/pi
end subroutine
  

subroutine get_aeob_given_fgwc0(m1,m2, c0,fgw, agw, egw)
	! fgw in physics units 
	use constant
	implicit none
	real(8) ain,ein,m1,m2,agw,egw 
	real(8) par(10),c0, fgw
	real(8),external::rtbis
	real(8),parameter:: x0=1d-10,x1=1d0-1d-10
	real(8) f0,f1
	integer ier
!---
   ! print*, ain,ein,m1,m2,c0
   ! call get_c0(ain,ein,m1,m2,c0)
   ! print*, "1"
	!par(1)=c0
	f0=ffun(x0,par)
	f1=ffun(x1,par)
	if(f0*f1>0)then
        print*, "can not find agw, egw", f0, f1
        print*, "c0=", c0
        agw=-1d0;egw=2d0
       ! stop
		return
	end if
!	print*, "start"
	egw= rtbis(ffun, x0,x1, 1d-12, par,ier, .false.)
	call get_a_given_ce(c0, egw, agw)

contains
	real(8) function ffun(x,par)
		implicit none
		real(8) x, fg,fg2,df, par(10)
		real(8) h0,f
        
		call get_fgw_orbfreq_given_ce(m1,m2, c0, x, fg)
		ffun=log(fg)-log(fgw)
		!print*, "ffun, x,fg,fgw=",ffun,x,fg/(86400d0*365.2425d0/2d0/pi),fgw/(86400d0*365.2425d0/2d0/pi)
	end function
end subroutine

subroutine get_fgw_orbfreq_given_ce( m1, m2, c0, ein,fgworb)
	use constant
	implicit none
    real(8) c0, ein, fgworb, m1,m2, sigma
    
    sigma=ein**(12d0/19d0)/(1-ein**2)*((1+121d0/304d0*ein**2)**(870d0/2299d0))
    fgworb=((m1+m2)/((c0*sigma)**3d0))**0.5d0/2d0/pi
    !print*, "ein, c0, sigma=",ein, c0, sigma, (c0*sigma)**3d0, m1+m2
end subroutine
subroutine get_a_given_ce(c0, e0,a0)
	implicit none
	real(8) c0, a0, e0
	a0= c0/(1-e0**2)*(1+121d0/304d0*e0**2)**(870d0/2299d0)*e0**(12d0/19d0)
end subroutine
subroutine get_e_given_ca(c0, ain,m1,m2, e0)
	use constant
	implicit none
	real(8) e0,m1,m2 
	real(8) par(10), cf
	real(8),intent(in):: c0, ain
	real(8),external::rtbis
	real(8),parameter:: x0=1d-10,x1=1d0-1d-10
	real(8) f0,f1
	integer ier
!--- 
	f0=ffun(x0,par)
	f1=ffun(x1,par)
	if(f0*f1>0)then
        print*, "can not find agw, egw", f0, f1
        print*, "c0=", c0
        e0=2d0
       ! stop
		return
	end if
!	print*, "start"
	e0= rtbis(ffun, x0,x1, 1d-12, par,ier, .false.)
contains
	real(8) function ffun(x,par)
		implicit none
		real(8) x, fg,fg2,df, par(10)
		real(8) h0,f
        
		call get_c0(ain,x,m1,m2,cf)
		ffun=log(c0)-log(cf)
		!print*, "ffun, x,fg,fgw=",ffun,x,fg/(86400d0*365.2425d0/2d0/pi),fgw/(86400d0*365.2425d0/2d0/pi)
	end function
end subroutine
 