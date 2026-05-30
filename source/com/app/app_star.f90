real(8) function star_Radius(mass)
use constant
IMPLICIT NONE
real(8) mass
if(mass.le.0.06d0)then
	star_Radius=0.08*rd_sun  ! Chabrier, G. \& Baraffe, I.\ 2000, \araa, 38, 337.
elseif(mass.gt.0.06d0.and.mass.le.1)then
	star_Radius=mass**0.8d0*rd_sun  ! Stellar Structure and Evolution, Rudolf Kippenhahn, Alfred Weigert, Achim Weiss P253
else
	star_Radius=mass**0.56d0*rd_sun  ! Stellar Structure and Evolution, Rudolf Kippenhahn, Alfred Weigert, Achim Weiss P253
end if
end function
subroutine prepare_table_of_brown_dwarf_radius(age,xs,ys,y2)
	!use constant
	implicit none
	real(8) age
	!integer no
	integer,parameter::n=14
	real(8) xs(n),ys(n),y2(n)
	xs=(/0.01d0, 0.012d0,0.015d0,0.02d0,0.03d0,0.04d0,0.05d0,0.06d0,0.07d0,0.072d0,0.075d0,0.08d0,0.09d0,0.1d0/)
	if(age<1000)then
		ys=(/0.120d0,0.129d0,0.124d0,0.122d0,0.126d0,0.132d0,0.14d0,0.150d0,0.160d0,0.162d0,0.166d0,0.170d0,0.180d0,0.189d0/)
	elseif(age>=1000.and.age<5000d0)then !1Gyr
		ys=(/0.107d0,0.106d0,0.103d0,0.1d0,0.096d0,0.093d0,0.092d0,0.092d0,0.094d0,0.095d0,0.098d0,0.102d0,0.113d0,0.125d0/)
	elseif(age<10000d0)then  ! 5Gyr
		ys=(/0.101d0,0.100d0,0.098d0,0.095d0,0.090d0,0.085d0,0.082d0,0.079d0,0.081d0,0.083d0,0.089d0,0.099d0,0.113d0,0.125d0/)
	else   ! 10Gyr
		ys=(/0.099d0,0.098d0,0.096d0,0.093d0,0.088d0,0.083d0,0.079d0,0.076d0,0.078d0,0.081d0,0.089d0,0.099d0,0.113d0,0.125d0/)
	end if
	call spline_mylib(xs,ys,n,1d30,1d30,y2)
	!no=n
end subroutine
subroutine get_brown_dwarf_radius(mass,xs, ys, y2, radius)
	implicit none
	real(8) mass, radius
	integer,parameter::n=14
	real(8) xs(n),ys(n),y2(n)
	integer ier
	if(mass<0.01d0.or.mass>0.1d0)then
		print*, "error! brown dwarf mass is 0.001-0.1"
		stop
	end if
	call splint_mylib(xs,ys,y2,n,mass,radius,ier)

end subroutine
real(8) function white_dwarf_radius(mass)
	use constant
	implicit none
	real(8) mass
	if(mass<1.44d0)then
		white_dwarf_radius=0.01*rd_sun*mass**(-1/3d0)
	else
		print*, "error, white dwarf mass should be smaller than 1.44 solar mass"
	end if
end function 
  
real(8) function triple_break_IMF(t,c,q,m1,m2,alpha,xb)
	implicit none
	real(8) fCBPowerLawN, fBpowerlawn_part_rnd
	real(8) alpha(3), xb(4),t(4),c(3),q
	real(8) m1, m2,y1,y2

	if(m1<xb(1).or.m2>xb(4).or. m1>=m2)then
		print*, "error! m1, m2 should be in the range of 0.05-150 and m1 != m2"
	end if
	y1= fCBPowerLawN(alpha, xb, 3, t,c,q, m1)
	y2= fCBPowerLawN(alpha, xb, 3, t,c,q, m2)
	triple_break_IMF=fBpowerlawn_part_rnd(alpha,xb,3,t,c,q,y1,y2)
end function

real(8) function triple_break_IMF_func(t,c,q,m,alpha,xb)
	implicit none
	real(8) fbpowerlawn
	real(8) alpha(3), xb(4),t(4),c(3),q
	real(8)  m

	triple_break_IMF_func=fbpowerlawn(alpha,xb,3,t,c,q,m)
end function

subroutine triple_break_IMF_prepare(t,c,q,alpha,xb) 
	implicit none
	real(8) alpha(3), xb(4),t(4),c(3),q

	call fCBPowerLawN_prepare(alpha,xb,3,1d0,t,c,q)
end subroutine



subroutine triple_break_IMF_prepare_n(t,c,q,alpha, xb)  !Kroupa (2001)
	implicit none
	real(8) alpha(3), xb(4),t(4),c(3),q

	call fCBPowerLawN_prepare(alpha,xb,3,1d0,t,c,q)
end subroutine

subroutine triple_break_IMF_prepare_m(t,c,q,alpha_in, xb)  
	implicit none
	real(8) alpha(3),alpha_in(3), xb(4),t(4),c(3),q
	alpha=alpha_in+1d0
	call fCBPowerLawN_prepare(alpha,xb,3,xb(1), t,c,q)
end subroutine

real(8) function triple_break_IMF_func_nfrac_at_bin(t,c,q,m1b,m2b,alpha,xb)
	implicit none
	!according to Kroupa, P. 2001, mnras, 322, 231
	!function of m in linear bin
	!not normalized!, cannonical_IMF_func=1 at m=m2=0.5d0
	real(8) m1b, m2b
	real(8) t(4),c(3),q,alpha(3),xb(4), f1,f2, fcbpowerlawn

	f1= fCBPowerLawN(alpha, xb, 3, t,c,q, m1b)
	f2= fCBPowerLawN(alpha, xb, 3, t,c,q, m2b)
	triple_break_IMF_func_nfrac_at_bin=f2-f1
	!print*, "m1b,m2b,f1,f2=",m1b,m2b,f1,f2
end function

real(8) function triple_break_IMF_func_mfrac_at_bin(t,c,q,m1b,m2b,alpha_in,xb)
	implicit none
	!according to Kroupa, P. 2001, mnras, 322, 231
	!function of m in linear bin
	!not normalized!, cannonical_IMF_func=1 at m=m2=0.5d0
	real(8) m1b, m2b
	real(8) t(4),c(3),q,alpha(3),alpha_in(3),xb(4), f1,f2, fcbpowerlawn
	alpha=alpha_in+1

	f1= fCBPowerLawN(alpha, xb, 3, t,c,q, m1b)
	f2= fCBPowerLawN(alpha, xb, 3, t,c,q, m2b)
	triple_break_IMF_func_mfrac_at_bin=f2-f1
	!print*, "m1b,m2b,f1,f2=",m1b,m2b,f1,f2
end function


real(8) function canonical_IMF(t,c,q,m1,m2)  !Kroupa (2001)
	!! Note that here the mass range is from 0.01-150Msun!!
	!! according to Lockmann et al. 2010, MNRAS, 402, 519
	implicit none
	real(8) fCBPowerLawN, fBpowerlawn_part_rnd
	real(8) alpha(3), xb(4),t(4),c(3),q
	real(8) m1, m2,y1,y2
	!canonical_IMF=fBPowerlaw_rnd(-1.3d0, -2.3d0,0.5d0, 0.08d0, 150d0)
	!alpha=(/-1.3d0,-2.3d0/)
	!xb=(/0.08d0, 0.5d0, 150d0/)
	!canonical_IMF=fBpowerlawn_rnd(alpha,xb,2)
	
	alpha=(/-0.3D0, -1.3d0,-2.3d0/)
	xb=(/0.01d0, 0.08d0, 0.5d0, 150d0/)
	if(m1<xb(1).or.m2>xb(4).or. m1>=m2)then
		print*, "error! m1, m2 should be in the range of 0.05-150 and m1 != m2"
	end if
	y1= fCBPowerLawN(alpha, xb, 3, t,c,q, m1)
	y2= fCBPowerLawN(alpha, xb, 3, t,c,q, m2)
	canonical_IMF=fBpowerlawn_part_rnd(alpha,xb,3,t,c,q,y1,y2)
end function

real(8) function canonical_IMF_func(t,c,q,m)  !this is not a random function!
	implicit none
	!according to Kroupa, P. 2001, mnras, 322, 231
	!function of m in linear bin
	!not normalized!, cannonical_IMF_func=1 at m=m2=0.5d0
	real(8) t(4),c(3),q,alpha(3),xb(4), f1,f2, fbpowerlawn, m
	alpha=(/-0.3D0, -1.3d0,-2.3d0/)
	xb=(/0.01d0, 0.08d0, 0.5d0, 150d0/)
	canonical_IMF_func=fbpowerlawn(alpha,xb,3,t,c,q,m)
end function

subroutine canonical_IMF_prepare_n(t,c,q)  !Kroupa (2001)
	!! Note that here the mass range is from 0.01-150Msun!!
	!! according to Lockmann et al. 2010, MNRAS, 402, 519
	implicit none
	!real(8) fBPowerlaw_rnd, fBpowerlawn_rnd
	real(8) alpha(3), xb(4),t(4),c(3),q
	!canonical_IMF=fBPowerlaw_rnd(-1.3d0, -2.3d0,0.5d0, 0.08d0, 150d0)
	!alpha=(/-1.3d0,-2.3d0/)
	!xb=(/0.08d0, 0.5d0, 150d0/)
	!canonical_IMF=fBpowerlawn_rnd(alpha,xb,2)
	alpha=(/-0.3D0, -1.3d0,-2.3d0/)
	xb=(/0.01d0, 0.08d0, 0.5d0, 150d0/)
	call fCBPowerLawN_prepare(alpha,xb,3,1d0,t,c,q)
end subroutine

subroutine canonical_IMF_prepare_m(t,c,q)  !Kroupa (2001)
	!! Note that here the mass range is from 0.01-150Msun!!
	!! according to Lockmann et al. 2010, MNRAS, 402, 519
	implicit none
	!!real(8) fBPowerlaw_rnd, fBpowerlawn_rnd
	real(8) alpha(3), xb(4),t(4),c(3),q
	!canonical_IMF=fBPowerlaw_rnd(-1.3d0, -2.3d0,0.5d0, 0.08d0, 150d0)
	!alpha=(/-1.3d0,-2.3d0/)
	!xb=(/0.08d0, 0.5d0, 150d0/)
	!canonical_IMF=fBpowerlawn_rnd(alpha,xb,2)
	alpha=(/-0.3D0, -1.3d0,-2.3d0/)+1d0
	xb=(/0.01d0, 0.08d0, 0.5d0, 150d0/)
	call fCBPowerLawN_prepare(alpha,xb,3,xb(1), t,c,q)
end subroutine



real(8) function canonical_IMF_func_nfrac_at_bin(t,c,q,m1b,m2b)
	implicit none
	!according to Kroupa, P. 2001, mnras, 322, 231
	!function of m in linear bin
	!not normalized!, cannonical_IMF_func=1 at m=m2=0.5d0
	real(8) m1b, m2b
	real(8) t(4),c(3),q,alpha(3),xb(4), f1,f2, fcbpowerlawn
	alpha=(/-0.3D0, -1.3d0,-2.3d0/)
	xb=(/0.01d0, 0.08d0, 0.5d0, 150d0/)

	f1= fCBPowerLawN(alpha, xb, 3, t,c,q, m1b)
	f2= fCBPowerLawN(alpha, xb, 3, t,c,q, m2b)
	canonical_IMF_func_nfrac_at_bin=f2-f1
	!print*, "m1b,m2b,f1,f2=",m1b,m2b,f1,f2
end function

real(8) function canonical_IMF_func_mfrac_at_bin(t,c,q,m1b,m2b)
	implicit none
	!according to Kroupa, P. 2001, mnras, 322, 231
	!function of m in linear bin
	!not normalized!, cannonical_IMF_func=1 at m=m2=0.5d0
	real(8) m1b, m2b
	real(8) t(4),c(3),q,alpha(3),xb(4), f1,f2, fcbpowerlawn
	alpha=(/-0.3D0, -1.3d0,-2.3d0/)+1
	xb=(/0.01d0, 0.08d0, 0.5d0, 150d0/)

	f1= fCBPowerLawN(alpha, xb, 3, t,c,q, m1b)
	f2= fCBPowerLawN(alpha, xb, 3, t,c,q, m2b)
	canonical_IMF_func_mfrac_at_bin=f2-f1
	!print*, "m1b,m2b,f1,f2=",m1b,m2b,f1,f2
end function