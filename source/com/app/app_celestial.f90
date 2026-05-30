!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!     MCO_EL2X.FOR    (ErikSoft  7 July 1999)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Calculates Cartesian coordinates and velocities given Keplerian orbital
! elements (for elliptical, parabolic or hyperbolic orbits).
!
!Based on a routine from Levison and Duncan's SWIFT integrator.
!
!  gm = grav const * (central + secondary mass)
! q = perihelion distance
!  e = eccentricity
! i = inclination                 )
!  p = longitude of perihelion !!! )   in
! n = longitude of ascending node ) radians
!  l = mean anomaly                )
!
!  x axis points to n=0
! 
!  x,y,z = Cartesian positions  ( units the same as a )
!  u,v,w =     "     velocities ( units the same as sqrt(gm/a) )
!
!------------------------------------------------------------------------------
!
      subroutine mco_el2x (gm,q,e,i,p,n,l,x,y,z,u,v,w)
!
      implicit none
!      include 'mercury.inc'
!
! Input/Output
      real*8 gm,q,e,i,p,n,l,x,y,z,u,v,w
!
! Local
      real*8 g,a,ci,si,cn,sn,cg,sg,ce,se,romes,temp
      real*8 z1,z2,z3,z4,d11,d12,d13,d21,d22,d23
      real*8 mco_kep, orbel_fhybrid, orbel_zget
!
!------------------------------------------------------------------------------
!
! Change from longitude of perihelion to argument of perihelion
      g = p - n
!
! Rotation factors
      call mco_sine (i,si,ci)
      call mco_sine (g,sg,cg)
      call mco_sine (n,sn,cn)
      z1 = cg * cn
      z2 = cg * sn
      z3 = sg * cn
      z4 = sg * sn
      d11 =  z1 - z4*ci
      d12 =  z2 + z3*ci
      d13 = sg * si
      d21 = -z3 - z2*ci
      d22 = -z4 + z1*ci
      d23 = cg * si
!
! Semi-major axis
      a = q / (1.d0 - e)
!
! Ellipse
      if (e.lt.1.d0) then
        romes = sqrt(1.d0 - e*e)
        temp = mco_kep (e,l)
        call mco_sine (temp,se,ce)
        z1 = a * (ce - e)
        z2 = a * romes * se
        temp = sqrt(gm/a) / (1.d0 - e*ce)
        z3 = -se * temp
        z4 = romes * ce * temp
      else
! Parabola
        if (e.eq.1.d0) then
          ce = orbel_zget(l)
          z1 = q * (1.d0 - ce*ce)
          z2 = 2.d0 * q * ce
          z4 = sqrt(2.d0*gm/q) / (1.d0 + ce*ce)
          z3 = -ce * z4
        else
! Hyperbola
          romes = sqrt(e*e - 1.d0)
          temp = orbel_fhybrid(e,l)
          call mco_sinh (temp,se,ce)
          z1 = a * (ce - e)
          z2 = -a * romes * se
          temp = sqrt(gm/abs(a)) / (e*ce - 1.d0)
          z3 = -se * temp
          z4 = romes * ce * temp
        end if
      endif
!
      x = d11 * z1  +  d21 * z2
      y = d12 * z1  +  d22 * z2
      z = d13 * z1  +  d23 * z2
      u = d11 * z3  +  d21 * z4
      v = d12 * z3  +  d22 * z4
      w = d13 * z3  +  d23 * z4
!
!------------------------------------------------------------------------------
!
      return
      end 
subroutine get_acc_of_orb(mg, r, x, acc)
      implicit none
      real(8) mg, r,x(3), acc(3)
      acc=-mg/r**3*x
end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MCO_SINE.FOR    (ErikSoft  17 April 1997)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Calculates sin and cos of an angle X (in radians).
!
!------------------------------------------------------------------------------
!
      subroutine mco_sine (x,sx,cx)
!
      implicit none
!
! Input/Output
      real*8 x,sx,cx
!
! Local
      real*8 pi,twopi
!
!------------------------------------------------------------------------------
!
      pi = 3.141592653589793d0
      twopi = 2.d0 * pi
!
      if (x.gt.0) then
        x = mod(x,twopi)
      else
        x = mod(x,twopi) + twopi
      end if
!
      cx = cos(x)
!
      if (x.gt.pi) then
        sx = -sqrt(1.d0 - cx*cx)
      else
        sx =  sqrt(1.d0 - cx*cx)
      end if
!
!------------------------------------------------------------------------------
!
      return
      end
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MCO_SINH.FOR    (ErikSoft  12 June 1998)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Calculates sinh and cosh of an angle X (in radians)
!
!------------------------------------------------------------------------------
!
      subroutine mco_sinh (x,sx,cx)
!
      implicit none
!
! Input/Output
      real*8 x,sx,cx
!
!------------------------------------------------------------------------------
!
      sx = sinh(x)
      cx = sqrt (1.d0 + sx*sx)
!
!------------------------------------------------------------------------------
!
      return
      end
!
!**********************************************************************
!                    ORBEL_FHYBRID.F
!**********************************************************************
!     PURPOSE:  Solves Kepler's eqn. for hyperbola using hybrid approach.  
!
!             Input:
!                           e ==> eccentricity anomaly. (real scalar)
!                           n ==> hyperbola mean anomaly. (real scalar)
!             Returns:
!               orbel_fhybrid ==>  eccentric anomaly. (real scalar)
!
!     ALGORITHM: For abs(N) < 0.636*ecc -0.6 , use FLON 
!	         For larger N, uses FGET
!     REMARKS: 
!     AUTHOR: M. Duncan 
!     DATE WRITTEN: May 26,1992.
!     REVISIONS: 
!     REVISIONS: 2/26/93 hfl
!**********************************************************************

	real*8 function orbel_fhybrid(e,n)

  !    include 'swift.inc'

!...  Inputs Only: 
	real*8 e,n

!...  Internals:
	real*8 abn
        real*8 orbel_flon,orbel_fget

!----
!...  Executable code 

	abn = n
	if(n.lt.0.d0) abn = -abn

	if(abn .lt. 0.636d0*e -0.6d0) then
	  orbel_fhybrid = orbel_flon(e,n)
	else 
	  orbel_fhybrid = orbel_fget(e,n)
	endif   

	return
	end  ! orbel_fhybrid
!**********************************************************************
!                    ORBEL_ZGET.F
!**********************************************************************
!     PURPOSE:  Solves the equivalent of Kepler's eqn. for a parabola 
!          given Q (Fitz. notation.)
!
!             Input:
!                           q ==>  parabola mean anomaly. (real scalar)
!             Returns:
!                  orbel_zget ==>  eccentric anomaly. (real scalar)
!
!     ALGORITHM: p. 70-72 of Fitzpatrick's book "Princ. of Cel. Mech."
!     REMARKS: For a parabola we can solve analytically.
!     AUTHOR: M. Duncan 
!     DATE WRITTEN: May 11, 1992.
!     REVISIONS: May 27 - corrected it for negative Q and use power
!	      series for small Q.
!**********************************************************************

	real*8 function orbel_zget(q)

   !   include 'swift.inc'

!...  Inputs Only: 
	real*8 q

!...  Internals:
	integer iflag
	real*8 x,tmp

!----
!...  Executable code 

	iflag = 0
	if(q.lt.0.d0) then
	  iflag = 1
	  q = -q
	endif

	if (q.lt.1.d-3) then
	   orbel_zget = q*(1.d0 - (q*q/3.d0)*(1.d0 -q*q))
	else
	   x = 0.5d0*(3.d0*q + sqrt(9.d0*(q**2) +4.d0))
	   tmp = x**(1.d0/3.d0)
	   orbel_zget = tmp - 1.d0/tmp
	endif

	if(iflag .eq.1) then
           orbel_zget = -orbel_zget
	   q = -q
	endif
	
	return
	end    ! orbel_zget
!----------------------------------------------------------------------
!**********************************************************************
!                    ORBEL_FGET.F
!**********************************************************************
!     PURPOSE:  Solves Kepler's eqn. for hyperbola using hybrid approach.  
!
!             Input:
!                           e ==> eccentricity anomaly. (real scalar)
!                        capn ==> hyperbola mean anomaly. (real scalar)
!             Returns:
!                  orbel_fget ==>  eccentric anomaly. (real scalar)
!
!     ALGORITHM: Based on pp. 70-72 of Fitzpatrick's book "Principles of
!           Cel. Mech. ".  Quartic convergence from Danby's book.
!     REMARKS: 
!     AUTHOR: M. Duncan 
!     DATE WRITTEN: May 11, 1992.
!     REVISIONS: 2/26/93 hfl
!     Modified by JEC
!**********************************************************************

	real*8 function orbel_fget(e,capn)

  !    include 'swift.inc'

!...  Inputs Only: 
	real*8 e,capn

!...  Internals:
	integer i,IMAX
	real*8 tmp,x,shx,chx
	real*8 esh,ech,f,fp,fpp,fppp,dx
	PARAMETER (IMAX = 10)

!----
!...  Executable code 

! Function to solve "Kepler's eqn" for F (here called
! x) for given e and CAPN. 

!  begin with a guess proposed by Danby	
	if( capn .lt. 0.d0) then
	   tmp = -2.d0*capn/e + 1.8d0
	   x = -log(tmp)
	else
	   tmp = +2.d0*capn/e + 1.8d0
	   x = log( tmp)
	endif

	orbel_fget = x

	do i = 1,IMAX
          call mco_sinh (x,shx,chx)
	  esh = e*shx
	  ech = e*chx
	  f = esh - x - capn
!	  write(6,*) 'i,x,f : ',i,x,f
	  fp = ech - 1.d0  
	  fpp = esh 
	  fppp = ech 
	  dx = -f/fp
	  dx = -f/(fp + dx*fpp/2.d0)
	  dx = -f/(fp + dx*fpp/2.d0 + dx*dx*fppp/6.d0)
	  orbel_fget = x + dx
!   If we have converged here there's no point in going on
	  if(abs(dx) .le. TINY) RETURN
	  x = orbel_fget
	enddo	

	!write(6,*) 'FGET : RETURNING WITHOUT COMPLETE CONVERGENCE' 
	return
	end   ! orbel_fget
!------------------------------------------------------------------
!***********************************************************************
!                    ORBEL_FLON.F
!**********************************************************************
!     PURPOSE:  Solves Kepler's eqn. for hyperbola using hybrid approach.  
!
!             Input:
!                           e ==> eccentricity anomaly. (real scalar)
!                        capn ==> hyperbola mean anomaly. (real scalar)
!             Returns:
!                  orbel_flon ==>  eccentric anomaly. (real scalar)
!
!     ALGORITHM: Uses power series for N in terms of F and Newton,s method
!     REMARKS: ONLY GOOD FOR LOW VALUES OF N (N < 0.636*e -0.6)
!     AUTHOR: M. Duncan 
!     DATE WRITTEN: May 26, 1992.
!     REVISIONS: 
!**********************************************************************

	real*8 function orbel_flon(e,capn)

 !     include 'swift.inc'

!...  Inputs Only: 
	real*8 e,capn

!...  Internals:
	integer iflag,i,IMAX
	real*8 a,b,sq,biga,bigb
	real*8 x,x2
	real*8 f,fp,dx
	real*8 diff
	real*8 a0,a1,a3,a5,a7,a9,a11
	real*8 b1,b3,b5,b7,b9,b11
	PARAMETER (IMAX = 10)
	PARAMETER (a11 = 156.d0,a9 = 17160.d0,a7 = 1235520.d0)
	PARAMETER (a5 = 51891840.d0,a3 = 1037836800.d0)
	PARAMETER (b11 = 11.d0*a11,b9 = 9.d0*a9,b7 = 7.d0*a7)
	PARAMETER (b5 = 5.d0*a5, b3 = 3.d0*a3)

!----
!...  Executable code 


! Function to solve "Kepler's eqn" for F (here called
! x) for given e and CAPN. Only good for smallish CAPN 

	iflag = 0
	if( capn .lt. 0.d0) then
	   iflag = 1
	   capn = -capn
	endif

	a1 = 6227020800.d0 * (1.d0 - 1.d0/e)
	a0 = -6227020800.d0*capn/e
	b1 = a1

!  Set iflag nonzero if capn < 0., in which case solve for -capn
!  and change the sign of the final answer for F.
!  Begin with a reasonable guess based on solving the cubic for small F	


	a = 6.d0*(e-1.d0)/e
	b = -6.d0*capn/e
	sq = sqrt(0.25*b*b +a*a*a/27.d0)
	biga = (-0.5*b + sq)**0.3333333333333333d0
	bigb = -(+0.5*b + sq)**0.3333333333333333d0
	x = biga + bigb
!	write(6,*) 'cubic = ',x**3 +a*x +b
	orbel_flon = x
! If capn is tiny (or zero) no need to go further than cubic even for
! e =1.
	if( capn .lt. TINY) go to 100

	do i = 1,IMAX
	  x2 = x*x
	  f = a0 +x*(a1+x2*(a3+x2*(a5+x2*(a7+x2*(a9+x2*(a11+x2))))))
	  fp = b1 +x2*(b3+x2*(b5+x2*(b7+x2*(b9+x2*(b11 + 13.d0*x2)))))   
	  dx = -f/fp
!	  write(6,*) 'i,dx,x,f : '
!	  write(6,432) i,dx,x,f
432	  format(1x,i3,3(2x,1p1e22.15))
	  orbel_flon = x + dx
!   If we have converged here there's no point in going on
	  if(abs(dx) .le. TINY) go to 100
	  x = orbel_flon
	enddo	

! Abnormal return here - we've gone thru the loop 
! IMAX times without convergence
	if(iflag .eq. 1) then
	   orbel_flon = -orbel_flon
	   capn = -capn
	endif
	!write(6,*) 'FLON : RETURNING WITHOUT COMPLETE CONVERGENCE' 
	  diff = e*sinh(orbel_flon) - orbel_flon - capn
	!  write(6,*) 'N, F, ecc*sinh(F) - F - N : '
	!  write(6,*) capn,orbel_flon,diff
	return

!  Normal return here, but check if capn was originally negative
100	if(iflag .eq. 1) then
	   orbel_flon = -orbel_flon
	   capn = -capn
	endif

	return
	end     ! orbel_flon
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MCO_KEP.FOR    (ErikSoft  7 July 1999)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Solves Kepler's equation for eccentricities less than one.
! Algorithm from A. Nijenhuis (1991) Cel. Mech. Dyn. Astron. 51, 319-330.
!
!  e = eccentricity
!  l = mean anomaly      (radians)
!  u = eccentric anomaly (   "   )
!
!------------------------------------------------------------------------------
!
      function mco_kep (e,oldl)
      implicit none
!
! Input/Outout
      real*8 oldl,e,mco_kep
!
! Local
      real*8 l,pi,twopi,piby2,u1,u2,ome,sign
      real*8 x,x2,sn,dsn,z1,z2,z3,f0,f1,f2,f3
      real*8 p,q,p2,ss,cc
      logical flag,big,bigg
!
!------------------------------------------------------------------------------
!
      pi = 3.141592653589793d0
      twopi = 2.d0 * pi
      piby2 = .5d0 * pi
!
! Reduce mean anomaly to lie in the range 0 < l < pi
      if (oldl.ge.0) then
        l = mod(oldl, twopi)
      else
        l = mod(oldl, twopi) + twopi
      end if
      sign = 1.d0
      if (l.gt.pi) then
        l = twopi - l
        sign = -1.d0
      end if
!
      ome = 1.d0 - e
!
      if (l.ge..45d0.or.e.lt..55d0) then
!
! Regions A,B or C in Nijenhuis
! -----------------------------
!
! Rough starting value for eccentric anomaly
        if (l.lt.ome) then
          u1 = ome
        else
          if (l.gt.(pi-1.d0-e)) then
            u1 = (l+e*pi)/(1.d0+e)
          else
            u1 = l + e
          end if
        end if
!
! Improved value using Halley's method
        flag = u1.gt.piby2
        if (flag) then
          x = pi - u1
        else
          x = u1
        end if
        x2 = x*x
        sn = x*(1.d0 + x2*(-.16605 + x2*.00761) )
        dsn = 1.d0 + x2*(-.49815 + x2*.03805)
        if (flag) dsn = -dsn
        f2 = e*sn
        f0 = u1 - f2 - l
        f1 = 1.d0 - e*dsn
        u2 = u1 - f0/(f1 - .5d0*f0*f2/f1)
      else
!
! Region D in Nijenhuis
! ---------------------
!
! Rough starting value for eccentric anomaly
        z1 = 4.d0*e + .5d0
        p = ome / z1
        q = .5d0 * l / z1
        p2 = p*p
        z2 = exp( log( dsqrt( p2*p + q*q ) + q )/1.5 )
        u1 = 2.d0*q / ( z2 + p + p2/z2 )
!
! Improved value using Newton's method
        z2 = u1*u1
        z3 = z2*z2
        u2 = u1 - .075d0*u1*z3 / (ome + z1*z2 + .375d0*z3)
        u2 = l + e*u2*( 3.d0 - 4.d0*u2*u2 )
      end if
!
! Accurate value using 3rd-order version of Newton's method
! N.B. Keep cos(u2) rather than sqrt( 1-sin^2(u2) ) to maintain accuracy!
!
! First get accurate values for u2 - sin(u2) and 1 - cos(u2)
      bigg = (u2.gt.piby2)
      if (bigg) then
        z3 = pi - u2
      else
        z3 = u2
      end if
!
      big = (z3.gt.(.5d0*piby2))
      if (big) then
        x = piby2 - z3
      else
        x = z3
      end if
!
      x2 = x*x
      ss = 1.d0
      cc = 1.d0
!
      ss = x*x2/6.*(1. - x2/20.*(1. - x2/42.*(1. - x2/72.*(1. - &
        x2/110.*(1. - x2/156.*(1. - x2/210.*(1. - x2/272.)))))))
      cc =   x2/2.*(1. - x2/12.*(1. - x2/30.*(1. - x2/56.*(1. - &
        x2/ 90.*(1. - x2/132.*(1. - x2/182.*(1. - x2/240.*(1. - &
        x2/306.))))))))
!
      if (big) then
        z1 = cc + z3 - 1.d0
        z2 = ss + z3 + 1.d0 - piby2
      else
        z1 = ss
        z2 = cc
      end if
!
      if (bigg) then
        z1 = 2.d0*u2 + z1 - pi
        z2 = 2.d0 - z2
      end if
!
      f0 = l - u2*ome - e*z1
      f1 = ome + e*z2
      f2 = .5d0*e*(u2-z1)
      f3 = e/6.d0*(1.d0-z2)
      z1 = f0/f1
      z2 = f0/(f2*z1+f1)
      mco_kep = sign*( u2 + f0/((f3*z1+f2)*z2+f1) )
!
!------------------------------------------------------------------------------
!
      return
      end
!
subroutine get_t_given_r_hypo(a,e,M, t_0, r, t)
      ! return time respect to t_0 for e>1 
      use constant
      implicit none
      real(8) a,e,M, t_0, r, t, cosf, mean_ano, ecc_ano
      real(8) coshE, mean_motion
      if(e<=1) then
            print*, 'get_t_given_r_hypo: e<=1'
            read(*,*)
      end if
     ! print*, "a,e,r=", a,e,r
      cosf=(a*(1-e**2)/r-1)/e
     ! print*, "cosf=",cosf
      coshE=(e+cosf)/(e*cosf+1)
      ecc_ano=acosh(coshE)
      mean_ano=e*sinh(ecc_ano)-ecc_ano
      mean_motion=sqrt(m/(-a)**3)
      t=t_0+mean_ano/mean_motion
   !   print*, "mean_motion, mean_ano, ecc_ano, coshE=", mean_motion, mean_ano, ecc_ano, coshE
end subroutine
subroutine get_ecc_ano_given_t(a,e,M,t,t_0,ecc_ano)
    use constant
	implicit none
	real(8) a,e,mean, ecc_ano
	real(8) mco_kep, n,t,t_0, M
	if(e.ge.1)then
		print*, "get_ano_given_t: e>1", e
        stop
	end if
	n=M**0.5/a**1.5d0
	mean=n*(t-t_0)
	ecc_ano=mco_kep (e,mean)

end subroutine
subroutine get_r_given_t(a,e,M,t,r)
    use constant
	implicit none
	real(8) a,e,mean, ecc_ano
	real(8) mco_kep, n,r,t, M
	if(e.ge.1)then
		print*, "get_r_given_t: e>1", e
        stop
	end if
	n=M**0.5/a**1.5d0
	mean=n*t
	ecc_ano=mco_kep (e,mean)
	r=a*(1-e*cos(ecc_ano))

end subroutine

subroutine get_r_given_mean(a,e,mean,r)
	implicit none
	real(8) a,e,mean, ecc_ano
	real(8) mco_kep, r

	if(e.ge.1)then
		print*, "error:get_r_given_mean: e>=1", e
		print*, "a,mean=", a,mean
 		stop        
	end if
	ecc_ano=mco_kep (e,mean)
	r=a*(1-e*cos(ecc_ano))
end subroutine

subroutine get_true_anomaly(mass, x,y,z,vx,vy,vz,ta)
	use constant
	implicit none
	!output ta in unit of rad
	real(8) x,y,z,vx,vy,vz,ta, vdot,rvdot,mass
	real(8) ve(3), vh(3), vr(3),vv(3),vtmp(3),vur(3),vue(3)
	vv=(/vx,vy,vz/)
	vr=(/x,y,z/)
!	print*, vr
!	print*, vv
	call vector_x(vr,vv, vh)
	call vector_x(vv,vh, vtmp)
	call vector_unit(vr,vur)
	ve=vtmp-vur*mass
	call vector_unit(ve,vue)
	call vector_dot(vue, vur,vdot)
	call vector_dot(vr,vv,rvdot)
!	print*, "vdot=",vdot, vue, vur
	if(vdot>1d0) vdot=1d0
	if(vdot<-1d0) vdot=-1d0
	if(rvdot>0)then
		ta=acos(vdot)
	else
		ta=2*pi-acos(vdot)
	end if
!	print*, ta, vdot
    if(isnan(ta))then
		print*, "vv, vr=", vv, vr
   	 	print*, "vur, vtmp=",vur,vtmp
    	print*,"vdot,vue,ve=", vdot, vue,ve
!       stop
    end if
end subroutine

      subroutine galactic_equatorial (Dir,Equinox,l,b, al,de)
!  If Dir: galactic -> equatorial coordinates transformation, if
!   .not.Dir: equatorial -> galactic transformation.
      implicit NONE
      real*8 pi, e,fa,fl, e2,fa2,fl2, e1,fa1,fl1, Equinox,l,b, al,de, w
      parameter ( pi = 3.141592653589793115997963468544185161591D+00 )
!  Parameters for Equinox=B1950.0:
      parameter ( e1 = 62.6d0 /180.d0*pi )
      parameter ( fa1 = 282.25d0 /180.d0*pi )
      parameter ( fl1 = 33.0d0 /180.d0*pi )
!  Parameters for Equinox=J2000.0:
      parameter ( e2 = 62.8717d0 /180.d0*pi )
      parameter ( fa2 = 282.8596d0 /180.d0*pi )
      parameter ( fl2 = 32.93192d0 /180.d0*pi )
      logical Dir
      
      if (Equinox .eq. 1950.0d0)  then
         e = e1
         fa = fa1
         fl = fl1
        elseif (Equinox .eq. 2000.0d0)  then
         e = e2
         fa = fa2
         fl = fl2
        else
         write (*,*) ' In galactic_equatorial only Equinox=1950.0'// &
                   ' or 2000.0 are supported!'
         stop         
        end if
         

      if (Dir)  then
         w = datan2(dsin(l-pi/2.d0-fl),(dcos(l-pi/2.d0-fl)*cos(e)- &
            dtan(b)*sin(e)))
         al = (w+fa-1.5d0*pi)
         de = dasin(sin(b)*cos(e)+cos(b)*sin(e)*cos(l-pi/2.d0-fl))
!  al - within interval [0...2*pi]
         al = al  - 2*pi*int(al/(2*pi))
         if (al .lt. 0.d0)  al = al + 2*pi
        else
		! print*, "here?"
         l = datan2(dcos(de)*dsin(al-fa)*dcos(e)+dsin(de)*dsin(e), &
                  dcos(de)*dcos(al-fa)) + fl
!  l - within interval [0...2*pi]
         l = l  - 2*pi*int(l/(2*pi))
         if (l .lt. 0.d0)  l = l + 2*pi
         b = dasin(dsin(de)*dcos(e)-dcos(de)*dsin(al-fa)*dsin(e))         
        end if

      return
      end
!============================================================================

      function Delta (Numbers,nn)
!  The function transforms three different formats of the equatorial coordinate
!  delta into degrees:
!   nn=1: dgr. -> dgr. (no change); {dgr.=Numbers(1)}
!   nn=2: dd mm.m -> dgr.; {dd=Numbers(1), mm.m=Numbers(2)}
!   nn=3: dd mm ss.s -> dgr. {dd=Numbers(1), mm=Numbers(2), ss.s=Numbers(3)}
!  The result is returned in Delta.
      implicit NONE
      integer nn
      real*8 Numbers(nn), Delta

      if (nn.gt.3 .or. nn.lt.1)  then
         write (*,*) ' Wrong nn in Delta!'
         stop
         end if

      if (nn .eq. 1)  then
         Delta = Numbers(1)
        else if (nn .eq. 2)  then
         Delta = dsign(1.d0,Numbers(1)) * (dabs(Numbers(1)) + &
                Numbers(2)/60.d0)
        else
         Delta = dsign(1.d0,Numbers(1)) * (dabs(Numbers(1)) + &
                      Numbers(2)/60.d0 + Numbers(3)/3600.d0)
        end if

      return
      end
!============================================================================


      INTEGER FUNCTION strlen_trail (STRING)
!
! Return the length of a string less any trailing blanks
!
!
!   STRING      CH*(*)  input   String to be measured
      implicit NONE
      CHARACTER*(*)     STRING
      INTEGER           I, N

      N = LEN (STRING)
      DO 1 I = N,1,-1
         IF(STRING(I:I).NE.' ') GOTO 2
    1 CONTINUE
    2 STRLEN_trail = I
      END
!--------------------------------------------------------------------------


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MCO_X2EL.FOR    (ErikSoft  23 January 2001)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Calculates Keplerian orbital elements given relative coordinates and
! velocities, and GM = G times the sum of the masses.
!
! The elements are: q = perihelion distance
!                   e = eccentricity
!                   i = inclination
!                   p = longitude of perihelion (NOT argument of perihelion!!)
!                   n = longitude of ascending node
!                   l = mean anomaly (or mean longitude if e < 1.e-8)
!
!------------------------------------------------------------------------------
!
      subroutine mco_x2el (gm,x,y,z,u,v,w,q,e,i,p,n,l)
!
      implicit none
!      include 'mercury.inc'
!
! Input/Output
      real*8 gm,q,e,i,p,n,l,x,y,z,u,v,w
!
! Local
      real*8 hx,hy,hz,h2,h,v2,r,rv,s,true
      real*8 ci,to,temp,tmp2,bige,f,cf,ce
      real*8 pi,twopi
!
!------------------------------------------------------------------------------
!
      pi = 3.141592653589793d0
      twopi = 2.d0 * pi      
!
!------------------------------------------------------------------------------
!
      hx = y * w  -  z * v
      hy = z * u  -  x * w
      hz = x * v  -  y * u
      h2 = hx*hx + hy*hy + hz*hz
      v2 = u * u  +  v * v  +  w * w
      rv = x * u  +  y * v  +  z * w
      r = sqrt(x*x + y*y + z*z)
      h = sqrt(h2)
      s = h2 / gm
!
! Inclination and node
      ci = hz / h
      if (abs(ci).lt.1) then
        i = acos (ci)
        n = atan2 (hx,-hy)
        if (n.lt.0) n = n + TWOPI
      else
        if (ci.gt.0) i = 0.d0
        if (ci.lt.0) i = PI
        n = 0.d0
      end if
!
! Eccentricity and perihelion distance
      temp = 1.d0  +  s * (v2 / gm  -  2.d0 / r)
      if (temp.le.0) then
        e = 0.d0
      else
        e = sqrt (temp)
      end if
      q = s / (1.d0 + e)
!
! True longitude
      if (hy.ne.0) then
        to = -hx/hy
        temp = (1.d0 - ci) * to
        tmp2 = to * to
        true = atan2((y*(1.d0+tmp2*ci)-x*temp),(x*(tmp2+ci)-y*temp))
      else
        true = atan2(y * ci, x)
      end if
      if (ci.lt.0) true = true + PI
!
      if (e.lt.3.d-8) then
        p = 0.d0
        l = true
      else
        ce = (v2*r - gm) / (e*gm)
!
! Mean anomaly for ellipse
        if (e.lt.1) then
          if (abs(ce).gt.1) ce = sign(1.d0,ce)
          bige = acos(ce)
          if (rv.lt.0) bige = TWOPI - bige
          l = bige - e*sin(bige)
        else
!
! Mean anomaly for hyperbola
          if (ce.lt.1) ce = 1.d0
          bige = log( ce + sqrt(ce*ce-1.d0) )
          if (rv.lt.0) bige = - bige
          l = e*sinh(bige) - bige
        end if
!
! Longitude of perihelion
        cf = (s - r) / (e*r)
        if (abs(cf).gt.1) cf = sign(1.d0,cf)
        f = acos(cf)
        if (rv.lt.0) f = TWOPI - f
        p = true - f
        p = mod (p + TWOPI + TWOPI, TWOPI)
      end if
!
      if (l.lt.0) l = l + TWOPI
      if (l.gt.TWOPI) l = mod (l, TWOPI)
!
!------------------------------------------------------------------------------
!
      return
      end

