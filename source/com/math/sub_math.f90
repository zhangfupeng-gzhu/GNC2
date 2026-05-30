subroutine math_xyz_sitafai(x,y,z,N,sita,fai)
integer n
real(8) x(n),y(n),z(n),sita(n),fai(n)
real(8),parameter:: PI=3.1415926535898
integer i
	!---------z=sin(sita)
	!---------y=cos(sita)*sin(fai)
	!---------x=cos(sita)*cos(fai)
	do i=1, n
		sita(i)=atan(z(i)/sqrt(x(i)*x(i)+y(i)*y(i)));
		if(x(i)>0..and.y(i)>=0.)then
			fai(i)=atan(y(i)/x(i))
		else if(x(i)<0..and.y(i)>=0.)then
			fai(i)=atan(y(i)/x(i))+PI
		else if(x(i)<0..and.y(i)<=0.)then
			fai(i)=atan(y(i)/x(i))+PI
		else if(x(i)>0..and.y(i)<=0.)then
			fai(i)=atan(y(i)/x(i))+PI*2
		end if
		if(x(i)==0..and.y(i)>=0.)then
			fai=PI/2
		end if
		if(x(i)==0..and.y(i)<=0.)then
			fai=PI*3/2
		end if
	end do
end subroutine

subroutine math_xyz_sitafai_one(x,y,z,sita,fai)
	real(8) x,y,z,sita,fai
	real(8),parameter:: PI=3.1415926535898
	integer i
	!---------z=sin(sita)
	!---------y=cos(sita)*sin(fai)
	!---------x=cos(sita)*cos(fai)
	
	sita=atan(z/sqrt(x*x+y*y));
	if(x>0..and.y>=0.)then
		fai=atan(y/x)
	else if(x<0..and.y>=0.)then
		fai=atan(y/x)+PI
	else if(x<0..and.y<=0.)then
		fai=atan(y/x)+PI
	else if(x>0..and.y<=0.)then
		fai=atan(y/x)+PI*2
	end if
	if(x==0..and.y>=0.)then
		fai=PI/2
	end if
	if(x==0..and.y<=0.)then
		fai=PI*3/2
	end if
	
end subroutine

subroutine math_xy_fai(x,y,fai)
	implicit none
	real(8) x,y,fai
	real(8),parameter:: PI=3.1415926535898d0
	integer i
	!---------y=sin(fai)
	!---------x=cos(fai)
		if(x>=0d0.and.y>=0d0)then
			fai=acos(x)
			return
		end if
		if(x>=0d0.and.y<=0d0)then
			fai=2*pi-acos(x)
			return
		end if
		if(x<=0d0)then
			fai=pi-asin(y)
		end if
end subroutine
subroutine math_sitafai_xyz(sita,fai,N,x,y,z)
! -pi/2<sita<pi/2, 0<fai<2pi
integer n
real(8) x(n),y(n),z(n),sita(n),fai(n)
real(8),parameter:: PI=3.1415926535898
integer i
	do i=1, n
		X(i)=cos(sita(i))*cos(fai(i))
		y(i)=cos(sita(i))*sin(fai(i))
		z(i)=sin(sita(i))
	end do
end subroutine

subroutine vector_x(v1,v2,vf)
        real(8) v1(3),v2(3),vf(3)
        vf(1)=v1(2)*v2(3)-v1(3)*v2(2)
        vf(2)=v1(3)*v2(1)-v1(1)*v2(3)
        vf(3)=v1(1)*v2(2)-v1(2)*v2(1)
end subroutine

subroutine vector_mag(v,mf)
		implicit none
        real(8) v(3),mf
		!if(v(1)**2+v(2)**2+v(3)**2.eq.0d0)then
		!	mf=0d0
		!	return
		!end if
        mf=sqrt(v(1)**2+v(2)**2+v(3)**2)
end subroutine

subroutine vector_m(v1,v2,vf)
        real(8) v1(3),v2(3),vf(3)
        vf(1)=v1(1)-v2(1)
        vf(2)=v1(2)-v2(2)
        vf(3)=v1(3)-v2(3)
end subroutine

subroutine vector_dot(v1,v2,mf)
        real(8) v1(3),v2(3),mf
        mf=v1(1)*v2(1)+v1(2)*v2(2)+v1(3)*v2(3)
end subroutine

subroutine vector_unit(v,vf)
	implicit none
	real(8) v(3),vf(3)
	real(8) vfm
	call vector_mag(v,vfm)
!	print*, "vfm=",vfm
	vf=v/vfm
end subroutine

subroutine vector2unit(v)
	implicit none
	real(8) v(3)
	real(8) vfm
	call vector_mag(v,vfm)
!	print*, "vfm=",vfm
	v=v/vfm
end subroutine

subroutine rotation_2D(vec1,vec2,theta)
	!vec2 is a rotation of vec1, theta is counter-clockwise
	real(8) vec1(2),vec2(2),theta
	vec2(1)=vec1(1)*cos(theta)-vec1(2)*sin(theta)
	vec2(2)=vec1(1)*sin(theta)+vec1(2)*cos(theta)
end subroutine

subroutine rotation_Eular(vecin, vecout, theta, phi, tau, direction)
	real(8) vecin(3), vecout(3), theta, phi,tau
	integer direction !-1 is the inverse transformation
	!Eular angle, if OO' is the line the two plane cross
	!theta: angle between zz';   phi: angle between  OO' and x;  tau angle between OO' and x'
	vecout(1)=cos(theta)*cos(phi)*vecin(2)+sin(phi)*vecin(1)+cos(phi)*sin(theta)*vecin(3)
	vecout(2)=cos(theta)*sin(phi)*vecin(2)-cos(phi)*vecin(1)+sin(phi)*sin(theta)*vecin(3)
	vecout(3)=-sin(theta)*vecin(2)+cos(theta)*vecin(3)
end subroutine

subroutine rotation_sph(vecin,vecout, theta,phi)
	real(8) vecin(3), vecout(3), theta, phi
	real(8) z_axis(3),vector_axis(3),vectmp(3)
!	z_axis=(/0,0,1/)
	call rotation_y(vecin,vecout,theta)
	vectmp=vecout
	call rotation_z(vectmp,vecout,phi)

end subroutine

subroutine rotation_around_vector(vecin, vecaxis, vecout, angle)
	implicit none
	real(8) vecin(3), vecaxis(3), vecout(3), angle	
	real(8) vnaxis(3), s, c
	real(8) x,y,z
!	print*, "-----------"
!	print*, vecaxis
	 call vector_unit(vecaxis, vnaxis)
		
	 s = sin(angle)
	 c = cos(angle)
	x=vnaxis(1);y=vnaxis(2);z=vnaxis(3)

	vecout(1)=(x**2*(1-c)+c)*vecin(1)+(x*y*(1-c)-z*s)*vecin(2) &
		+(x*z*(1-c)+y*s)*vecin(3)
	vecout(2)=(y*x*(1-c)+z*s)*vecin(1)+(y**2*(1-c)+c)*vecin(2) &
		+(y*z*(1-c)-x*s)*vecin(3)
	vecout(3)=(x*z*(1-c)-y*s)*vecin(1)+(y*z*(1-c)+x*s)*vecin(2) &
		+(z**2*(1-c)+c)*vecin(3)

!	print*, "----------"
end subroutine

subroutine rotation_x(vecin, vecout, angle)
	implicit none
	real(8) vecin(3), vecout(3), angle
	real(8) rx(3, 3)
	rx(1,1:3)=(/1d0, 0d0, 0d0/)
	rx(2,1:3)=(/0d0, cos(angle), sin(angle)/)
	rx(3,1:3)=(/0d0, -sin(angle),cos(angle)/)
	call matrix_13_33(vecin, rx, vecout)
end subroutine

subroutine rotation_y(vecin, vecout, angle)
	implicit none
	real(8) vecin(3), vecout(3), angle
	real(8) rx(3, 3)
	rx(1,1:3)=(/cos(angle), 0d0, -sin(angle)/)
	rx(2,1:3)=(/0d0, 1d0, 0d0/)
	rx(3,1:3)=(/sin(angle), 0d0,cos(angle)/)
	call matrix_13_33(vecin, rx, vecout)
	
end subroutine
subroutine rotation_z(vecin, vecout, angle)
	implicit none
	real(8) vecin(3), vecout(3), angle
	real(8) rx(3, 3)
	rx(1,1:3)=(/cos(angle), sin(angle), 0d0/)
	rx(2,1:3)=(/-sin(angle), cos(angle), 0d0/)
	rx(3,1:3)=(/0d0, 0d0,1d0/)
	call matrix_13_33(vecin, rx, vecout)
end subroutine
	
subroutine Matrix_33_33(a,b,c)
	real(8) a(3,3),b(3,3)
	real(8) c(3,3)
	integer i,j,k
	c=0
	do i=1,3
		do j=1,3
			do k=1,3
			c(i,j)=c(i,j)+a(i,k)*b(k,j)
			end do
		end do
	end do
end subroutine

subroutine Matrix_13_33(a,b,c)
!   j i---->
!   |
!   |
!   v
	real(8) a(3),b(3,3)
	real(8) c(3)
	integer i,j
	c=0
	do i=1,3
		do j=1,3
		c(i)=c(i)+a(j)*b(j,i)
		end do
	end do		
end subroutine

subroutine double_equation(A,X)
	implicit none
	real(8) A(3)	!a,b,c
	real(8) x(2)    !two solves
	real(8) sq
	sq=A(2)**2-4*A(1)*A(3)
	if(sq.ge.0)then
	X(1)=(-A(2)+sqrt(sq))/A(1)/2d0
	X(2)=(-A(2)-sqrt(sq))/A(1)/2d0
	else
		write(*,*) "sq<0",sq
	end if
end subroutine


subroutine coordtrasform_v_sph_car(vsph, sph, vcar)
	implicit none
	real(8) vsph(3), sph(3), vcar(3)
	real(8) mx(3,3)
	mx(1:3,1)=(/sin(sph(2))*cos(sph(3)),cos(sph(2))*cos(sph(3)), -sin(sph(3))/)
	mx(1:3,2)=(/sin(sph(2))*sin(sph(3)),cos(sph(2))*sin(sph(3)), cos(sph(3))/)
	mx(1:3,3)=(/cos(sph(2)),-sin(sph(2)), 0d0/)
	
	call matrix_13_33(vsph, mx, vcar)
	
end subroutine

subroutine coordtrasform_v_car_sph(vcar, sph, vsph)
	implicit none
	! it's vr, vtheta, vphi, not dot r, dot theta, dot phi
	real(8) vsph(3), sph(3), vcar(3)
	real(8) mx(3,3)
	mx(1:3,1)=(/sin(sph(2))*cos(sph(3)),sin(sph(2))*sin(sph(3)), cos(sph(2))/)
	mx(1:3,2)=(/cos(sph(2))*cos(sph(3)),cos(sph(2))*sin(sph(3)), -sin(sph(2))/)
	mx(1:3,3)=(/-sin(sph(3)),cos(sph(3)), 0d0/)
	
	call matrix_13_33(vcar, mx, vsph)
	
end subroutine
 
  

real(8) function fBPowerlawN(alpha,xb,n, t, c, q, x)   ! get the value of broken powerlaw function at x
	implicit none
	integer n, i
	real(8) alpha(n), c(n), xb(n+1), a(n), b(n), s(n+1),t(n+1)
	real(8) Q, x
	do i=1, n
		if(x<=xb(i+1).and.x>xb(i))then
			fBPowerlawN=c(i)*(x/xb(i))**alpha(i)
		end if
	end do
	
end function
real(8) function fCBPowerLawN(alpha,xb,n, t, c, q, x) ! get cumulative broken powerlaw function at x
	implicit none
	integer n, i
	real(8) alpha(n), c(n), xb(n+1), a(n), b(n), s(n+1),t(n+1)
	real(8) ymin, ymax, Q, x
	if(x.eq.xb(1))fCBPowerLawN=0
	do i=1, n
		if(x<=xb(i+1).and.x>xb(i))then
			if(alpha(i).ne.-1d0)then
				fCBPowerLawN=t(i)+c(i)/q*xb(i)/(alpha(i)+1)*((x/xb(i))**(alpha(i)+1)-1)
			else
				fCBPowerLawN=t(i)+c(i)/q*xb(i)*log(x/xb(i))
			end if 
		end if
	end do
end function

subroutine fCBPowerLawN_prepare(alpha,xb,n, c1, t,c,q) ! get cumulative broken powerlaw function at x
	implicit none
	integer n, i
	real(8) alpha(n), c(n), xb(n+1), a(n), b(n), s(n+1),t(n+1),c1
	real(8)  Q, x
	c(1)=c1
	do i=1, n-1
		c(i+1)=c(i)*(xb(i+1)/xb(i))**(alpha(i))     !c=f(x)
		if(xb(i+1)<xb(i))then
			print*, "fBpowerlawn_rnd:xb should be ascending"
			stop
		end if
	end do
	!print*, "c=", c
	Q=0d0
	s(1)=0
	do i=1, n
		if(alpha(i).ne.-1d0)then
			a(i)=c(i)*xb(i)/(1+alpha(i))*((xb(i+1)/xb(i))**(alpha(i)+1)-1)   ! a=int c dx
		else
			a(i)=c(i)*xb(i)*log(xb(i+1)/xb(i))
		end if
		!Q=Q+a(i)
		s(i+1)=s(i)+a(i)           ! s is for normalization
	end do
	!print*, "a=",a
	Q=S(n+1)
	!b=a/Q                          !
	t=S/Q							! t=F_cum(x)
	
end subroutine