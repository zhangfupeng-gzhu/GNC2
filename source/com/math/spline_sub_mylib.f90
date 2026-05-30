      SUBROUTINE spline_mylib(x,y,n,yp1,ypn,y2)
!	  to avoid name conflict      
		implicit none
      INTEGER n!,NMAX
      !PARAMETER (NMAX=40000000)
      real(8) yp1,ypn,x(n),y(n),y2(n)
      INTEGER i,k
      real(8) p,qn,sig,un,u(n)
      if (yp1.gt..99d30) then
        y2(1)=0.
        u(1)=0.
      else
        y2(1)=-0.5
        u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
	  if(x(2)<x(1))then
			print*, "in routine spline_mylib"
		   	print*,"x is ascending,x(2)<x(1), x(2), x(1)=",x(2),x(1) 
			pause
	  end if
      do 11 i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
		if(x(i+1)==x(i-1)) then
			print*, "in routine spline_mylib"
			print*, "x=",x
		   	print*,"x(i+1)=x(i-1), the x is not well ascending,i=",i 
			pause
		endif
        p=sig*y2(i-1)+2.
        y2(i)=(sig-1.)/p
        u(i)=(6d0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
11    continue
      if (ypn.gt..99d30) then
        qn=0.
        un=0.
      else
        qn=0.5
        un=(3d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
      do 12 k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
12    continue
      return
      END      
      
      SUBROUTINE splint_mylib(xa,ya,y2a,n,x,y,ier)
!	  to avoid the name conflict      
	implicit none
      INTEGER n
      real(8) x,y,xa(n),y2a(n),ya(n)
      INTEGER k,khi,klo,ier
      real(8) a,b,h
      ier=0
      klo=1
      khi=n
1     if (khi-klo.gt.1) then
        k=(khi+klo)/2d0
        if(xa(k).gt.x)then
          khi=k
	     ! write(*,*) "x,xa(k),khi=",x,xa(k),khi
		!elseif(xa(k).eq.x)then
		!	y=ya(k)
		!	return
        else
          klo=k
	      !write(*,*) "x,xa(k),klo=",x,xa(k),klo
        endif
      goto 1
      endif
     
      h=xa(khi)-xa(klo)
      !write(*,*) khi,klo
      if (h.eq.0.) then

      	!print*,"xa(khi)=xa(klo),khi,klo=", khi,klo
      	!print*,"xa(khi),xa(klo)=", xa(khi),xa(klo)
        !print*, "input x=",x," is not in the range of the data"
        !print*, "or input data is not ascending"
        !print*, "xa=",xa
      !	pause 'bad xa input in splint'
       !ier=1
        y=ya(khi)
        print*, "x,y=",x,y, ya(khi),ya(klo),khi,klo
        return 
      end if
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6d0
	!	print*, khi,klo
	!	print*, a, h, b
	!	print*, y, ya(khi), ya(klo)
		!print*, "ya(klo)*a+b*ya(khi)=",ya(klo)*a+b*ya(khi)
	!	print*, "rest=",y2a(khi)
      return
      END       
	   
	   
	   
	         
