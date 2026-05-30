      FUNCTION rtnewt(funcd,x1,x2,xacc,par)
      INTEGER JMAX
      REAL(8) rtnewt,x1,x2,xacc,par(50)
      EXTERNAL funcd
      PARAMETER (JMAX=100)
      INTEGER j
      REAL(8) df,dx,f
      rtnewt=.5*(x1+x2)
      do 11 j=1,JMAX
        call funcd(rtnewt,f,df,par)
   		!write(*,*) "rtnewt=",rtnewt,"r=",par(1),"df=",df,"f=",f 
       dx=f/df
        rtnewt=rtnewt-dx
        if((x1-rtnewt)*(rtnewt-x2).lt.0.)then
				rtnewt=-1d99
				return
!	     		write(*,*) "rtnewt=",rtnewt,"r=",par(1),"df=",df,"f=",f
!				pause 'rtnewt jumped out of brackets'
		 end if
        if(abs(dx).lt.xacc) return
11    continue
	     write(*,*) "rtnewt=",rtnewt,"r=",par(1), "f=",f
      pause 'rtnewt exceeded maximum iterations'
      END
!  (C) Copr. 1986-92 Numerical Recipes Software :)z%+.

!      FUNCTION find_root_nwt(x,y,dy,n,x1,x2,xacc,par)
!      INTEGER JMAX
!	  integer n
!      REAL(8) find_root_nwt,x1,x2,xacc,x(n),y(n),dy(n),par(50)
!      PARAMETER (JMAX=100)
!      INTEGER j
!      REAL(8) df,dx,f
!	  real(8)::y2a(:),y2b(:)
!	  real(8),parameter::yp1=1d30;ypn=1d30
!
!      rtnewt=.5*(x1+x2)
!	  allocate(yp1(n))
!	  call spline_mylib(x,y,n,yp1,ypn,y2a)
!
!      do 11 j=1,JMAX
!        call funcd(rtnewt,f,df,par)
!		
!   		write(*,*) "rtnewt=",rtnewt,"r=",par(1),"df=",df,"f=",f 
!       dx=f/df
!        rtnewt=rtnewt-dx
!        if((x1-rtnewt)*(rtnewt-x2).lt.0.)then
!	     		write(*,*) "rtnewt=",rtnewt,"r=",par(1),"df=",df,"f=",f
!				pause 'rtnewt jumped out of brackets'
!		 end if
!        if(abs(dx).lt.xacc) return
!11    continue
!	     write(*,*) "rtnewt=",rtnewt,"r=",par(1)
!      pause 'rtnewt exceeded maximum iterations'
!      END
!!  (C) Copr. 1986-92 Numerical Recipes Software :)z%+.
