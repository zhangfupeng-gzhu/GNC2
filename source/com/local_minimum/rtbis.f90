      FUNCTION rtbis(func,x1,x2, xacc,par,ier,silent)
!	 using bisection, find the root of a function func known to lie between x1 and x2. 
!    the root, returned as rtbis, will be refined unitl its accuracy +-xacc.
      INTEGER JMAX
      REAL(8) rtbis,x1,x2,xacc,func,par(50)
      EXTERNAL func
      PARAMETER (JMAX=200)
      INTEGER j,ier
      logical silent
      REAL(8) dx,f,fmid,xmid
      fmid=func(x2,par)
      f=func(x1,par)
      ier=0
      if(f*fmid.ge.0.) then
		write(*,*) "x1,x2,f(x1,x2)=",x1,x2,f, fmid
    if(.not.silent)then
		  pause 'root must be bracketed in rtbis'
    else
      ier=1
    end if
!		error=-1
!		rtbits=-1d99
	  end if
      if(f.lt.0.)then
        rtbis=x1
        dx=x2-x1
      else
        rtbis=x2
        dx=x1-x2
      endif
      do 11 j=1,JMAX
        dx=dx*0.5
        xmid=rtbis+dx
        fmid=func(xmid,par)
        if(fmid.le.0.)rtbis=xmid
        if(abs(dx).lt.xacc .or. fmid.eq.0.) return
11    continue
      pause 'too many bisections in rtbis'
      END
!  (C) Copr. 1986-92 Numerical Recipes Software :)z%+.


FUNCTION rtbis_xrelacc(func,x1,x2, xracc,par)
  !	 using bisection, find the root of a function func known to lie between x1 and x2. 
  !    the root, returned as rtbis, will be refined unitl its accuracy +-yacc.
        implicit none
        INTEGER JMAX
        REAL(8) rtbis_xrelacc,x1,x2,xracc,func,par(50)
        EXTERNAL func
        PARAMETER (JMAX=200)
        INTEGER j
        REAL(8) dx,f,fmid,xmid
        fmid=func(x2,par)
        f=func(x1,par)
        if(f*fmid.ge.0.) then
      write(*,*) "x1,x2,f(x1,x2)=",x1,x2,f, fmid
      pause 'root must be bracketed in rtbis_xrelacc'
  !		error=-1
  !		rtbits=-1d99
      end if
        if(f.lt.0.)then
          rtbis_xrelacc=x1
          dx=x2-x1
        else
          rtbis_xrelacc=x2
          dx=x1-x2
        endif
        do 11 j=1,JMAX
          dx=dx*0.5
          xmid=rtbis_xrelacc+dx
          fmid=func(xmid,par)
          if(fmid.le.0.)rtbis_xrelacc=xmid
          print*, "dx/xmid=",dx,xmid,xracc
          if(abs(dx/xmid).lt.xracc .or. fmid.eq.0.) return
  11    continue
        pause 'too many bisections in rtbis_xrelacc'
END
  !  (C) Copr. 1986-92 Numerical Recipes Software :)z%+.
        
FUNCTION rtbis_yacc(func,x1,x2, yacc,par,niter,nmax, ier,silent)
  !	 using bisection, find the root of a function func known to lie between x1 and x2. 
  !    the root, returned as rtbis, will be refined unitl its accuracy +-yacc.
        implicit none
        INTEGER nmax
        REAL(8) rtbis_yacc,x1,x2,yacc,func,par(50)
        EXTERNAL func
        !PARAMETER (JMAX=1000)
        INTEGER j,ier,niter
        REAL(8) dx,f,fmid,xmid
        logical silent
        ier=0
        fmid=func(x2,par)
        f=func(x1,par)
        niter=0
        if(f*fmid.ge.0.) then
          if(f.eq.0)then
            rtbis_yacc=x1
            return
          end if
          if(fmid.eq.0)then
            rtbis_yacc=x2
            return
          end if
          if(.not.silent)then
            write(*,*) "x1,x2,f(x1,x2)=",x1,x2,f, fmid
            pause 'root must be bracketed in rtbis_yacc'
          else
            !write(*,*) "x1,x2,f(x1,x2)=",x1,x2,f, fmid
            ier=2
            return
          end if
  !		error=-1
  !		rtbits=-1d99
      end if
        if(f.gt.0.)then
          rtbis_yacc=x1
          dx=x1-x2
        else
          rtbis_yacc=x2
          dx=x2-x1
        endif
		!print*, "dx=",dx
        do 11 j=1,nMAX
          dx=dx*0.5d0
          niter=niter+1
          xmid=rtbis_yacc-dx
          fmid=func(xmid,par)
		  !print*, "xmid,xmid+dx=",xmid-dx*2,rtbis_yacc, func(xmid-dx*2,par), func(rtbis_yacc,par)
          if(fmid.ge.0.)then
            rtbis_yacc=xmid
            if(fmid.lt.yacc) then
              !print*, "fmid,yacc=",fmid,yacc
              return
            end if
          end if
  11    continue
        ier=1
        if(.not.silent)then
          print*, 'too many bisections in rtbis_xrelacc'
        else
          rtbis_yacc=xmid
          par(50)=dx
		  !print*, "xmid, dx=",xmid,dx
        end if
END
  !  (C) Copr. 1986-92 Numerical Recipes Software :)z%+.