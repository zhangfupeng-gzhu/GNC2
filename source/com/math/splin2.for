      SUBROUTINE splin2(x1a,x2a,ya,y2a,m,n,x1,x2,y)
      INTEGER m,n,NN,ier
      REAL(8) x1,x2,y,x1a(m),x2a(n),y2a(m,n),ya(m,n)
      PARAMETER (NN=1000000)
CU    USES spline,splint
!	  Given x1a, x2a, ya, m,n as described in splie2 and y2a
!     as produced by that routine; and given a desired interpolating point
!     x1,x2; this routine returns an interpolated function value y by
!     bicubic spline interpolation
      INTEGER j,k
      REAL(8) y2tmp(NN),ytmp(NN),yytmp(NN)
      do 12 j=1,m
        do 11 k=1,n
          ytmp(k)=ya(j,k)
          y2tmp(k)=y2a(j,k)
11      continue
        call splint_mylib(x2a,ytmp,y2tmp,n,x2,yytmp(j),ier)
12    continue
      call spline_mylib(x1a,yytmp,m,1.d30,1.d30,y2tmp)
      call splint_mylib(x1a,yytmp,y2tmp,m,x1,y,ier)
      return
      END
