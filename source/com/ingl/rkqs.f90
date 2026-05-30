      SUBROUTINE rkqs(y,dydx,n,x,htry,eps,yscal,hdid,hnext,derivs,errormethod, par,ipar)
!        Fifth-order Runge-Kutta step with monitoring of local truncation error to ensure accuracy
!        and adjust stepsize. Input are the dependent variable vector y(1:n) and its derivative
!        dydx(1:n) at the starting value of the independent variable x. Also input are the stepsize
!        to be attempted htry, the required accuracy eps, and the vector yscal(1:n) against
!        which the error is scaled. On output, y and x are replaced by their new values, hdid is the
!        stepsize that was actually accomplished, and hnext is the estimated next stepsize. derivs
!        is the user-supplied subroutine that computes the right-hand side derivatives.
      implicit none
      INTEGER n,NMAX
      REAL(8) eps,hdid,hnext,htry,x,dydx(n),y(n),yscal(n)
      real(8) par(100)
      integer ipar(100)
      PARAMETER (NMAX=50)
      external:: derivs,rkck
      INTEGER i, errormethod
      REAL(8) errmax,h,htemp,xnew,yerr(NMAX),ytemp(NMAX)
      real(8),parameter::tiny=1d-20
      real(8):: SAFETY=0.9,PGROW=-.2,PSHRNK=-.25,ERRCON=1.89e-4

      h=htry

1     call rkck(y,dydx,n,x,h,ytemp,yerr,derivs,par,ipar)
      errmax=0.
      select case(errormethod)
      case(0)
        do i=1,n
          errmax=max(errmax,abs(yerr(i)/yscal(i)))
        end do
        !print*, "yerr, yscal, errmax=", yerr(1:n), yscal(1:n), errmax
        errmax=errmax/eps
      case(1)
        do i=1, n
          yscal(i)=y(i)
          errmax=errmax+(yerr(i)/(eps*abs(yscal(i))+TINY))**2
        end do
        errmax=sqrt(errmax/n)
        print*, "errmax=", errmax,yerr(1:n), yscal(1:n)
       ! read(*,*)
      end select
      if(errmax.gt.1.)then
        htemp=SAFETY*h*(errmax**PSHRNK)
        h=sign(max(abs(htemp),0.1*abs(h)),h)
        !write(*,*) h
        xnew=x+h
        if(xnew.eq.x)pause 'stepsize underflow in rkqs'
        goto 1
      else
        if(errmax.gt.ERRCON)then
          hnext=SAFETY*h*(errmax**PGROW)
        else
          hnext=5.*h
        endif
        hdid=h
        x=x+h
        do 12 i=1,n
          y(i)=ytemp(i)
12      continue
          print*, "step success!", "x,y=", x,y(1:n)
          print*, "hdid, errmax, yerr, yscal=", hdid, errmax, yerr(1:2), yscal(1:2)
        return
      endif
      END
