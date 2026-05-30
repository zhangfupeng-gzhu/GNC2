      SUBROUTINE odeint(ystart,nvar,x1,x2,eps,h1,hmin,nok,nbad,derivs, &
     rkqs,my_solout, iwork,  iout, par,  ipar)

      implicit none
      INTEGER nbad,nok,nvar,KMAXX,MAXSTP,NMAX
      REAL(8) eps,h1,hmin,x1,x2,ystart(nvar),TINY
      EXTERNAL derivs,rkqs,my_solout
      PARAMETER (MAXSTP=100000000,NMAX=50,KMAXX=200,TINY=1.e-20)
      INTEGER i,kmax,kount,nstp, iwork(10)
      ! iwork(1) specify the method of calculating yscal(i)
      ! iwork(1)=0  good for equations when y(i) is always positive or negative.

      REAL(8) dxsav,h,hdid,hnext,x,xsav,dydx(NMAX),xp(KMAXX),y(NMAX)
      real(8) yp(NMAX,KMAXX),yscal(NMAX)
      real(8) par(100),xout
      integer ipar(100),IRTRN
      integer,parameter::nd=100
      COMMON /path/ kmax,kount,dxsav,xp,yp
      real(8) CON(5*ND),ICOMP(ND)
      integer::IOUT

      x=x1
      h=sign(h1,x2-x1)
      nok=0
      nbad=0
      kount=0
      do i=1,nvar
        y(i)=ystart(i)
      end do
      select case(iout)
        case(0)
          if (kmax.gt.0) xsav=x-2.*dxsav
        case(1)
          xsav=x
      end select

      do nstp=1,MAXSTP
        call derivs(nvar, x,y, dydx,par,ipar)
        select case (iwork(1))
        case(0)
          do i=1,nvar
            yscal(i)=abs(y(i))+abs(h*dydx(i))+TINY
          end do
        case(1)
         ! do i=1,nvar
         !   yscal(i)=abs(dydx(i))
          !end do
        end select
        select case(iout)
        case(0)
          if(kmax.gt.0)then
            if(abs(x-xsav).gt.abs(dxsav)) then
              if(kount.lt.kmax-1)then
                kount=kount+1
                xp(kount)=x
                do i=1,nvar
                  yp(i,kount)=y(i)
                end do
                xsav=x
              endif
            endif
          endif
        case(1)
          if(x.ne.xsav) then
            call my_solout(nok, xsav,X,Y,Nvar,CON, ICOMP,nd, PAR,IPAR,IRTRN,xout)
            if(irtrn.eq.-1)then
              return
            end if
            xsav=x
          end if
        end select
        if((x+h-x2)*(x+h-x1).gt.0.) h=x2-x
        call rkqs(y,dydx,nvar,x,h,eps,yscal,hdid,hnext,derivs,iwork(1), par,ipar)
        if(hdid.eq.h)then
          nok=nok+1
        else
          nbad=nbad+1
        endif
        if((x-x2)*(x2-x1).ge.0.)then
          do i=1,nvar
            ystart(i)=y(i)
          end do
          if(kmax.ne.0)then
            kount=kount+1
            xp(kount)=x
            do i=1,nvar
              yp(i,kount)=y(i)
            end do
          endif
          return
        endif
        if(abs(hnext).lt.hmin) then
          print*, "hmin, hnext=", hmin, hnext
          print*, 'stepsize smaller than minimum in odeint'
          read(*,*)
        end if
        h=hnext     
      end do
      print*, "nok, nbad=", nok, nbad
      pause 'too many steps in odeint'
      return
END subroutine
