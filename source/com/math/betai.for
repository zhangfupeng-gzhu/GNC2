      FUNCTION betai(a,b,x)
      REAL(8) betai,a,b,x
CU    USES betacf,gammln
      REAL(8) bt,betacf,gammln
      if(x.lt.0..or.x.gt.1.)pause 'bad argument x in betai'
      if(x.eq.0..or.x.eq.1.)then
        bt=0.
      else
        bt=exp(gammln(a+b)-gammln(a)-gammln(b)+a*log(x)+b*log(1.-x))
C       print*, gammln(a+b)-gammln(a)-gammln(b)+a*log(x)+b*log(1.-x)
      endif
      if(x.lt.(a+1.)/(a+b+2.))then
        betai=bt*betacf(a,b,x)/a
        return
      else
        betai=1.-bt*betacf(b,a,1.-x)/b
        return
      endif
      END
C  (C) Copr. 1986-92 Numerical Recipes Software :)z%+.
