subroutine set_range(x,n,xmin,xmax,flag)
!	flag=0
!   |=x=|=x=|=x=|
!   flag=1
!   x=|=x=|=x=|=x  
	implicit none
	integer n,i,flag
	real(8) x(n),xmin,xmax
	real(8) xstep
	
	select case (flag)
	case (0)
		xstep=(xmax-xmin)/real(n)
		do i=1, n
			x(i)=xmin+xstep*(i-0.5d0)
		end do
	case (1)
		xstep=(xmax-xmin)/real(n-1)
		do i=1, n
			x(i)=xmin+xstep*(i-1)
		end do
	end select
end subroutine
subroutine get_dstr_num_in_each_bin(x, n, xbg, xstep, nbin, fx,  n_num)
	implicit none
	integer i, n, nbin
	real(8) x(n), w(n), xbg, xstep
	integer fx(nbin)
	integer n_num
	integer indx
	do i=1, n
		indx=int((x(i)-xbg)/xstep+1)
		if(indx>0.and.indx<=nbin)then
			fx(indx)=fx(indx)+1
			n_num=n_num+1
		end if
	end do
end subroutine
subroutine get_dstr_num_in_each_bin_weight(x, w, n, xbg, xstep, nbin, fxw, n_numw)
	implicit none
	integer i, n, nbin
	real(8) x(n), w(n), xbg, xstep
	real(8) fxw(nbin),  n_numw
	integer indx

	do i=1, n
		indx=int((x(i)-xbg)/xstep+1)
		if(indx>0.and.indx<=nbin)then
			fxw(indx)=fxw(indx)+w(i)
			n_numw=n_numw+w(i)
		end if
	end do

end subroutine



subroutine return_idx_ird(x,xb,xsteps,nx,idx,flag)
	!	flag=0
	!   |=x=|=x=|=x=|
	!   flag=1
	!   x=|=x=|=x=|=x  
		implicit none
		integer idx,idy
		integer flag
		integer nx,ny,i
		real(8) xsteps(nx),xb(nx), x,y
		real(8) xmin,xmax
	!	flag=0   sts_type_dstr
	!   |=x=|=x=|=x=|
	!   flag=1   sts_type_grid
	!   x=|=x=|=x=|=x  
	
		select case(flag)
		case (1)
			xmin=xb(1)-xsteps(1)/2d0
			xmax=xb(nx)+xsteps(nx)/2d0
			if(x.eq.xmax)then
				idx=nx
			elseif(x>xmax)then
				idx=-9999
			elseif(x<xmin)then
				idx=-9999
			elseif(x.eq.xmin)then
				idx=1
			else
				do i=1, nx
					if(xb(i)-xsteps(i)/2d0<x.and.xb(i)+xsteps(i)>x)then
						idx=i
					end if
				end do
			end if
		case (0)
			xmin=xb(1)
			xmax=xb(nx)
			if(x>xmax-xsteps(nx)/2d0.and.x.le.xmax)then
				idx=nx
			elseif(x>xmax)then
				idx=-9999
			elseif(x<xmin)then
				idx=-9999
			elseif(x.ge.xmin.and.x<xmin+xsteps(1)/2d0)then
				idx=1
			else
				do i=2, nx-1
					if(xb(i)-xsteps(i)/2d0<x.and.xb(i)+xsteps(i)/2d0>x)then
						idx=i
					end if
				end do
			end if
		end select
	end subroutine
subroutine return_idx_ir(x,n, xev,idx,flag)
	!	flag=0
	!   |=x=|=x=|=x=|
	!   flag=1
	!   x=|=x=|=x=|=x 
	implicit none
	integer n,flag
	real(8)x(n)
	real(8) xev
	integer idx
	integer ib,ie,mid
	ib=1; ie=n
10	if (ie-ib>1)then
		mid=(ib+ie)/2
		if(x(mid)>xev)then
		   ie=mid	
		else
		   ib=mid
		end if
		goto 10
	end if
	if(flag.eq.0)then
		if(xev<(x(ie)+x(ib))/2d0)then
			idx=ib
		else
			idx=ie
		end if
	else
		idx=ie
	end if

end subroutine

subroutine return_idx_ir_dstr(x,xstep, n, xev,idx,flag)
	!	flag=0
	!   |=x=|=x=|=x=|
	!   flag=1
	!   x=|=x=|=x=|=x 
	implicit none
	integer n,flag
	real(8)x(n), xstep(n)
	real(8) xev
	integer idx
	integer ib,ie,mid
	ib=1; ie=n
10	if (ie-ib>1)then
		mid=(ib+ie)/2
		if(x(mid)>xev)then
		   ie=mid	
		else
		   ib=mid
		end if
		goto 10
	end if
	if(flag.eq.1)then
		if(xev<x(ib)+xstep(ib)/2d0)then
			idx=ib
		else
			idx=ie
		end if
	else
		idx=ie
	end if

end subroutine

subroutine return_idx(x,xmin,xmax,nx,idx,flag)
!	flag=0
!   |=x=|=x=|=x=|
!   flag=1
!   x=|=x=|=x=|=x  
	implicit none
	integer idx,idy
	integer flag
	integer nx,ny
	real(8) xstep,ystep,x,y
	real(8) xmin,xmax,ymin,ymax
!	flag=0   sts_type_dstr
!   |=x=|=x=|=x=|
!   flag=1   sts_type_grid
!   x=|=x=|=x=|=x  

	select case(flag)
	case (1)
		xstep=(xmax-xmin)/real(nx-1)
		if(x.ge.xmin.and.x<xmax)then
			idx=nint((x-xmin)/xstep)+1
		else if(x.eq.xmax)then
			idx=nx
		else
			idx=-99999;
		end if
	case (0)
		xstep=(xmax-xmin)/real(nx)

		if(x>xmin.and.x<xmax)then
			idx=int((x-xmin)/xstep+1)
		elseif(x.eq.xmin)then
			idx=1
		elseif(x.eq.xmax)then
			idx=nx
		else
			idx=-99999;
		end if
	end select
end subroutine

subroutine return_rdx_ir(xb,xstep,x,idx,rdx,flag)
!	flag=0
!   |=x=|=x=|=x=|
!   flag=1
!   x=|=x=|=x=|=x  
	implicit none
	integer flag, n,idx
	real(8) xb,x,rdx,x0,x1
	real(8) xstep
!	flag=0   sts_type_dstr
!   |=x=|=x=|=x=|
!   flag=1   sts_type_grid
!   x=|=x=|=x=|=x  

	select case(flag)
	case (1)
		print*, "rdx:finish the code here"
		stop
	case (0)
		x0=xb-xstep/2d0
		!x1=xb+xstep/2d0
		rdx=-1+(x-x0)/xstep
	end select
end subroutine
	

subroutine Frequency_Count(x, weights, n, xbg, xend, rn, bin, pb, cb, fc, nfc, cc ,nca, ncum)
IMPLICIT NONE
integer n,rn
real(8) xbg, xend, incr
real(8) x(n)
real(8) pb(rn),cb(rn) !probability distribution, cumulative distribution
real(8) fc(rn),cc(rn), weights(n) !number of samples in each bin
real(8) bin(rn)
integer i,indx
real(8) nca
integer nfc(rn), ncum
	fc=0;nfc=0
	cc=0; ncum=0
	if(rn<1)then 
		pause 'error in FC'
	end if
	incr=(xend-xbg)/real(rn)
	do i=1, rn
		bin(i)=xbg+incr*(i-0.5d0)
	end do
	nca=0
	do i=1, n
		indx=int((x(i)-xbg)/incr+1)
		!print*, x(i),indx, xbg, incr
		if(indx>0.and.indx<=rn)then
			fc(indx)=fc(indx)+weights(i)
		    nfc(indx)=nfc(indx)+1
			nca=nca+weights(i)
            ncum=ncum+1
		end if
		!print*
	end do
	do i=1,rn
		pb(i)=fc(i)/nca/incr
	end do
	cc(1)=fc(1)
	cb(1)=pb(1)*incr
	do i=1,rn-1
		cc(i+1)=cc(i)+fc(i+1)
		cb(i+1)=cb(i)+pb(i+1)*incr
	end do
end subroutine

subroutine cal_freq_arr(x, n, xbg, xend, rn, bins, freqs,rfcs,nfcs,flagfc,flaglog)
	implicit none
	integer n,rn,i,flagfc,flaglog
	real(8) xbgl, xendl,nca
	real(8) x(n)
	real(8),allocatable:: weights(:)
	real(8) xbg,xend, rfcs(rn)
	real(8) bins(rn),freqs(rn)
	integer  nfcs(rn)
	
	allocate(weights(n))
	weights=1d0
	call cal_freq_arr_weights(x,weights, n, xbg, xend, rn, bins,freqs, rfcs, nfcs,flagfc,flaglog)

end subroutine
subroutine cal_freq_arr_weights(x,weights, n, xbg, xend, rn, bins, freqs, fcs, nfcs ,flagfc,flaglog)
	IMPLICIT NONE
	integer n,rn,i,flagfc,flaglog
	real(8) xbgl, xendl,nca
	real(8) x(n)
	real(8) xbg,xend
	real(8) bins(rn), freqs(rn)
	real(8) weights(n)
	real(8) fcs(rn)
    integer nfcs(rn), ncum
	real(8),allocatable:: fc(:),cc(:)
	real(8),allocatable::bin(:),xl(:),pb(:),cb(:)
    integer,allocatable:: nfc(:)

	allocate(bin(rn))
	allocate(fc(rn))
	allocate(cc(rn))
	allocate(xl(n))
	allocate(pb(rn))
	allocate(cb(rn))
    allocate(nfc(rn))
	
	if(flaglog.eq.1)then
		xbgl=log10(xbg)
		xendl=log10(xend)	
		xl=log10(x)	
	else
		xl=x;xbgl=xbg;xendl=xend
	end if
	
	call frequency_count(xl,weights, n, xbgl,xendl,rn, bin,pb,cb, fc,nfc, cc,nca, ncum)

	select case(flaglog)
	case (0)
		bins=bin
	case (1)
		bins=10**bin
	end select

	select case (flagfc)
	case (0)
		freqs=pb
		fcs=fc
	case (1)
		freqs=cb
		fcs=cc
	case (2)
		freqs=1-cb
		fcs=nca-cc
	end select
    nfcs=nfc

end subroutine
subroutine cal_freq(x, n, xbg, xend, rn, fout ,flagfc,flaglog )
	implicit none
	character*(*) fout
	integer n,rn,i,flagfc,flaglog
	real(8) xbgl, xendl,nca
	real(8) x(n)
	real(8),allocatable:: weights(:)
	real(8) xbg,xend
	
	allocate(weights(n))
	weights=1d0
	call cal_freq_weights(x,weights, n, xbg, xend, rn, fout ,flagfc,flaglog)
end subroutine

subroutine cal_freq_weights(x,weights, n, xbg, xend, rn, fout ,flagfc,flaglog)
	IMPLICIT NONE
	character*(*) fout
	integer n,rn,i,flagfc,flaglog, ncum
	real(8) xbgl, xendl,nca
	real(8) x(n)
	real(8) xbg,xend
	real(8),allocatable:: fc(:),cc(:)
	real(8) weights(n)
	real(8),allocatable::bin(:),xl(:),pb(:),cb(:)
    integer,allocatable::nfc(:)

	allocate(bin(rn))
	allocate(fc(rn))
	allocate(cc(rn))
	allocate(xl(n))
	allocate(pb(rn))
	allocate(cb(rn))
    allocate(nfc(rn))
	
	if(flaglog.eq.1)then
		xbgl=log10(xbg)
		xendl=log10(xend)	
		xl=log10(x)	
	else
		xl=x;xbgl=xbg;xendl=xend
	end if
	
	call frequency_count(xl,weights, n, xbgl,xendl,rn, bin,pb,cb, fc,nfc,cc,nca, ncum)
	open(unit=999,file=fout)
	if(flaglog.eq.1)then
		bin=10**bin
	end if

	do i=1,rn
	select case (flagfc)
	case (0)
		write(unit=999,fmt="(1P3E27.12)") bin(i),pb(i)
	case (1)
		write(unit=999,fmt="(1P3E27.12)") bin(i),cb(i)
	case (2)
		write(unit=999,fmt="(1P3E27.12)") bin(i),1-cb(i)
!---------------------------------------------------------
	case (3)
		write(unit=999,fmt="(1P3E27.12 )") bin(i),fc(i)
	case (4)
		write(unit=999,fmt="(1P3E27.12)") bin(i),cc(i)
	case (5)
		write(unit=999,fmt="(1P3E27.12)") bin(i),nca-cc(i)
!---------------------------------------------------------
	case (10)
		write(unit=999,fmt="(1P3E27.12E4, I10)") bin(i),pb(i),fc(i), nfc(i)
	case (11)
		write(unit=999,fmt="(1P3E27.12E4, I10)") bin(i),cb(i),cc(i), nfc(i)
	end select
	end do
	close(unit=999)
end subroutine
 

subroutine cal_freq_auto(arr,weights, n, rx, fileout,wflag,lflag)
	implicit none
	integer n,rx
	real(8) minx,maxx
	real(8) arr(n)
	real(8) weights(n)
	integer wflag, lflag
	character*(*) fileout
	
	minx=minval(arr)
	maxx=maxval(arr)
	minx=(1-abs(minx)/minx*0.2)*minx
	maxx=(1+abs(maxx)/maxx*0.2)*maxx
	
	call cal_freq_weights(arr,weights, n, minx,maxx,rx, fileout, wflag,lflag)
end subroutine
