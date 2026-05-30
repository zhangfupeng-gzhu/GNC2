subroutine root_finding_2(func,par, xmin,xmax,r1,r2)
	! func is known to have two roots between xmin and xmax
	implicit none
	real(8),external::func
	real(8) xmin,xmax,xminloc,xmaxloc
	integer,parameter::n=20
	integer i, nb
	real(8) xb(n),yb(n),par(50)
	real(8) xb1(2),xb2(2)
	integer minlocate(1),minidx,maxidx,ier
	real(8) r1,r2,rtbis

	xminloc=xmin
	xmaxloc=xmax
100	do i=1, n
		xb(i)=(xmaxloc-xminloc)*real(i-1)/real(n-1)+xminloc
		yb(i)=func(xb(i),par)
		!print*, "xb(i),yb(i)=",xb(i),yb(i)
	end do
	!print*, "==="
	nb=0
	do i=2,n
		!print*, "yb(i-1),yb(i)=",xb(i), yb(i-1),yb(i)
		if(yb(i-1)*yb(i)<0)then
			nb=nb+1
			xb1(nb)=xb(i-1)
			xb2(nb)=xb(i)
			!print*, "yb(i-1),yb(i)=",yb(i-1),yb(i)
		end if
	end do
	if(nb.eq.0) then
		minlocate=minloc(yb**2)
		minidx=minlocate(1)-1
		maxidx=minlocate(1)+1
		xminloc=xb(max(minidx,1))
		xmaxloc=xb(min(maxidx,n))
		!print*, "xminloc,xmaxloc=",xminloc,xmaxloc
		goto 100
	else
		!print*, "r1,xb=",xb1(1),xb2(1)
		r1=rtbis(func,xb1(1),xb2(1),1d-10,par,ier,.true.)
		!print*, "r2,xb=",xb1(2),xb2(2)
		r2=rtbis(func,xb1(2),xb2(2),1d-10,par,ier,.true.)
	end if

end subroutine

subroutine root_finding_2_known_peak(func,par, xmin,xmid, xmax,inidx,r1,r2)
	! func is known to have two roots between xmin and xmax
	implicit none
	real(8),external::func
	real(8) xmin,xmax,xminloc,xmaxloc
	integer,parameter::n=20
	integer i, nb
	real(8) xmid, ymid,par(50), xtry, ytry
	real(8) xb1(2),xb2(2), inidx
	integer minlocate(1),minidx,maxidx,ier
	real(8) r1,r2,rtbis

	xminloc=xmin
	xmaxloc=xmax
	ymid=func(xmid,par)
	xtry=xmid-inidx
	ytry=func(xtry, par)
	do while (ymid*ytry>0)

	end do
	xb1(1)=xtry
	print*, "not yet finished"
	stop
	xtry=xmid+inidx
	do while (ymid*ytry>0)
		
	end do
	xb1(2)=xtry

	r1=rtbis(func,xb1(1),xmid,1d-10,par,ier,.true.)	
	r2=rtbis(func,xmid,xb2(2),1d-10,par,ier,.true.)

end subroutine

real(8) function root_finding_1(func, xmin,xmax,yacc,itermax,par,ier,silent)
	! func is known to have one root between xmin and xmax
	implicit none
	real(8),external::func
	real(8) xmin,xmax,xminloc,xmaxloc,yacc
	logical silent
	integer,parameter::n=20!, itermax=10
	integer i, nb,iter,itermax,niter,j
	real(8) xb(n),yb(n),par(50)
	real(8) xb1(2),xb2(2)
	integer minlocate(1),minidx,maxidx,  ier, ier_in
	real(8) rtbis_yacc

	xminloc=xmin
	xmaxloc=xmax
	ier=0
	iter=0
100	do i=1, n
		xb(i)=(xmaxloc-xminloc)*real(i-1)/real(n-1)+xminloc
		yb(i)=func(xb(i),par)
		!print*, "xb(i),yb(i)=",xb(i),yb(i)
	end do
	nb=0
	do i=2,n
		!print*, "yb(i-1),yb(i)=",xb(i), yb(i-1),yb(i)
		if(yb(i-1)*yb(i).le.0)then
			if(nb<=1)then
			nb=nb+1
			xb1(nb)=xb(i-1)
			xb2(nb)=xb(i)
			else
				do j=1, n
					xb(j)=(xmax-xmin)*real(j-1)/real(n-1)+xmin
					yb(j)=func(xb(j),par)				
					print*, "xb,yb=",xb(j),yb(j)
				end do
				stop
			end if
			!print*, "yb(i-1),yb(i)=",yb(i-1),yb(i)
		end if
	end do
	if(nb.eq.0) then
		iter=iter+1
		
		minlocate=minloc(yb**2)
		minidx=minlocate(1)-1
		maxidx=minlocate(1)+1
		xminloc=xb(max(minidx,1))
		xmaxloc=xb(min(maxidx,n))
		!print*, "xminloc,xmaxloc=",xminloc,xmaxloc
		if(iter<itermax)then
			goto 100
		else
			print*, "root_finding_1:error! maximum iter"
			print*, "can not find local root!,xminloc,xmaxloc=",xminloc,xmaxloc
			do i=1, n
				xb(i)=(xmax-xmin)*real(i-1)/real(n-1)+xmin
				yb(i)=func(xb(i),par)				
				print*, "xb,yb=",xb(i),yb(i)
			end do
			print*, "------------------------"
			do i=1, n
				xb(i)=(xmaxloc-xminloc)*real(i-1)/real(n-1)+xminloc
				yb(i)=func(xb(i),par)				
				print*, "xb,yb=",xb(i),yb(i)
			end do
			ier=1
			return
		end if
	else
		!print*, "r1,xb=",xb1(1),xb2(1)
		root_finding_1=rtbis_yacc(func,xb1(1),xb2(1),yacc,par,niter,1000,ier_in,silent)
		if(ier_in.eq.1)then
			if(iter<itermax)then
				iter=iter+1
				xminloc=xb1(1); xmaxloc=xb2(1)
				
				goto 100
			end if
		end if
		if(ier_in.eq.2.or.iter.eq.itermax)then
			if(.not.silent)then
				print*, "iter,ier,xlow,xup=",iter,ier_in,xb1(1),xb2(1)
				print*, "root_finding_1:error!, can not find local root!",xminloc,xmaxloc
				print*, "best=", root_finding_1, func(root_finding_1,par)
				do i=1, n
					xb(i)=(xmax-xmin)*real(i-1)/real(n-1)+xmin
					yb(i)=func(xb(i),par)				
					print*, "xb,yb=",xb(i),yb(i)
				end do
				print*, "------------------------"
				xminloc=xb1(1); xmaxloc=xb2(1)
				do i=1, n
					xb(i)=(xmaxloc-xminloc)*real(i-1)/real(n-1)+xminloc
					yb(i)=func(xb(i),par)				
					print*, "xb,yb=",xb(i),yb(i)
				end do
			end if
			ier=2
			return
		end if

		!print*, "r2,xb=",xb1(2),xb2(2)
		!r2=rtbis(func,xb1(2),xb2(2),1d-10,par)
	end if

end function
