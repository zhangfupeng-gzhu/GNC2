
subroutine bin2(x,y,n, xmin,xmax, abinx, ymin, ymax, abiny,  abin2D,flag)
!!����������2άbin��ͳ�Ƹ��������ڵ�������
!flag=0,the bin are in the right-top side
!flag=1,the bin are in the middle
implicit none
integer i, n,abinx, abiny,flag, idx, idy
real(8) x(n),y(n)
integer abin2D(abinx, abiny)
real(8) xmin,xmax, ymin, ymax
	abin2D=0
	do i=1, n
		call return_idxy(x(i),y(i),xmin,xmax,ymin,ymax,abinx,abiny,idx,idy,flag)
		if(idx>=1.and.idy>=1.and.idx<=abinx.and.idy<=abiny)then
			abin2D(idx,idy)=abin2D(idx,idy)+1
		end if
	end do
end subroutine

subroutine return_idxy(x,y,xmin,xmax,ymin,ymax,nx,ny,idx,idy,flag)
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
	if(flag==1)then
		xstep=(xmax-xmin)/real(nx-1)
		ystep=(ymax-ymin)/real(ny-1)
		!print*, "x, xmin, xmax=",x,xmin,xmax
		if(x>=xmin.and.x<xmax)then
			idx=nint((x-xmin)/xstep)+1
		!	print*, "idx=",idx
		else if(x.eq.xmax)then
			idx=nx
		else
			idx=-9999
		end if
		if(y>=ymin.and.y<ymax)then
			idy=nint((y-ymin)/ystep)+1
		else if(y.eq.ymax)then
			idy=ny
		else 
			idy=-9999
		end if
		!print*, "idx,idy=",idx,idy
		return
	end if

	if(flag==0)then
		xstep=(xmax-xmin)/real(nx)
		ystep=(ymax-ymin)/real(ny)
		if(x>xmin.and.x<xmax)then
			idx=int((x-xmin)/xstep+1)
		elseif(x.eq.xmax)then
			idx=nx
		elseif(x.eq.xmin)then
			idx=1
		else
			idx=-9999
		end if
		if(y>ymin.and.y<ymax)then
			idy=int((y-ymin)/ystep+1)
		elseif(y.eq.ymax)then
			idy=ny
		elseif(y.eq.ymin)then
			idy=1
		else 
			idy=-9999
		end if
	end if

end subroutine


subroutine return_idxy_ird(x,y,xb,xsteps,yb,ysteps,nx,ny,idx,idy,flag)
	implicit none
	integer idx,idy
	integer flag
	integer nx,ny
	real(8) xb(nx),yb(ny)
	real(8) xsteps(nx),ysteps(ny),x,y
	integer i,j
	real(8) xmin,xmax

!	flag=0   sts_type_dstr
!   |=x=|=x=|=x=|
!   flag=1   sts_type_grid
!   x=|=x=|=x=|=x  
	
	call return_idx_ird(x,xb,xsteps,nx,idx,flag)
	call return_idx_ird(y,yb,ysteps,ny,idy,flag)

end subroutine

subroutine cal_bin2_arr(x, y, n, xmin, xmax, rxn, ymin, ymax, ryn, xx,yy,fabin2D,bflag, mflag)
	implicit none
 	integer n, rxn,ryn,i,j,mflag,bflag
 	integer,allocatable:: abin2d(:,:)
	real(8) fabin2D(rxn,ryn)
	real(8) xx(rxn),yy(rxn), sumi
 	real(8) x(n),y(n),xmin,xmax,ymin,ymax,xstep,ystep
 	character*6  tmp
	!print*, rxn,ryn
 	allocate(abin2d(rxn,ryn))

	 call set_range(xx,rxn,xmin,xmax,bflag)
	 call set_range(yy,ryn,ymin,ymax,bflag)

 	call bin2(x, y, n, xmin,xmax, rxn, ymin, ymax, ryn, abin2d, bflag)	 	

 	sumi=0
 	do i=1, rxn
 		do j=1, ryn
 			sumi=sumi+abin2d(i,j)
 		end do
	end do 		

	xstep=xx(2)-xx(1)
	ystep=yy(2)-yy(1)
 	select case (mflag)
	case (0)
		fabin2D=dble(abin2D)
	case (1)		!normalized distribution function
		select case(bflag)
		case(0)
			fabin2D=dble(abin2D)/sumi/xstep/ystep
		case(1)

			do i=1, rxn
				do j=1, ryn
					if(i.eq.1.or.i.eq.rxn) then
						xstep=(xx(2)-xx(1))/2d0
					else
						xstep=(xx(2)-xx(1))
					end if
					if(j.eq.1.or.j.eq.ryn)then
						ystep=(yy(2)-yy(1))/2d0
					else
						ystep=yy(2)-yy(1)
					end if

					fabin2D(i,j)=dble(abin2d(i,j))/sumi/xstep/ystep
				end do
			end do
			
		end select
	case(2)    ! distribution function
		select case(bflag)
		case(0)
			fabin2D=dble(abin2D)/xstep/ystep
		case(1)
			do i=1, rxn
				do j=1, ryn
					if(i.eq.1.or.i.eq.rxn) then
						xstep=(xx(2)-xx(1))/2d0
					else
						xstep=(xx(2)-xx(1))
					end if
					if(j.eq.1.or.j.eq.ryn)then
						ystep=(yy(2)-yy(1))/2d0
					else
						ystep=yy(2)-yy(1)
					end if
					fabin2D(i,j)=dble(abin2d(i,j))/xstep/ystep
				end do
			end do
		end select
	end select

end subroutine

subroutine bin2_weight(x,y,w, n, xmin,xmax, abinx, ymin, ymax, abiny, abin2D,flag)
implicit none
integer n,abinx, abiny,flag, i, idx, idy
real(8) x(n),y(n), w(n)
real(8) abin2D(abinx, abiny)
real(8) xmin,xmax, ymin, ymax

	abin2D=0

	do i=1, n
		call return_idxy(x(i),y(i),xmin,xmax,ymin,ymax,abinx,abiny,idx,idy,flag)
		if(idx>=0.and.idy>=0.and.idx<=abinx.and.idy<=abiny)then
			abin2D(idx,idy)=abin2D(idx,idy)+w(i)
		end if
	end do
end subroutine

subroutine cal_bin2_arr_weight(x, y, weight, n, xmin, xmax, rxn, ymin, ymax, ryn, xx, yy, &
    fabin2D, bflag, mflag)
	implicit none
 	integer n, rxn,ryn,i,j,mflag, bflag
 	real(8),allocatable:: abin2d(:,:)
	real(8) fabin2D(rxn,ryn)
	real(8) xx(rxn),yy(rxn), sumi
 	real(8) x(n),y(n),weight(n),xmin,xmax,ymin,ymax,xstep,ystep
 	character*6  tmp
	!print*, rxn,ryn
 	allocate(abin2d(rxn,ryn))

	call set_range(xx,rxn,xmin,xmax,bflag)
	call set_range(yy,ryn,ymin,ymax,bflag)
 	call bin2_weight(x, y, weight, n, xmin,xmax, rxn, ymin, ymax, ryn,  abin2d, bflag)	 
	!print*, "abin2d"	
	!do i=1, rxn
	!	print*, abin2d(i,:)
	!end do
 	sumi=0
 	do i=1, rxn
 		do j=1, ryn
 			sumi=sumi+abin2d(i,j)
 		end do
	end do 		
	xstep=xx(2)-xx(1)
	ystep=yy(2)-yy(1)
 	select case (mflag)
	case (0)
		fabin2D=dble(abin2D)
	case (1)	!normalized distribution function
		select case(bflag)
		case(0)
			fabin2D=dble(abin2D)/sumi/xstep/ystep
		case(1)
			do i=1, rxn
				do j=1, ryn
					if(i.eq.1.or.i.eq.rxn) then
						xstep=(xx(2)-xx(1))/2d0
					else
						xstep=(xx(2)-xx(1))
					end if
					if(j.eq.1.or.j.eq.ryn)then
						ystep=(yy(2)-yy(1))/2d0
					else
						ystep=yy(2)-yy(1)
					end if

					fabin2D(i,j)=dble(abin2d(i,j))/sumi/xstep/ystep
				end do
			end do
		end select
	case(2)   ! distribution function
		select case(bflag)
		case(0)
			fabin2D=dble(abin2D)/xstep/ystep
		case(1)

			do i=1, rxn
				do j=1, ryn
					if(i.eq.1.or.i.eq.rxn) then
						xstep=(xx(2)-xx(1))/2d0
					else
						xstep=(xx(2)-xx(1))
					end if
					if(j.eq.1.or.j.eq.ryn)then
						ystep=(yy(2)-yy(1))/2d0
					else
						ystep=yy(2)-yy(1)
					end if

					fabin2D(i,j)=dble(abin2d(i,j))/xstep/ystep
				end do
			end do
		end select
	end select
	!print*, "fabin2d"	
	!do i=1, rxn
	!	print*, fabin2d(i,:)
	!end do

end subroutine

subroutine fcal_bin2(x, y, n, xmin, xmax, rxn, ymin, ymax, ryn, fout,bflag, mflag)
	implicit none
 	character*(*) fout
 	integer n, rxn,ryn,i,j, sumi,mflag,bflag
 	integer,allocatable:: abin2d(:,:)
	real(8),allocatable:: xbin(:),ybin(:),fabin2d(:,:)
 	real(8) x(n),y(n),xmin,xmax,ymin,ymax,xstep,ystep
 	character*6  tmp
 	
 	allocate(abin2d(rxn,ryn),fabin2d(rxn,ryn))
	allocate(xbin(rxn),ybin(ryn))

	call set_range(xbin,rxn,xmin,xmax,bflag)
	call set_range(ybin,ryn,ymin,ymax,bflag)

	call bin2(x, y, n, xmin,xmax, rxn, ymin, ymax, ryn, abin2d, bflag)	 	

	sumi=0
	do i=1, rxn
		do j=1, ryn
			sumi=sumi+abin2d(i,j)
		end do
   end do 		

   xstep=xbin(2)-xbin(1)
   ystep=ybin(2)-ybin(1)
	select case (mflag)
   case (0)
	   fabin2D=dble(abin2D)
   case (1)		
	   select case(bflag)
	   case(0)
		   fabin2D=dble(abin2D)/sumi/xstep/ystep
	   case(1)
		   do i=1, rxn
			   do j=1, rxn
				   if(i.eq.1.or.i.eq.rxn) then
					   xstep=(xbin(2)-xbin(1))/2d0
				   else
					   xstep=(xbin(2)-xbin(1))
				   end if
				   if(j.eq.1.or.j.eq.rxn)then
					   ystep=(ybin(2)-ybin(1))/2d0
				   else
					   ystep=ybin(2)-ybin(1)
				   end if

				   fabin2D(i,j)=dble(abin2d(i,j))/sumi/xstep/ystep
			   end do
		   end do
		   
	   end select
	   !if(isnan(fabin2D(1,1)))then
	   !	print*, "warnning:fabin2d nan:", fabin2D(1,1), abin2D(1,1), sumi, xstep, ystep
	   !end if
   end select

 	
 	open(unit=99,file=trim(adjustl(fout))//"_con.txt") 	

	write(tmp,fmt="(I4)") rxn	
	select case (mflag)
	case (0)
		do i=1,ryn
			write(unit=99,fmt="("//tmp//"I8)") abin2D(:,i)
		end do	
	case (1)
		do i=1,ryn
			write(unit=99,fmt="("//tmp//"F10.4)") fabin2D(:,i)
		end do	
	end select
	close(unit=99)
 	open(unit=999,file=trim(adjustl(fout))//"_x.txt") 	
	do i=1,rxn
		write(unit=999, fmt="(F20.10)") xbin(i)		
	end do
 	open(unit=999,file=trim(adjustl(fout))//"_y.txt") 	
	do i=1,ryn
		write(unit=999, fmt="(F20.10)") ybin(i)		
	end do

end subroutine



subroutine cal_bin2_arr_ird_weight(x, y, weight, n, xb, xsteps, nx, yb, ysteps, ny,  fabin2D, bflag, mflag)
	implicit none
 	integer n, nx,ny,i,j,mflag, bflag
 	real(8),allocatable:: abin2d(:,:)
	real(8) fabin2D(nx,ny)
	real(8) xb(nx),yb(ny), sumi
	real(8) xsteps(nx),ysteps(ny)
 	real(8) x(n),y(n),weight(n)
 	character*6  tmp
	integer idx,idy

 	allocate(abin2d(nx,ny))

	do i=1, n
		call return_idxy_ird(x(i),y(i),xb,xsteps,yb,ysteps,nx,ny,idx,idy,bflag)
	end do
	print*, "finished the code : cal_bin2_arr_ird_weight"
	stop

end subroutine
