module md_dm2
	type dm2
        integer,allocatable::idx(:)
		real(8), allocatable::y(:), x(:)
		real(8) ymin, ymax, center, avg, sct
		integer n
	end type
contains
	subroutine init_dm21(dm,n)
		implicit none
		type(dm2)::dm
		integer n
		if (allocated(dm%x))then
			deallocate(dm%x)
		end if
		if (allocated(dm%y))then
			deallocate(dm%y)
		endif
		if (allocated(dm%idx))then
			deallocate(dm%idx)
		endif
		allocate(dm%x(n),dm%y(n),dm%idx(n))
		dm%n=n
	end subroutine
	subroutine write_dm2_to_file_unit_bin(dm, file_unit)
		implicit none
		integer n,file_unit
		type(dm2)::dm

		write(unit=file_unit) dm%n
		n=dm%n
!		print*, "dm%n=",n
		write(unit=file_unit) dm%y(1:n),dm%x(1:n),dm%idx(1:n)
		write(unit=file_unit) dm%ymin,dm%ymax,dm%center
	end subroutine
	subroutine read_dm2_from_file_unit_bin(dm, file_unit)
		implicit none
		integer n,file_unit
		type(dm2)::dm
		
		read(unit=file_unit) dm%n
!		print*, "dm%n=",n
		n=dm%n
		call init_dm21(dm,dm%n)
		read(unit=file_unit) dm%y(1:n),dm%x(1:n),dm%idx(1:n)
		read(unit=file_unit) dm%ymin,dm%ymax,dm%center
	end subroutine

end module

subroutine d2to1(x,y,n, xmin,xmax,rxn, ymin, ymax, dy)
use md_dm2
implicit none
integer n,rxn
real(8) x(n),y(n)
integer idx,idy
integer,allocatable::abin2D(:, :)
real(8) xstep,ystep,xmin,xmax, ymin, ymax
integer i
real(8),allocatable:: xbin(:),ybin(:)
integer,allocatable::pn(:)
type(dm2) dy(rxn)

 	xstep=(xmax-xmin)/real(rxn)
!	ystep=(ymax-ymin)/real(abiny)

	allocate(abin2D(rxn, 1))
	!print*,"x=", x(1:10)
	call bin2(x,y,n, xmin,xmax,rxn,ymin,ymax, 1, abin2D,0)
	!print*, "xmin,xmax,ymin,ymax=",xmin,xmax,ymin,ymax
	do i=1, rxn
		!print*, "i,n=", i,abin2D(i,1)
		call init_dm21(dy(i), abin2D(i,1))
	end do

!	write(*,*) "pn allocate ..."
	allocate(pn(rxn))
!	write(*,*) "pn allocate success"
	pn=0
!	write(*,*) "rxn=",rxn
!	dy(1:rxn)%ymin=1d99
!	dy(1:rxn)%ymax=0
!	print*, "???4"
	do i=1, n
		call return_idxy(x(i),y(i),xmin,xmax,ymin,ymax,rxn,1,idx,idy,0)
		if(idy.eq.1)then
			if(idx>=1.and.idx<=rxn)then
				pn(idx)=pn(idx)+1
				dy(idx)%x(pn(idx))=x(i)
				dy(idx)%y(pn(idx))=y(i)				
				dy(idx)%idx(pn(idx))=i
			end if
		end if
	end do
	! checking 
	do i=1, rxn
		if(pn(i)<dy(i)%n)then
			print*, "d2to1: error! samples seems less than expected", i, pn(i),dy(i)%n
			stop
		end if
	end do
	!
	do i=1, rxn
		dy(i)%center=xmin+xstep*(i-0.5)
		dy(i)%ymin=minval(dy(i)%y)
		dy(i)%ymax=maxval(dy(i)%y)
	end do

end subroutine

module md_dm1
	type dm1
		real(8), allocatable::x(:)
		integer,allocatable::idx(:)
		real(8) center
        integer n
	end type
contains
	subroutine init_dm1(dm,n)
		implicit none
		type(dm1)::dm
		integer n

		if (allocated(dm%x))then
			deallocate(dm%x,dm%idx)
		end if
		allocate(dm%x(n),dm%idx(n))
		dm%n=n
	end subroutine
end module

subroutine d1to1(xin,n, xmin,xmax,rn,dy, flaglog)
use md_dm1
implicit none
integer n,rn
real(8) xin(n)
integer indx,nca,flaglog
real(8) incr,xbg,xed,xmin,xmax
integer i,j
integer,allocatable::pn(:)
real(8),allocatable::x(:)
type(dm1) dy(rn)

	allocate(x(n))
	if(rn<1)then 
		pause 'error in FC'
	end if
	
	if(flaglog.eq.1)then
		xbg=log10(xmin)
		xed=log10(xmax)	
		x=log10(xin)	
	else
		x=xin;xbg=xmin;xed=xmax
	end if
	
 	incr=(xed-xbg)/real(rn)
	
	nca=0
	dy(:)%n=0
	do i=1, n
		indx=int((x(i)-xbg)/incr)+1
		if(indx>0.and.indx<=rn)then
			dy(indx)%n=dy(indx)%n+1
			nca=nca+1
		end if
	end do

	do i=1, rn
		print*, "i=",i,dy(i)%n
		call init_dm1(dy(i), dy(i)%n)
	end do

	allocate(pn(rn))
	!write(*,*) "pn allocate success"
	pn=1

	do i=1, n
		if(x(i)>xbg.and.x(i)<xed)then
			indx=int((x(i)-xbg)/incr)+1
			if(indx>0.and.indx<=rn)then
				print*, "indx=",indx
				dy(indx)%x(pn(indx))=x(i)
				dy(indx)%idx(pn(indx))=i
				pn(indx)=pn(indx)+1
			end if
		end if
	end do

	do i=1, rn
		dy(i)%center=xbg+incr*(i-0.5)
	end do
	
	if(flaglog.eq.1)then
		do i=1, rn
			dy(i)%x=10**(dy(i)%x)	
			dy(i)%center=10**(dy(i)%center)
		end do
	end if	

end subroutine


