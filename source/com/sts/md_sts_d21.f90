module md_sts_d21
	use md_dm2
	type sts_d21
		type(dm2),allocatable:: dy(:)
		real(8) xmin,xmax, ymin,ymax
        integer rxn
		logical ylimit
        contains
            procedure::init=>init_stsd21
	end type
    private::init_stsd21
contains
	subroutine init_stsd21_ylimit(sd,rxn,xmin,xmax, ymin,ymax)
		implicit none
		type(sts_d21) sd
		real(8) xmin,xmax, ymin,ymax
		integer rxn

		sd%ymin=ymin; sd%ymax=ymax
		sd%xmin=xmin; sd%xmax=xmax
		sd%rxn=rxn
		allocate(sd%dy(rxn))
		sd%ylimit=.true.
	end subroutine
	subroutine init_stsd21(sd,rxn,xmin,xmax)
		implicit none
		class(sts_d21) sd
		real(8) xmin,xmax, ymin,ymax
		integer rxn

		!sd%ymin=ymin; sd%ymax=ymax
		sd%xmin=xmin; sd%xmax=xmax
		sd%rxn=rxn
		if(allocated(sd%dy))then
			deallocate(sd%dy)
		end if
		allocate(sd%dy(rxn))
		sd%ylimit=.false.
		sd%dy(:)%n=0
	end subroutine
	subroutine get_d21_sts(x, y, n,sd)
		implicit none
		integer n,i
		real(8) x(n),y(n)
		type(sts_d21) sd

		!if(.not.sd%ylimit)then
		!	sd%ymin=minval(y)-0.1; sd%ymax=maxval(y)+0.1
		!endif
		!call d2to1(x,y,n,sd%xmin,sd%xmax,sd%rxn, sd%ymin,sd%ymax,sd%dy)
		do i=1, sd%rxn		
			call arravg(sd%dy(i)%y,sd%dy(i)%n, sd%dy(i)%avg)
			call arrsct(sd%dy(i)%y,sd%dy(i)%n, sd%dy(i)%avg, sd%dy(i)%sct)
		end do
	end subroutine
	subroutine get_d21(x, y, n,sd)
		implicit none
		integer n,i
		real(8) x(n),y(n)
		type(sts_d21) sd

		if(.not.sd%ylimit)then
			sd%ymin=minval(y); sd%ymax=maxval(y)
			sd%ymin=sd%ymin-abs(sd%ymin)*0.1
			sd%ymax=sd%ymax+abs(sd%ymax)*0.1
		endif
		call d2to1(x,y,n,sd%xmin,sd%xmax,sd%rxn, sd%ymin,sd%ymax,sd%dy)
		
	end subroutine
end module
