module md_dim2_cells
	type cell_type
		real(8), allocatable::y(:)
		real(8), allocatable::x(:)
		real(8),allocatable::w(:)   !weight
		integer,allocatable::idx(:)
        real(8) nw   !number using weight
		integer n
	end type
	type dm2s_type
		type(cell_type),allocatable::cell(:,:)
		real(8),allocatable::xedge(:),yedge(:)   !xbins(nx+1),ybins(ny+1)
		real(8),allocatable::xcenter(:),ycenter(:) !xcenter(nx), ycenter(ny)
		real(8) xstep,ystep, ntot
		real(8) ymin, ymax, xmin,xmax
		integer nx,ny
		integer flag
		contains
		procedure::get_from_data=>dm2s_get_from_data
		procedure::get_from_data_weight=>dm2s_get_from_data_weight
	end type
	integer,parameter::flag_along_xaxis=1, flag_along_yaxis=2
contains
	subroutine dm2s_init(dm,nx,ny, flag)
		implicit none
		type(dm2s_type)::dm
		integer nx,ny,flag

		if (allocated(dm%cell))then
			deallocate(dm%cell)
		end if
		allocate(dm%cell(nx,ny), dm%xedge(nx+1), &
				dm%yedge(ny+1), dm%xcenter(nx),dm%ycenter(ny))
		dm%nx=nx; dm%ny=ny
		if(flag.ne.1)then
			print*, "warning, flag.ne.1"
		end if
		dm%flag=flag
		dm%ntot=0
	end subroutine
	subroutine dm2s_init_cells(dm)
		implicit none
		type(dm2s_type)::dm
		integer i,j
		do i=1, dm%nx
			do j=1, dm%ny
				call dm2s_init_cell(dm%cell(i,j),dm%cell(i,j)%n)
			end do
		end do
	end subroutine	
	subroutine dm2s_init_cell(cell,n)
		implicit none
		type(cell_type)::cell
		integer n
		if (allocated(cell%x))then
			deallocate(cell%x,cell%y,cell%idx)
		end if
		allocate(cell%x(n),cell%y(n),cell%idx(n), cell%w(n))
		cell%n=n
	end subroutine		
	subroutine dm2s_get_frames(dm)
		implicit none
		type(dm2s_type)::dm
		integer i

		dm%xstep=(dm%xmax-dm%xmin)/real(dm%nx)
		dm%ystep=(dm%ymax-dm%ymin)/real(dm%ny)

		do i=1, dm%nx
			dm%xcenter(i)=dm%xmin+dm%xstep*(i-0.5d0)
		end do

		do i=1, dm%ny
			dm%ycenter(i)=dm%ymin+dm%ystep*(i-0.5d0)
		end do
		
		do i=1, dm%nx+1
			dm%xedge(i)=(dm%xmax-dm%xmin)/real(i-1)*real(dm%nx-1)+dm%xmin
		end do	

		do i=1, dm%ny+1
			dm%yedge(i)=(dm%ymax-dm%ymin)/real(i-1)*real(dm%ny-1)+dm%ymin
		end do	

	end subroutine

	subroutine dm2s_get_num_in_each_cell(dm,x,y,w,n)
		implicit none
		type(dm2s_type)::dm
		integer i,n,idx,idy
		real(8) x(n),y(n), w(n)
		dm%cell(:,:)%n=0
		do i=1, n
				if(x(i)>dm%xmin.and.x(i)<dm%xmax.and.y(i)>dm%ymin.and.y(i)<dm%ymax)then
					if(dm%flag==0)then
						idx=(x(i)-dm%xmin)/dm%xstep+2
						idy=(y(i)-dm%ymin)/dm%ystep+2
					else
						idx=int((x(i)-dm%xmin)/dm%xstep+1)
						idy=int((y(i)-dm%ymin)/dm%ystep+1)
					end if
					if(idx/=0.and.idy/=0.and.idx<=dm%nx.and.idy<=dm%ny)then
						dm%cell(idx,idy)%n=dm%cell(idx,idy)%n+1
						dm%cell(idx,idy)%nw=dm%cell(idx,idy)%nw+w(i)
						dm%ntot=dm%ntot+1						
					end if
				end if
		end do
	end subroutine
	subroutine dm2s_put_sam_in_cell(dm, x,y,w, n)
		type(dm2s_type)::dm
		integer i,n,idx,idy
		real(8) x(n),y(n),w(n)
		integer,allocatable:: ncell(:,:)
		allocate(ncell(dm%nx,dm%ny))
		ncell=0
		do i=1, n
			if(x(i)>dm%xmin.and.x(i)<dm%xmax.and.y(i)>dm%ymin.and.y(i)<dm%ymax)then
				if(dm%flag==0)then
					idx=(x(i)-dm%xmin)/dm%xstep+2
					idy=(y(i)-dm%ymin)/dm%ystep+2
				else
					idx=int((x(i)-dm%xmin)/dm%xstep+1)
					idy=int((y(i)-dm%ymin)/dm%ystep+1)
				end if
				if(idx/=0.and.idy/=0.and.idx<=dm%nx.and.idy<=dm%ny)then
					ncell(idx,idy)=ncell(idx,idy)+1
					dm%cell(idx,idy)%x(ncell(idx,idy))=x(i)
					dm%cell(idx,idy)%y(ncell(idx,idy))=y(i)
					dm%cell(idx,idy)%w(ncell(idx,idy))=w(i)
					dm%cell(idx,idy)%idx(ncell(idx,idy))=i
				end if
			end if
		end do

	end subroutine
	subroutine dm2s_get_from_data(dm,x,y,n, xmin,xmax,rxn, ymin, ymax, ryn)
		implicit none
		integer n,rxn,ryn
		real(8) x(n),y(n)
		integer abinx, abiny
		integer,allocatable::abin2D(:, :)
		real(8) xstep,ystep,xmin,xmax, ymin, ymax
		integer i,j
		integer,allocatable::pn(:)
		integer,allocatable::xbin(:),ybin(:)
		class(dm2s_type) dm
		real(8),allocatable:: w(:)
		allocate(w(n))
		w=1
		call dm2s_init(dm, rxn,ryn,1)
		dm%xmin=xmin;dm%xmax=xmax
		dm%ymin=ymin;dm%ymax=ymax

		call dm2s_get_frames(dm)
		call dm2s_get_num_in_each_cell(dm,x,y,w,n)
		call dm2s_init_cells(dm)
		call dm2s_put_sam_in_cell(dm,x,y,w, n)	
	end subroutine
	subroutine dm2s_get_from_data_weight(dm,x,y,w, n, xmin,xmax,rxn, ymin, ymax, ryn)
		implicit none
		integer n,rxn,ryn
		real(8) x(n),y(n), w(n)
		integer abinx, abiny
		integer,allocatable::abin2D(:, :)
		real(8) xstep,ystep,xmin,xmax, ymin, ymax
		integer i,j
		integer,allocatable::pn(:)
		integer,allocatable::xbin(:),ybin(:)
		class(dm2s_type) dm

		call dm2s_init(dm, rxn,ryn,1)
		dm%xmin=xmin;dm%xmax=xmax
		dm%ymin=ymin;dm%ymax=ymax

		call dm2s_get_frames(dm)
		call dm2s_get_num_in_each_cell(dm,x,y,w,n)
		call dm2s_init_cells(dm)
		call dm2s_put_sam_in_cell(dm,x,y,w,n)	
	end subroutine
	subroutine dm_get_pdf2d(dm, pdf2d, nx, ny)
		type(dm2s_type)::dm
		integer nx,ny,i,j
		real(8) pdf2d(nx,ny)
		do i=1, nx
			do j=1, ny
				pdf2d(i,j)=dm%cell(i,j)%n/real(dm%ntot)/dm%xstep/dm%ystep
			end do
		end do
	end subroutine

!	subroutine dm_d2tod1(dm, flag_direction, x, y, n)
!		type(dm2s_type)::dm
!		integer nx,ny,i,j
!		real(8) x(n),y(n)
!	
!		do i=1, nx
!			do j=1, ny
!				pdf2d(i,j)=dm%cell(i,j)%n/real(dm%ntot)/dm%xstep/dm%ystep
!			end do
!		end do
!
!	end subroutine


end module


