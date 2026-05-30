module md_event_datas
    use com_sts_type 
	
    type event_data_series
		!logical::same_mass=.true.
        integer n
		! real(8) mass_value
        real(8),allocatable:: x(:), jm(:)
		real(8),allocatable:: rp(:),ra(:),ac(:),ec(:), r_ls(:) ! AU
		real(8),allocatable:: w(:), mass(:),radius(:),inc(:)
        integer,allocatable:: exit_flag(:),nlvl(:)
    contains 
        procedure::init=>init_event_data_series
		procedure::read=>read_event_data_series
		procedure::write=>write_event_data_series
		!procedure::get_at_same_freq
    end type
	type event_data
		type(event_data_series)::emris, td
	contains
		procedure::read=>read_event_data
		procedure::write=>write_event_data
	end type
	type(event_data)::data_sbh,data_wd,data_ns,data_bd
	type(event_data)::data_star,data_rg
	integer,parameter:: nint_em=2, nreal_em=11
	private::init_event_data_series,read_event_data_series,write_event_data_series
	private::read_event_data,write_event_data	
contains
    subroutine init_event_data_series(em, n)
		implicit none
		class(event_data_series)::em
		integer n
		if(allocated(em%x))then
			deallocate(em%x,em%jm,em%rp,em%ra,em%ac,em%ec,&
				em%exit_flag,em%r_ls,em%mass, em%inc,em%radius, em%w,em%nlvl)
		end if
		allocate(em%x(n),em%jm(n),em%rp(n),em%ra(n),em%ac(n),em%ec(n),&
			em%r_ls(n),em%mass(n),em%w(n), em%radius(n), em%inc(n), em%exit_flag(n),em%nlvl(n))
		em%n=n
		!em%same_mass=same_mass
	end subroutine
    subroutine write_event_data_series(em,funit)
		implicit none
		class(event_data_series)::em
		!character*(*) fl
		integer i
		integer funit

		!open(unit=funit,file=trim(adjustl(fl))//"_emdata.bin",form='unformatted', &
		!	access='stream')
		write(funit) em%n
		
		if(em%n>0)then
			do i=1, em%n
				write(unit=funit) em%x(i), em%jm(i),em%rp(i), em%ra(i), &
					em%ac(i),em%ec(i),em%r_ls(i), em%w(i), em%mass(i), em%radius(i), em%inc(i), &
					em%exit_flag(i),em%nlvl(i)
			end do
		end if


		!close(unit=funit)
	end subroutine
	subroutine read_event_data_series(em,funit)
		implicit none
		class(event_data_series)::em
		!character*(*) fl
		integer i,funit
		!open(unit=funit,file=trim(adjustl(fl))//"_emdata.bin",form='unformatted', &
		!	access='stream',status="old")
		read(funit) em%n

		
		if(em%n>0)then
			call em%init(em%n)
			do i=1, em%n
				read(unit=funit) em%x(i), em%jm(i),em%rp(i), em%ra(i), &
					em%ac(i),em%ec(i),em%r_ls(i),em%w(i), em%mass(i), em%radius(i), em%inc(i), &
					em%exit_flag(i),em%nlvl(i)
			end do
		endif
		!close(unit=funit)
	end subroutine
	subroutine read_event_data(ed,funit)
		implicit none
		class(event_data)::ed
		integer funit
		call ed%emris%read(funit)
		call ed%td%read(funit)
	end subroutine
	subroutine write_event_data(ed,funit)
		implicit none
		class(event_data)::ed
		integer funit
		call ed%emris%write(funit)
		call ed%td%write(funit)
	end subroutine
	subroutine copy_i_event_data_series(em,idx,em_copy_to,idx_copy_to)
		implicit none
		type(event_data_series)::em,em_copy_to
		integer idx,idx_copy_to
		em_copy_to%ac(idx_copy_to)=em%ac(idx)
		em_copy_to%ec(idx_copy_to)=em%ec(idx)
		em_copy_to%x(idx_copy_to)=em%x(idx)
		em_copy_to%jm(idx_copy_to)=em%jm(idx)
		em_copy_to%rp(idx_copy_to)=em%rp(idx)
		em_copy_to%ra(idx_copy_to)=em%ra(idx)
		em_copy_to%r_ls(idx_copy_to)=em%r_ls(idx)
		em_copy_to%w(idx_copy_to)=em%w(idx)
		em_copy_to%mass(idx_copy_to)=em%mass(idx)
		em_copy_to%radius(idx_copy_to)=em%radius(idx)
		em_copy_to%inc(idx_copy_to)=em%inc(idx)
		em_copy_to%exit_flag(idx_copy_to)=em%exit_flag(idx)
		em_copy_to%nlvl(idx_copy_to)=em%nlvl(idx)
	end subroutine
	subroutine conv_em_int_real_arrays(em, intarr,realarr)
		implicit none
		type(event_data_series)::em
		real(8) realarr(nreal_em,em%n)
		integer i,intarr(nint_em,em%n)
		
		do i=1, em%n
			realarr(1:nreal_em,i)=(/em%x(i),em%jm(i),em%rp(i),em%ra(i), &
			em%ac(i),em%ec(i),em%r_ls(i),em%w(i), em%mass(i),em%radius(i),em%inc(i) /)
			intarr(1:nint_em,i)=(/em%exit_flag(i),em%nlvl(i)/)
		end do
	end subroutine
	
	subroutine conv_int_real_arrays_em(em, intarr,realarr)
		implicit none
		type(event_data_series)::em
		real(8) realarr(nreal_em,em%n)
		integer i,intarr(nint_em,em%n)
		do i=1, em%n
			em%x(i)=realarr(1,i)
			em%jm(i)=realarr(2,i)
			em%rp(i)=realarr(3,i)
			em%ra(i)=realarr(4,i)
			em%ac(i)=realarr(5,i)
			em%ec(i)=realarr(6,i)
			em%r_ls(i)=realarr(7,i)
			em%w(i)=realarr(8,i)
			em%mass(i)=realarr(9,i)
			em%radius(i)=realarr(10,i)
			em%inc(i)=realarr(11,i)
			em%exit_flag(i)=intarr(1,i)
			em%nlvl(i)=intarr(2,i)
			!em%star_flag(i)=intarr(2,i)
			!intarr(1:nint_em,i)=(/em%exit_flag(i)/)
		end do
	end subroutine
	
end module
