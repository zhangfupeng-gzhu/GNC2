module md_hdf5
	use md_hdf5_table
	use md_hdf5_file
	use md_hdf5_group
	use md_hdf5_attr
	type hdf5_type
		type(hdf5_file_type)::h5f
		type(hdf5_group_type)::h5g
		type(hdf5_table_type)::h5t
	!	type(HID_T):: group_id_current
		contains
		procedure::open_file=>open_file_hdf5_type
		procedure::close_file=>close_file_hdf5_type
		procedure::create_group=>create_group_hdf5_type
		procedure::close_group=>close_group_hdf5_type
		!procedure::create_table=>create_table_hdf5_type		
	end type
	private::open_file_hdf5_type, create_group_hdf5_type, close_file_hdf5_type
	! create_table_hdf5_type
	contains
	subroutine open_file_hdf5_type(this, fd)
		implicit none
		class(hdf5_type)::this
		character*(*) fd
		call this%h5f%open(fd)
	end subroutine
	subroutine close_file_hdf5_type(this)
		implicit none
		class(hdf5_type)::this
		call this%h5f%close()
	end subroutine
	subroutine create_group_hdf5_type(this, gname, IDin)
		implicit none
		class(hdf5_type)::this
		integer(HID_T),optional:: IDin
		integer(HID_T) ID
		character*(*) gname

		if(present(IDin))then
			ID=IDin
		else
			ID=this%h5f%file_id
		end if
		call this%h5g%create(ID, gname)
	end subroutine
	subroutine close_group_hdf5_type(this)
		implicit none
		class(hdf5_type)::this
		call this%h5g%close()
	end subroutine
!	subroutine create_table_hdf5_type(this, nfin, nrin, tname)
!		implicit none
!		integer nfin, nrin
!		integer(hsize_t)::nf, nr
!		character*(*) tname
!		class(hdf5_type)::this
!		call this%h5t%create_table(nfin, nrin,tname)
!	end subroutine
    subroutine write_2d_arr(group_id, arr, nx,ny, aname)
        implicit none
        integer nx, ny
        integer(HID_T)::group_id
        character*(*) aname
        real(8) arr(nx,ny)
        integer error
        INTEGER(HSIZE_T), DIMENSION(2) :: dims2  ! Dataset dimensions

        dims2=(/nx,ny/)
        call h5ltmake_dataset_double_f(group_id, trim(adjustl(aname)), 2, dims2, &
             arr(:,:), error)
    end subroutine
end module