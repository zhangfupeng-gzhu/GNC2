module md_hdf5_file
	use hdf5
	use h5lt
	use h5tb
	type hdf5_file_type
		INTEGER(HID_T) :: file_id      ! Group identifier		
		integer error
		character*(200) fd
	contains
		procedure::open=>hdf5_file_open
		procedure::close=>hdf5_file_close
	end type
	private::hdf5_file_open, hdf5_file_close
	integer,parameter::hdf5_file_flag_read=1,hdf5_file_flag_write=2,hdf5_file_flag_rw=3
contains
	subroutine hdf5_file_open(this, fd,flag)
		implicit none
		class(hdf5_file_type)::this
		character*(*) fd
		integer, optional::flag
		if(present(flag))then
			select case(flag)
			case(hdf5_file_flag_read)
				call h5open_file_read_hdf5(fd, this%file_id, H5F_ACC_RDONLY_F,this%error)
			case(hdf5_file_flag_rw)
				call h5open_file_read_hdf5(fd, this%file_id, H5F_ACC_RDWR_F,this%error)
			case(hdf5_file_flag_write)
				call h5create_file_open_hdf5(fd, this%file_id,this%error)
			end select
		else
			call h5create_file_open_hdf5(fd, this%file_id,this%error)
		end if
		
		this%fd=fd
	end subroutine
	subroutine hdf5_file_close(this)
		implicit none
		class(hdf5_file_type)::this
		call h5create_file_close_hdf5(this%file_id,this%error)
	end subroutine
end module
subroutine h5open_file_read_hdf5(fd, file_id,flag,error)
	use hdf5 
	use h5lt
	implicit none
	character*(*) fd
	INTEGER(HID_T) :: file_id
	integer error, flag
	CALL h5open_f(error)

	CALL h5fopen_f(trim(adjustl(fd)), flag, file_id, error)
end subroutine
subroutine h5create_file_open_hdf5(fd, file_id,error)
	use hdf5 
	use h5lt
	implicit none
	character*(*) fd
	INTEGER(HID_T) :: file_id
	integer error
	CALL h5open_f(error)
	CALL h5fcreate_f(trim(adjustl(fd)), H5F_ACC_TRUNC_F, file_id, error)
end subroutine
subroutine h5create_file_close_hdf5(file_id,error)
	use hdf5 
	use h5lt
	implicit none
	INTEGER(HID_T) :: file_id
	integer error
	CALL h5fclose_f(file_id, error)
	CALL h5close_f(error)
end subroutine