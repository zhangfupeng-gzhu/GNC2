module md_hdf5_table
	use hdf5
	use h5lt
	use h5tb
	type hdf5_table_type
		INTEGER(HSIZE_T) :: nfields 
		INTEGER(HSIZE_T) :: nrecords 
		integer ::si
		INTEGER(HID_T) :: group_id      ! Group identifier		
		CHARACTER*(200),allocatable :: field_names(:)  ! field names
		INTEGER(SIZE_T),allocatable :: field_sizes(:)  ! field sizes
		INTEGER(SIZE_T),allocatable :: field_offset(:) ! field offset	
		character*(200) tablename
		INTEGER(HID_T), allocatable:: field_types(:) ! field types
		integer,allocatable::type_sizes(:)
	contains
		procedure::prepare_write_table
		procedure::prepare_read_table
		procedure::init_table
        procedure::write_column_real
        procedure::write_column_int
        procedure::write_column_str
		procedure::read_column_real
		procedure::read_column_int
		procedure::read_column_str
	end type

contains
	subroutine prepare_write_table(this, group_id)
		implicit none
		class(hdf5_table_type)::this
		INTEGER(HID_T) :: group_id      ! Group identifier
		!intger icheck	
		integer i
		this%group_id=group_id
		do i=1, this%nfields
			if(scan(this%field_names(i), ",").ne.0)then
				print*, "i, field names=",i, trim(adjustl(this%field_names(i)))
				print*, "field names should not contains ',', stoped"
				stop
			end if
		end do
		
		call prepare_table_write_hdf5(this%group_id,trim(adjustl(this%tablename)), &
		this%nrecords, this%nfields, this%field_types, this%field_sizes,this%type_sizes, this%field_names)
	end subroutine
	subroutine prepare_read_table(this, group_id)
		implicit none
		class(hdf5_table_type)::this
		INTEGER(HID_T) :: group_id      ! Group identifier
		!intger icheck	
		integer i
		this%group_id=group_id
		do i=1, this%nfields
			if(scan(this%field_names(i), ",").ne.0)then
				print*, "i, field names=",i, trim(adjustl(this%field_names(i)))
				print*, "field names should not contains ',', stoped"
				stop
			end if
		end do
		
		call prepare_table_info_hdf5(this%group_id,trim(adjustl(this%tablename)), &
		this%nrecords, this%nfields, this%field_types, this%field_sizes,this%type_sizes, this%field_names)
	end subroutine
	subroutine init_table(this,nfin, nrin,tname)
		implicit none
		integer nfin, nrin
		integer(hsize_t)::nf, nr
		character*(*) tname
		class(hdf5_table_type)::this
		nf=int(nfin,kind=hsize_t)
		nr=int(nrin,kind=hsize_t)
		if(allocated(this%field_names))then
			deallocate(this%field_names, this%field_sizes, this%field_offset, &
				this%field_types, this%type_sizes)
		end if
		allocate(this%field_names(nf),this%field_sizes(nf), this%field_offset(nf), &
				this%field_types(nf), this%type_sizes(nf))
		this%nfields=nf
		if(nr<=0)then
			print*, "init_table: error, number of records must be >=0",  " "//trim(adjustl(tname)), nrin
			stop
		end if
		this%nrecords=nr
		this%si=1
		this%type_sizes=1
		this%tablename=trim(adjustl(tname))
	end subroutine
	subroutine write_column_real(this, da)
		implicit none
		class(hdf5_table_type)::this
		real(8) da(this%nrecords)
		integer errcode
		INTEGER(HSIZE_T)   :: start = 0                        ! start record
		
		!print*, this%group_id
		!print*, trim(adjustl(this%tablename))
		!print*, this%field_names(this%si)
		!print*,	start
		!print*, this%nrecords
		!!print*, this%field_sizes(this%si)
		!print*, da(18000:18140)
		CALL h5tbwrite_field_name_f(this%group_id,this%tablename, this%field_names(this%si), &
		 start,this%nrecords,this%field_sizes(this%si), da,errcode)
		this%si=this%si+1
	end subroutine
	subroutine write_column_int(this, da)
		implicit none
		class(hdf5_table_type)::this
!		integer n
		integer da(this%nrecords)
		integer errcode
		INTEGER(HSIZE_T)   :: start = 0                        ! start record
		!print*, "this%si=", this%si
		CALL h5tbwrite_field_name_f(this%group_id,this%tablename,this%field_names(this%si), &
		 start,this%nrecords,this%field_sizes(this%si), da,errcode)
		this%si=this%si+1
	end subroutine
	subroutine write_column_str(this, da, str_size)
		implicit none
		class(hdf5_table_type)::this
		integer str_size
		character*(str_size) da(this%nrecords)
		integer errcode
		INTEGER(HSIZE_T)   :: start = 0                        ! start record
		!print*, "this%si=", this%si
		CALL h5tbwrite_field_name_f(this%group_id,this%tablename, this%field_names(this%si), &
		 start,this%nrecords,this%field_sizes(this%si)*str_size, da,errcode)
		this%si=this%si+1
	end subroutine

	subroutine read_column_int(this, da, field_name)
		implicit none
		class(hdf5_table_type)::this
		CHARACTER(LEN=*),optional :: field_name  ! field names
		INTEGER da(this%nrecords)
		integer errcode
		INTEGER(HSIZE_T)   :: start = 0  
		if(present(field_name))then
			call h5tbread_field_name_f(this%group_id,this%tablename,field_name,start,this%nrecords,this%field_sizes(this%si),&
			da,errcode)
		else
			call h5tbread_field_name_f(this%group_id,this%tablename,this%field_names(this%si),&
				start,this%nrecords,this%field_sizes(this%si),	da,errcode)
			this%si=this%si+1
		end if
	
	end subroutine
	subroutine read_column_real(this, da, field_name)
		implicit none
		class(hdf5_table_type)::this
		CHARACTER(LEN=*),optional :: field_name  ! field names
		real*8 da(this%nrecords)
		integer errcode
		INTEGER(HSIZE_T)   :: start = 0  
		if(present(field_name))then
			call h5tbread_field_name_f(this%group_id,this%tablename,field_name,start,this%nrecords,this%field_sizes(this%si),&
			da,errcode)
		else
			call h5tbread_field_name_f(this%group_id,this%tablename,this%field_names(this%si),&
				start,this%nrecords,this%field_sizes(this%si),	da,errcode)
			this%si=this%si+1
		end if
	
	end subroutine
	subroutine read_column_str(this, da, str_size, field_name)
		implicit none
		class(hdf5_table_type)::this
		CHARACTER(LEN=*),optional :: field_name  ! field names
		integer str_size
		character*(str_size) da(this%nrecords)
		integer errcode
		INTEGER(HSIZE_T)   :: start = 0  

		if(present(field_name))then
			call h5tbread_field_name_f(this%group_id,this%tablename,field_name,start,this%nrecords,this%field_sizes(this%si),&
			da,errcode)
		else
			call h5tbread_field_name_f(this%group_id,this%tablename,this%field_names(this%si),&
				start,this%nrecords,this%field_sizes(this%si)*str_size,	da,errcode)
			this%si=this%si+1
		end if
	
	end subroutine
end module
subroutine prepare_table_info_hdf5(group_id,tablename, nrecords, nfields, field_types, field_sizes,type_sizes, field_names)
	use hdf5
	use h5lt
	use h5tb
	implicit none
	character*(*) tablename
	INTEGER(HSIZE_T) :: nfields 
	INTEGER(HSIZE_T) :: nrecords 
	INTEGER(HID_T) :: group_id      ! Group identifier
	CHARACTER(LEN=*), DIMENSION(1:nfields) :: field_names  ! field names
	INTEGER(SIZE_T),  DIMENSION(1:nfields) :: field_sizes  ! field sizes
	INTEGER(SIZE_T),  DIMENSION(1:nfields) :: field_offset ! field offset	
	integer::type_sizes(nfields)
	INTEGER(SIZE_T)    :: offset                           ! Member's offset
	INTEGER(SIZE_T)  type_size
	INTEGER(HSIZE_T) :: chunk_size 
	INTEGER, PARAMETER :: compress = 0                     ! compress
	integer errcode,i
	INTEGER(HID_T),   DIMENSION(1:nfields) :: field_types  ! field types

	offset=0; type_size=0
	chunk_size=nrecords
	do i=1, nfields
		 field_offset(i) = offset
		 !field_sizes(i)=type_sizes(i)
		! print*, "i=",i
		 CALL h5tget_size_f(field_types(i),field_sizes(i) , errcode)
		  offset = offset + field_sizes(i)*type_sizes(i)
		  !print*, offset, field_sizes(i), type_sizes(i)
	end do
	type_size= offset
end subroutine

subroutine prepare_table_write_hdf5(group_id,tablename, nrecords, nfields, field_types, field_sizes,type_sizes, field_names)
	use hdf5
	use h5lt
	use h5tb
	implicit none
	character*(*) tablename
	INTEGER(HSIZE_T) :: nfields 
	INTEGER(HSIZE_T) :: nrecords 
	INTEGER(HID_T) :: group_id      ! Group identifier
	CHARACTER(LEN=*), DIMENSION(1:nfields) :: field_names  ! field names
	INTEGER(SIZE_T),  DIMENSION(1:nfields) :: field_sizes  ! field sizes
	INTEGER(SIZE_T),  DIMENSION(1:nfields) :: field_offset ! field offset	
	integer::type_sizes(nfields)
	INTEGER(SIZE_T)    :: offset                           ! Member's offset
	INTEGER(SIZE_T)  type_size
	INTEGER(HSIZE_T) :: chunk_size 
	INTEGER, PARAMETER :: compress = 0                     ! compress
	integer errcode,i
	INTEGER(HID_T),   DIMENSION(1:nfields) :: field_types  ! field types

	offset=0; type_size=0
	chunk_size=nrecords
	do i=1, nfields
		 field_offset(i) = offset
		 !field_sizes(i)=type_sizes(i)
		! print*, "i=",i
		 CALL h5tget_size_f(field_types(i),field_sizes(i) , errcode)
		  offset = offset + field_sizes(i)*type_sizes(i)
		  !print*, offset, field_sizes(i), type_sizes(i)
	end do

	type_size= offset
	!print*, "..."
	call h5tbmake_table_f(tablename, group_id, tablename, nfields,&
       nrecords,   type_size, field_names,  field_offset, field_types,&
       chunk_size,  compress, errcode )

end subroutine
