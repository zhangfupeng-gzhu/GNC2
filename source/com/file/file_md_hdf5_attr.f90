module md_hdf5_attr
    use hdf5
	use h5lt
	use h5tb
    implicit none
    !type hdf5_attr
    !    integer(HID_T) id
    !end type
contains 
    subroutine add_attr_int(group_id,attr_id, attr_name, attr_value)
        implicit none
        integer(HID_T)::group_id, aspace_id,attr_id
        character*(*) attr_name
        integer attr_value, error
        INTEGER(HSIZE_T), DIMENSION(1) :: adims = (/1/) ! Attribute dimension
        INTEGER(HSIZE_T), DIMENSION(1) :: data_dims=(/1/)
        CALL h5screate_simple_f(1, adims, aspace_id, error)
        CALL h5acreate_f(group_id, attr_name, H5T_NATIVE_INTEGER, aspace_id, attr_id, error)
        CALL h5awrite_f(attr_id, H5T_NATIVE_INTEGER, attr_value, data_dims, error)
        CALL h5aclose_f(attr_id, error)
        CALL h5sclose_f(aspace_id, error)
    end subroutine    
    subroutine add_attr_dble(group_id,attr_id, attr_name, attr_value)
        implicit none
        integer(HID_T)::group_id, aspace_id,attr_id
        character*(*) attr_name
        integer  error
        real(8) attr_value
        INTEGER(HSIZE_T), DIMENSION(1) :: adims = (/1/) ! Attribute dimension
        INTEGER(HSIZE_T), DIMENSION(1) :: data_dims=(/1/)
        CALL h5screate_simple_f(1, adims, aspace_id, error)
        CALL h5acreate_f(group_id, attr_name, H5T_NATIVE_DOUBLE, aspace_id, attr_id, error)
        CALL h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, attr_value, data_dims, error)
        CALL h5aclose_f(attr_id, error)
        CALL h5sclose_f(aspace_id, error)
    end subroutine   
    
    subroutine add_attr_int_arr(group_id,attr_id, attr_name, attr_value, n_attr)
        implicit none
        integer(HID_T)::group_id, aspace_id,attr_id
        character*(*) attr_name
        integer  error
        integer n_attr
        integer attr_value(n_attr)
        INTEGER(HSIZE_T), DIMENSION(1) :: adims  ! Attribute dimension
        INTEGER(HSIZE_T), DIMENSION(1) :: data_dims
        adims=n_attr
        data_dims=n_attr

        CALL h5screate_simple_f(1, adims, aspace_id, error)
        CALL h5acreate_f(group_id, attr_name, H5T_NATIVE_INTEGER, aspace_id, attr_id, error)
        CALL h5awrite_f(attr_id, H5T_NATIVE_INTEGER, attr_value, data_dims, error)
        CALL h5aclose_f(attr_id, error)
        CALL h5sclose_f(aspace_id, error)
    end subroutine  
    subroutine add_attr_dble_arr(group_id,attr_id, attr_name, attr_value, n_attr)
        implicit none
        integer(HID_T)::group_id, aspace_id,attr_id
        character*(*) attr_name
        integer  error
        integer n_attr
        real(8) attr_value(n_attr)
        INTEGER(HSIZE_T), DIMENSION(1) :: adims  ! Attribute dimension
        INTEGER(HSIZE_T), DIMENSION(1) :: data_dims
        adims=n_attr
        data_dims=n_attr

        CALL h5screate_simple_f(1, adims, aspace_id, error)
        CALL h5acreate_f(group_id, attr_name, H5T_NATIVE_DOUBLE, aspace_id, attr_id, error)
        CALL h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, attr_value, data_dims, error)
        CALL h5aclose_f(attr_id, error)
        CALL h5sclose_f(aspace_id, error)
    end subroutine       
    subroutine read_attr_dble_arr(file_id,dir_name,attr_id, attr_name, attr_value, n_attr)
        implicit none
        integer(HID_T)::file_id, aspace_id,attr_id
        character*(*) attr_name,dir_name
        integer  error
        integer n_attr
        real(8) attr_value(n_attr)
        INTEGER(HSIZE_T), DIMENSION(1) :: adims  ! Attribute dimension
        INTEGER(HSIZE_T), DIMENSION(1) :: data_dims
        adims=n_attr
        data_dims=n_attr

        call h5aopen_by_name_f(file_id,dir_name,attr_name,attr_id,error)
        CALL h5aread_f(attr_id, H5T_NATIVE_DOUBLE, attr_value, data_dims, error)
        CALL h5aclose_f(attr_id, error)
    end subroutine    

    subroutine read_attr_dble(file_id,dir_name,attr_id, attr_name, attr_value)
        implicit none
        integer(HID_T)::file_id, aspace_id,attr_id
        character*(*) attr_name,dir_name
        integer  error
        real(8) attr_value
        INTEGER(HSIZE_T), DIMENSION(1) :: adims = (/1/) ! Attribute dimension
        INTEGER(HSIZE_T), DIMENSION(1) :: data_dims=(/1/)
        call h5aopen_by_name_f(file_id,dir_name,attr_name,attr_id,error)
        CALL h5aread_f(attr_id, H5T_NATIVE_DOUBLE, attr_value, data_dims, error)
        CALL h5aclose_f(attr_id, error)
    end subroutine        
end module