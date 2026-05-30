! Returns the inverse of a matrix calculated by finding the LU
! decomposition.  Depends on LAPACK.
function Matrix_inv(A, nx) 
    integer nx
    real(8), intent(in) :: A(nx,nx)
    real(8) Matrix_inv(nx,nx)
  
    real(8) work(nx) ! work array for LAPACK
    integer ipiv(nx)  ! pivot indices
    integer :: n, info
  
    ! External procedures defined in LAPACK
    external DGETRF
    external DGETRI
  
    ! Store A in Ainv to prevent it from being overwritten by LAPACK
    Matrix_inv = A
    n = nx
   ! print*, "in matrix_inv 1"
    ! DGETRF computes an LU factorization of a general M-by-N matrix A
    ! using partial pivoting with row interchanges.
    call DGETRF(n, n, Matrix_inv, n, ipiv, info)
   ! print*, "in matrix_inv 2"
    if (info /= 0) then
       stop 'Matrix is numerically singular!'
    end if
  
    ! DGETRI computes the inverse of a matrix using the LU factorization
    ! computed by DGETRF.
    call DGETRI(n, Matrix_inv, n, ipiv, work, n, info)
  
    if (info /= 0) then
       stop 'Matrix inversion failed!'
    end if
  end function Matrix_inv