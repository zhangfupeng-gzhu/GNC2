module fun_math
implicit none
contains
	function fvector_x(v1,v2)
			real(8) v1(3),v2(3),fvector_x(3)
			fvector_x(1)=v1(2)*v2(3)-v1(3)*v2(2)
			fvector_x(2)=v1(3)*v2(1)-v1(1)*v2(3)
			fvector_x(3)=v1(1)*v2(2)-v1(2)*v2(1)
	end function
	
	function fvector_mag(v)
			real(8) v(3),fvector_mag
			fvector_mag=sqrt(v(1)**2+v(2)**2+v(3)**2)	
	end function
	
	function fvector_unit(v)
		real(8) fvector_unit(3),v(3)
		fvector_unit=v/fvector_mag(v)
	end function

	function fvector_dot(v1,v2)
			real(8) fvector_dot,v1(3),v2(3)
			fvector_dot=v1(1)*v2(1)+v1(2)*v2(2)+v1(3)*v2(3)
	end function

    function fmatrix_33_33(a,b)
		real(8) a(3,3),b(3,3)
		real(8) fMatrix_33_33(3,3)
		integer i,j,k
		fMatrix_33_33=0
		do i=1,3
			do j=1,3
				do k=1,3
				fMatrix_33_33(i,j)=fMatrix_33_33(i,j)+a(i,k)*b(k,j)
				end do
			end do
		end do
	end function

	function fmatrix_13_33(a,b)
		real(8) a(3),b(3,3)
		real(8) fMatrix_13_33(3)
		integer i,j
		fMatrix_13_33=0
		do i=1,3
			do j=1,3
			fMatrix_13_33(i)=fMatrix_13_33(i)+a(j)*b(j,i)
			end do
		end do		
	end function
	
	function farrmax(arr,n)
		integer i,n
		real(8) farrmax,arr(n)
		farrmax=arr(1)
		do i=2,n
			if(farrmax<arr(i)) farrmax=arr(i)
		end do
	end function
	
	function farrmin(arr,n)
		integer i,n
		real(8) farrmin,arr(n)
		farrmin=arr(1)
		do i=2,n
			if(farrmin>arr(i)) farrmin=arr(i)
		end do
	end function
	
	real(8) function farravg(arr,n)
		integer n
		real(8) arr(n)
		farravg=sum(arr)/dble(n)
	end function
	real(8) function farrsct(arr, n)
		integer n
		real(8) arr(n),mean
		mean=farravg(arr,n)
		farrsct=(sum((arr-mean)**2)/dble(n-1))**0.5d0
	end function
	real(8) function farrsct2(arr,mean, n)
		integer n
		real(8) arr(n), mean
		farrsct2=(sum((arr-mean)**2)/dble(n-1))**0.5d0
	end function
	
!	function fgaussian()
		
!	end function
end module
