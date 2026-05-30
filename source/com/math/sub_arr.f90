
subroutine find_arr_closest(arr, n, val, idx_closest)  ! find the position of an array element that is closest to the given value
	implicit none
	integer n, idx_closest
	real(8) arr(n), val
	integer i, j, mid
	! arr has to be sorted small to big 1 to n
	if(val<=arr(1))then
		idx_closest=1
		return
	end if
	if(val>=arr(n))then
		idx_closest=n
        return
	end if
	
    if(n.eq.2)then
        if(val-arr(1)>arr(2)-val)then
            idx_closest=2
        else
            idx_closest=1
        end if
        return
    endif
    i=1; j=n; mid=0
	do while(i<j)
		mid=(i+j)/2
		if(val.eq.arr(mid))then
			idx_closest=mid
			return
		end if
		if(val<arr(mid))then
			if(mid>=2.and.val>arr(mid-1))then
				if(val-arr(mid-1)>=arr(mid)-val)then
					idx_closest=mid
				else
					idx_closest=mid-1
				end if
				return
			end if
			j=mid
		else
			if(mid<=n-1.and.val<arr(mid+1))then
				if(arr(mid+1)-val>val-arr(mid))then
					idx_closest=mid
				else
					idx_closest=mid+1
				end if
				return
			endif
			i=mid+1
		end if
	end do
end subroutine


subroutine arrmax(arr,n, max)
	integer i,n
	real(8) max,arr(n)
	max=arr(1)
	do i=2,n
		if(max<arr(i)) max=arr(i)
	end do
end subroutine
subroutine arrmin(arr,n,min)
	integer i,n
	real(8) min,arr(n)
	min=arr(1)
	do i=2,n
		if(min>arr(i)) min=arr(i)
	end do
end subroutine
subroutine arravg(arr,n,avg)
	integer n
	real(8) avg,arr(n),sum

	avg=sum(arr(1:n))/real(n)

end subroutine
subroutine arrsum_weight(arr,weight, n,sums)
	integer n
	real(8) sums,arr(n)
	real(8) weight(n)
	sums=sum(arr(1:n)*weight(1:n))
end subroutine
subroutine arravg_weight(arr,weight, n,avg)
	integer n
	real(8) avg,arr(n)
	real(8) weight(n)
	avg=sum(arr(1:n)*weight(1:n))/sum(weight(1:n))
end subroutine
subroutine get_val_percent_of_arr_weight(arr, w, n, value, percent)
	implicit none
	integer n, i
	real(8) nc
	real(8) arr(n), value, percent(2)   !percenter(1):percentage of number of members in array less than value
										!percenter(2):percentage of number of members in array larger than value
	real(8) w(n), wtot
	nc=0
	wtot=0
	do i=1, n
		if(arr(i)<value)then
			nc=nc+w(i)
		end if
	end do
	wtot=sum(w(1:n))
	percent(1)=nc/wtot
	percent(2)=1-percent(1)
end subroutine
subroutine arrsct(arr,n,mean,sct)
	implicit none
	integer n,i
	real(8)	arr(n),mean,sct
	sct=0
	do i=1,n
		sct=sct+(arr(i)-mean)**2
	end do
	sct=sqrt(sct/real(n-1))

end subroutine

subroutine arrmedian(arr,n,median)
	implicit none
	integer n,i
	real(8) arr(n),median
	real(8),allocatable:: brr(:)
	allocate(brr(n))
	brr=arr
	call sort(n,brr)
	if(n>1)then
		median=brr(int(n/2))
	else
		median=brr(n)
	endif
end subroutine

subroutine arrpos(arr,n,px,x)
	!return the number in an array where it (closest) 
	!corresponds to px% in position
	implicit none
	integer n,i
	real(8) arr(n),px,x
	real(8),allocatable:: brr(:)
	allocate(brr(n))
	brr=arr
	call sort(n,brr)
	if(n>1)then
		x=brr(int(px*n))
	else
		x=brr(n)
	endif
end subroutine

subroutine arr2_remove_zeros(x,y,n, xmd,ymd, nout)
	implicit none
	integer n,i,nout,j
	real(8) x(n),y(n),xmd(n),ymd(n)
	j=0
	do i=1, n
		if(y(i).ne.0d0)then
			j=j+1
			xmd(j)=x(i); ymd(j)=y(i)
		end if
	end do
	nout=j
end subroutine

subroutine arr2_remove_zeros_replace(x,y,n)
	implicit none
	integer n,i,nout,j
	real(8) x(n),y(n)
	real(8),allocatable::xmd(:),ymd(:)
	allocate(xmd(n),ymd(n))
	xmd=0;ymd=0
	call arr2_remove_zeros(x,y, n, xmd, ymd, nout)
	x=xmd
	y=ymd
	n=nout
end subroutine

subroutine print_array_real(a, n)
	implicit none
	integer n, i
	real(8) a
	write(unit=*, fmt="(2A20)") "N", "X"
	do i=1, n
		write(unit=*, fmt="(2A20)") i, a(i)
	end do
end subroutine