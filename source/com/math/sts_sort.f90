      SUBROUTINE sort(n,arr)
      INTEGER n,M,NSTACK
      REAL(8) arr(n)
      PARAMETER (M=7,NSTACK=50)
      INTEGER i,ir,j,jstack,k,l,istack(NSTACK)
      REAL(8) a,temp
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 12 j=l+1,ir
          a=arr(j)
          do 11 i=j-1,1,-1
            if(arr(i).le.a)goto 2
            arr(i+1)=arr(i)
11        continue
          i=0
2         arr(i+1)=a
12      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        temp=arr(k)
        arr(k)=arr(l+1)
        arr(l+1)=temp
        if(arr(l+1).gt.arr(ir))then
          temp=arr(l+1)
          arr(l+1)=arr(ir)
          arr(ir)=temp
        endif
        if(arr(l).gt.arr(ir))then
          temp=arr(l)
          arr(l)=arr(ir)
          arr(ir)=temp
        endif
        if(arr(l+1).gt.arr(l))then
          temp=arr(l+1)
          arr(l+1)=arr(l)
          arr(l)=temp
        endif
        i=l+1
        j=ir
        a=arr(l)
3       continue
          i=i+1
        if(arr(i).lt.a)goto 3
4       continue
          j=j-1
        if(arr(j).gt.a)goto 4
        if(j.lt.i)goto 5
        temp=arr(i)
        arr(i)=arr(j)
        arr(j)=temp
        goto 3
5       arr(l)=arr(j)
        arr(j)=a
        jstack=jstack+2
        if(jstack.gt.NSTACK)pause 'NSTACK too small in sort'
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
      END
!  (C) Copr. 1986-92 Numerical Recipes Software :)z%+.
      SUBROUTINE sort2(n,arr,brr)
      INTEGER n,M,NSTACK
      REAL(8) arr(n),brr(n)
      PARAMETER (M=7,NSTACK=50)
      INTEGER i,ir,j,jstack,k,l,istack(NSTACK)
      REAL(8) a,b,temp
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 12 j=l+1,ir
          a=arr(j)
          b=brr(j)
          do 11 i=j-1,l,-1
            if(arr(i).le.a)goto 2
            arr(i+1)=arr(i)
            brr(i+1)=brr(i)
11        continue
          i=l-1
2         arr(i+1)=a
          brr(i+1)=b
12      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        temp=arr(k)
        arr(k)=arr(l+1)
        arr(l+1)=temp
        temp=brr(k)
        brr(k)=brr(l+1)
        brr(l+1)=temp
        if(arr(l).gt.arr(ir))then
          temp=arr(l)
          arr(l)=arr(ir)
          arr(ir)=temp
          temp=brr(l)
          brr(l)=brr(ir)
          brr(ir)=temp
        endif
        if(arr(l+1).gt.arr(ir))then
          temp=arr(l+1)
          arr(l+1)=arr(ir)
          arr(ir)=temp
          temp=brr(l+1)
          brr(l+1)=brr(ir)
          brr(ir)=temp
        endif
        if(arr(l).gt.arr(l+1))then
          temp=arr(l)
          arr(l)=arr(l+1)
          arr(l+1)=temp
          temp=brr(l)
          brr(l)=brr(l+1)
          brr(l+1)=temp
        endif
        i=l+1
        j=ir
        a=arr(l+1)
        b=brr(l+1)
3       continue
          i=i+1
        if(arr(i).lt.a)goto 3
4       continue
          j=j-1
        if(arr(j).gt.a)goto 4
        if(j.lt.i)goto 5
        temp=arr(i)
        arr(i)=arr(j)
        arr(j)=temp
        temp=brr(i)
        brr(i)=brr(j)
        brr(j)=temp
        goto 3
5       arr(l+1)=arr(j)
        arr(j)=a
        brr(l+1)=brr(j)
        brr(j)=b
        jstack=jstack+2
        if(jstack.gt.NSTACK)pause 'NSTACK too small in sort2'
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
      END

      SUBROUTINE sort3_i(n,ra,rb,rc,wksp,iwksp)
      INTEGER n,iwksp(n)
      REAL(8) ra(n),rb(n),rc(n),wksp(n)
	  external:: indexx
      INTEGER j
     ! call indexx(n,ra,iwksp)
      print*, "indexx may have some problem"
      stop
      do 11 j=1,n
        wksp(j)=ra(j)
11    continue
      do 12 j=1,n
        ra(j)=wksp(iwksp(j))
12    continue
      do 13 j=1,n
        wksp(j)=rb(j)
13    continue
      do 14 j=1,n
        rb(j)=wksp(iwksp(j))
14    continue
      do 15 j=1,n
        wksp(j)=rc(j)
15    continue
      do 16 j=1,n
        rc(j)=wksp(iwksp(j))
16    continue
      return
      END
!  (C) Copr. 1986-92 Numerical Recipes Software :)z%+.

SUBROUTINE sort3(n,ra,rb,rc)
	integer n
	real(8) ra(n),rb(n),rc(n)
	real(8),allocatable:: wksp(:)
	integer,allocatable::iwksp(:)
	allocate(wksp(n))
	allocate(iwksp(n))
	call sort3_i(n,ra,rb,rc, wksp, iwksp)
end subroutine
subroutine sort_arr_index(n,a,iwksp)
	integer n,i
	real(8) wksp(n)
	real(8) a(n)
	integer iwksp(n)
	
	!allocate(wksp(n))
	do i=1, n
		wksp(i)=a(i)
	end do
	do i=1, n
		a(i)=wksp(iwksp(i))
	end do
end subroutine

SUBROUTINE sortn(n,m,arr)
	integer n,i,j
	real(8) arr(n,m)
	real(8),allocatable:: wksp(:)
	integer,allocatable::iwksp(:)
	allocate(wksp(n))
	allocate(iwksp(n))
	call sort3_i(n,arr(:,1),arr(:,2),arr(:,3), wksp, iwksp)
	do i=4, m
		do j=1,n
		wksp(j)=arr(j,i)
		end do
		do j=1,n
		arr(j,i)=wksp(iwksp(j))
		end do
	end do
end subroutine

