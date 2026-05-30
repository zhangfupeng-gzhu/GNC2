!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MIO_SPL.FOR    (ErikSoft  14 November 1999)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Given a character string STRING, of length LEN bytes, the routine finds 
! the beginnings and ends of NSUB substrings present in the original, and 
! delimited by spaces. The positions of the extremes of each substring are 
! returned in the array DELIMIT.
! Substrings are those which are separated by spaces or the = symbol.
!
!------------------------------------------------------------------------------
!
      subroutine mio_spl (len,str,nsub,delimit)
      !subroutine mio_spl (str,nsub)
      implicit none
      character*(1) str(len)
      integer delimit(2,100)
      integer nsub,len
      integer j,k
      character*1 c
!
!------------------------------------------------------------------------------
!
      nsub = 0
      j = 0
      c = ' '
      delimit(1,1) = -1
!
! Find the start of string
10    j = j + 1
      if (j.gt.len) goto 99
      c = str(j)
      if (c.eq.' '.or.c.eq.'=') goto 10
!
! Find the end of string
      k = j
20    k = k + 1
      if (k.gt.len) goto 30
      c = str(k)
      if (c.ne.' '.and.c.ne.'=') goto 20
!
! Store details for this string
30    nsub = nsub + 1
      delimit(1,nsub) = j
      delimit(2,nsub) = k - 1
!
      if (k.lt.len) then
        j = k
        goto 10
      end if
!
  99  continue
!
!------------------------------------------------------------------------------
!
      return
      end
!

!------------------------------------------------------------------------------
!
      subroutine mio_spl_sp(len,string,nsub,delimit,separator)
!
      implicit none
!
! Input/Output
      integer len,nsub,delimit(2,100)
      character*1 string(len)
      character*1 separator
!
! Local
      integer j,k
      character*1 c
!
!------------------------------------------------------------------------------
!
      nsub = 0
      j = 0
      c = ' '
      delimit(1,1) = -1
!
! Find the start of string
  10  j = j + 1
      if (j.gt.len) goto 99
      c = string(j)
      if (c.eq.separator) goto 10
!
! Find the end of string
      k = j
  20  k = k + 1
      if (k.gt.len) goto 30
      c = string(k)
      if (c.ne.separator) goto 20
!
! Store details for this string
  30  nsub = nsub + 1
      delimit(1,nsub) = j
      delimit(2,nsub) = k - 1
!
      if (k.lt.len) then
        j = k
        goto 10
      end if
!
  99  continue
!
!------------------------------------------------------------------------------
!
      return
      end
!
