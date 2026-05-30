        subroutine READPAR_INT_AUTO(PAR,N, file_unit, comstr)
        	implicit none
			character*200 string, fpara, str
			character*(1) comstr
			integer n
			integer lim(2,100),nsub,i,j,error,k,par(N)
			integer file_unit
	
			!----------read in the number steps	
20      	read(file_unit,fmt='(a200)') string
			call mio_spl(200,string,nsub,lim) 
			if(string(lim(1,1):lim(1,1))==comstr)goto 20
			read(unit=string(lim(1,nsub):lim(2,nsub)),fmt=*) par
			return
		end subroutine
		
        subroutine READPAR_DBL_AUTO(PAR,N,file_unit,comstr)
            implicit none
			character*200 string, fpara, str
			character*1 comstr
			integer N
			integer lim(2,100),nsub,i,j,error,k
			integer file_unit,num
			real(8) par(N)

			!----------read in the number steps	
20      	read(file_unit,fmt='(a200)') string
			call mio_spl(200,string,nsub,lim) 
			if(string(lim(1,1):lim(1,1))==comstr)goto 20	
			do i=1, nsub
           ! 	print*,string(lim(1,i):lim(2,i))
				read(unit=string(lim(1,i):lim(2,i)),fmt=*) par(i)
			end do
			return
		end subroutine		

        subroutine READPAR_DBL_ONE(PAR,file_unit,comstr)
            implicit none
			character*200 string, fpara, str
			character*1 comstr
			integer lim(2,100),nsub,i,j,error,k
			integer file_unit,num
			real(8) par

			!----------read in the number steps	
20      	read(file_unit,fmt='(a200)') string
			call mio_spl(200,string,nsub,lim) 
			if(string(lim(1,1):lim(1,1))==comstr)goto 20			
			read(unit=string(lim(1,nsub):lim(2,nsub)),fmt=*) par
			return
		end subroutine		

        subroutine READPAR_STR(str,file_unit,comstr)
        	implicit none
			character*200 string, fpara,str
			character*1 comstr,sp
			integer lim(2,100),nsub,i,j,error,k
			integer file_unit,num

			!----------read in the number steps	
20      	read(file_unit,fmt='(a200)') string
			call mio_spl(200,string,nsub,lim) 
			if(string(lim(1,1):lim(1,1))==comstr)goto 20
			str=string(lim(1,nsub):lim(2,nsub))
			return
		end subroutine	

		subroutine READPAR_STR_SPLIT_TWO(str,file_unit,comstr,ier)
        	implicit none
			character(len=200) string, strbg, stred,str(2)
			character*1 comstr,sp
			integer file_unit, j, ispt, ib, id,ier, flag
			!----------read in the number steps	
20      	read(file_unit,fmt='(a200)',err=100,iostat=flag) string
			!print*, "flag=",flag, string
			
			ib=0; id=0
			!str_(1:200)=string
			ispt=index(string,"=")
			strbg=string(1:ispt-1)
			stred=string(ispt+1:)
			sp=strbg(1:1)
			!print*,"str=",sp, sp.eq."#", ispt
			if((ispt.eq.0.or.sp.eq."#").and.flag==0)then 
				!print*, "sp=#"!,flag
				goto 20
			end if
			str(1)=strbg
			str(2)=stred
			ier=0
			if(flag.ne.0)then
100				ier=-1
			end if
		end subroutine	

        subroutine READPAR_INT(PAR,file_unit,comstr)
        	implicit none
			character*200 string, fpara, str
			integer lim(2,100),nsub,i,j,error,k,par
			integer file_unit
				character*1 comstr,sp
			!----------read in the number steps	
20      	read(file_unit,fmt='(a200)') string
			call mio_spl(200,string,nsub,lim) 
			if(string(lim(1,1):lim(1,1))==comstr)goto 20
			read(unit=string(lim(1,nsub):lim(2,nsub)),fmt=*) par
			return
		end subroutine

        subroutine READPAR_LOG(PAR,file_unit,comstr)
        	implicit none
			character*200 string, fpara, str
			integer lim(2,100),nsub,i,j,error,k
			logical par
			integer file_unit
				character*1 comstr,sp
			!----------read in the number steps	
20      	read(file_unit,fmt='(a200)') string
			call mio_spl(200,string,nsub,lim) 
			if(string(lim(1,1):lim(1,1))==comstr)goto 20
			read(unit=string(lim(1,nsub):lim(2,nsub)),fmt=*) par
			return
		end subroutine
		
        subroutine READPAR_DBL(PAR,file_unit,comstr)
            implicit none
			character*200 string, fpara
			integer lim(2,100),i,j,error,k
			integer file_unit,num,nsub
			real(8) par
			character*1 comstr
			!----------read in the number steps	
20      	read(file_unit,fmt='(a200)') string
			call mio_spl(200,string,nsub,lim)  
			if(string(lim(1,1):lim(1,1))==comstr)goto 20
			read(unit=string(lim(1,nsub):lim(2,nsub)),fmt=*) par
			return
		end subroutine			
 
!!!-------------------------------------------------------------------

       subroutine READPAR_INT_AUTO_SP(PAR,N, file_unit, comstr,sp)
        	implicit none
			character*200 string, fpara, str
			character*(1) comstr,sp
			integer n
			integer lim(2,100),nsub,i,j,error,k,par(N)
			integer file_unit
	
			!----------read in the number steps	
20      	read(file_unit,fmt='(a200)') string
			call mio_spl_sp(200,string,nsub,lim,sp) 
			if(string(lim(1,1):lim(1,1))==comstr)goto 20
			read(unit=string(lim(1,nsub):lim(2,nsub)),fmt=*) par
			return
		end subroutine
		
        subroutine READPAR_DBL_AUTO_SP(PAR,N,file_unit,comstr,sp)
			character*200 string, fpara, str
			character*1 comstr,sp
			integer N
			integer lim(2,100),nsub,i,j,error,k
			integer file_unit,num
			real(8) par(N)

			!----------read in the number steps	
20      	read(file_unit,fmt='(a200)') string
			call mio_spl_sp(200,string,nsub,lim,sp) 
			if(string(lim(1,1):lim(1,1))==comstr)goto 20	
           ! print*,string(lim(1,nsub):lim(2,nsub))
			read(unit=string(lim(1,nsub):lim(2,nsub)),fmt=*) par
			return
		end subroutine		

        subroutine READPAR_DBL_ONE_SP(PAR,file_unit,comstr,sp)
			character*200 string, fpara, str
			character*1 comstr,sp
			integer lim(2,100),nsub,i,j,error,k
			integer file_unit,num
			real(8) par

			!----------read in the number steps	
20      	read(file_unit,fmt='(a200)') string
			call mio_spl_sp(200,string,nsub,lim,sp) 
			if(string(lim(1,1):lim(1,1))==comstr)goto 20			
			read(unit=string(lim(1,nsub):lim(2,nsub)),fmt=*) par
			return
		end subroutine		

        subroutine READPAR_STR_SP(str,file_unit,comstr,sp)
        	implicit none
			character*200 string, fpara, str
			character*1 comstr,sp
			integer lim(2,100),nsub,i,j,error,k
			integer file_unit,num

			!----------read in the number steps	
20      	read(file_unit,fmt='(a200)') string
			call mio_spl_sp(200,string,nsub,lim,sp) 
			if(string(lim(1,1):lim(1,1))==comstr)goto 20
			str=string(lim(1,nsub):lim(2,nsub))
			return
		end subroutine			

        subroutine READPAR_STR_AUTO_SP(strall,str,N,Nnum,file_unit,comstr,sp,sepstr)
        	implicit none
            integer n, nnum
			character*200 string, fpara, strall,str(N)
			character*1 comstr,sepstr,sp
			integer lim(2,100),nsub,i,j,error,k
			integer file_unit,num

			!----------read in the number steps	
20      	read(file_unit,fmt='(a200)') string
			call mio_spl_sp(200,string,nsub,lim,sp)  
			if(string(lim(1,1):lim(1,1))==comstr)goto 20
			strall=string(lim(1,nsub):lim(2,nsub))
            call mio_spl_sp(200,strall,nnum,lim,sepstr)
            do i=1, nnum
                str(i)=strall(lim(1,i):lim(2,i))
            end do
			return  
		end subroutine		
        subroutine READPAR_INT_SP(PAR,file_unit,comstr,sp,ier)
        	implicit none
			character*200 string, fpara, str
			integer lim(2,100),nsub,i,j,error,k,par
			integer file_unit, ier
            character*1 comstr,sp
			!----------read in the number steps	
20      	read(file_unit,fmt='(a200)') string
			call mio_spl_sp(200,string,nsub,lim,sp) 
			if(string(lim(1,1):lim(1,1))==comstr)goto 20
			read(unit=string(lim(1,nsub):lim(2,nsub)),fmt=*,err=100) par
            ier=0
			return
100         ier=-1           
		end subroutine

        subroutine READPAR_LOG_SP(PAR,file_unit,comstr,sp)
        	implicit none
			character*200 string, fpara, str
			integer lim(2,100),nsub,i,j,error,k
			logical par
			integer file_unit
				character*1 comstr,sp
			!----------read in the number steps	
20      	read(file_unit,fmt='(a200)') string
			call mio_spl_sp(200,string,nsub,lim,sp) 
			if(string(lim(1,1):lim(1,1))==comstr)goto 20
			read(unit=string(lim(1,nsub):lim(2,nsub)),fmt=*) par
			return
		end subroutine
		
        subroutine READPAR_DBL_SP(PAR,file_unit,comstr,sp,ier)
			implicit none
			character*200 string
			integer lim(2,100),nsub
			integer file_unit,ier
			real(8) par
			character*1 comstr,sp
			!----------read in the number steps	
20      	read(file_unit,fmt='(a200)') string
			call mio_spl_sp(200,string,nsub,lim,sp)  
			if(string(lim(1,1):lim(1,1))==comstr)goto 20
			read(unit=string(lim(1,nsub):lim(2,nsub)),fmt=*,ERR=100) par
			IER=0
			RETURN
100         IER=-1			
		end subroutine			
		subroutine skip_comments(file_unit,comstr)
			implicit none
			integer file_unit
			character*200 string
			character*1 comstr
20      	read(file_unit,fmt='(a200)') string			
			if(string(1:1)==comstr)goto 20
			backspace(file_unit)
			return
		end subroutine

		subroutine number_formatter(d,num_fix, str, num_str)  ! convert double precision into a number with X digits
			implicit none
			real(8) d
			integer num_fix, num_str
			character*(num_str) str
			character*200 formatter
			integer num_f, num_b, digits, num_char_tot
			character*(10) str_num_tot, str_num_back, str_zeros
			if (num_str<num_fix+1) then
				print*, "increase the size of num_str! stoped!"
				stop
			end if
			if(d<0)then
				print*, "d should be larger than zero! stoped!"
			end if
			if(d.eq.0d0)then
				num_f=1; num_b=num_fix-1
				num_char_tot=num_fix+1
				write(unit=str_num_tot,fmt="(I10)") num_char_tot
				formatter="(F"//trim(adjustl(str_num_tot))
				formatter=trim(adjustl(formatter))//"."
				write(unit=str_num_back,fmt="(I10)") num_b
				formatter=trim(adjustl(formatter))//trim(adjustl(str_num_back))//")"
				write(unit=str, fmt=formatter) d
				return
			end if
			
			if(d.ge.1)then
				digits=int(log10(d))+1
				if(num_fix>digits)then
					num_f=digits
					num_b=num_fix-digits
					num_char_tot=num_f+num_b+2
					write(unit=str_num_tot,fmt="(I10)") num_char_tot
					formatter="(F"//trim(adjustl(str_num_tot))
					formatter=trim(adjustl(formatter))//"."
					write(unit=str_num_back,fmt="(I10)") num_b
					formatter=trim(adjustl(formatter))//trim(adjustl(str_num_back))//")"
					write(unit=str, fmt=formatter) d
				elseif(num_fix.eq.digits)then
					!num_f=
					num_char_tot=digits+1
					write(unit=str_num_tot,fmt="(I10)") num_char_tot
					formatter="(I"//trim(adjustl(str_num_tot))//")"
					write(unit=str, fmt=formatter) nint(d)
				else
					num_char_tot=digits+1
					write(unit=str_num_tot,fmt="(I10)") num_char_tot
					formatter="(I"//trim(adjustl(str_num_tot))//")"
					write(unit=str, fmt=formatter) int(nint(d/10**(digits-num_fix))*10**(digits-num_fix))
				end if
				!print*, d, digits, num_fix, num_char_tot
			else
				digits=abs(int(log10(d))-1)
				if(num_fix>digits)then
					num_f=digits
					num_b=num_fix-digits
					num_char_tot=num_f+num_b+1
					write(unit=str_num_tot,fmt="(I10)") num_char_tot
					formatter="(F"//trim(adjustl(str_num_tot))
					formatter=trim(adjustl(formatter))//"."
					write(unit=str_num_back,fmt="(I10)") num_b
					formatter=trim(adjustl(formatter))//trim(adjustl(str_num_back))//")"
					write(unit=str, fmt=formatter) d
				elseif(num_fix.eq.digits)then
					!num_f=
					num_char_tot=digits
					write(unit=str_num_tot,fmt="(I10)") num_char_tot
					formatter="(I"//trim(adjustl(str_num_tot))//")"
					write(unit=str, fmt=formatter) nint(d)
				else
					num_char_tot=digits
					write(unit=str_num_tot,fmt="(I10)") num_char_tot
					formatter="(I"//trim(adjustl(str_num_tot))//")"
					write(unit=str, fmt=formatter) int(nint(d/10**(digits-num_fix))*10**(digits-num_fix))
					!str_zeros=repeat("0", digits-num_fix)

				end if
			end if
			
			
		end subroutine
		subroutine remove_tab_in_string(str_in,str_out)
			implicit none
			character*(*),intent(in):: str_in
			character*(*),intent(out)::str_out
			integer i,j,length
			length=len(str_in)
			
			j=0
			str_out=''
			do i=1, len(str_in)
				if(iachar(str_in(i:i)).ne.9 )then
					j=j+1
					str_out(j:j)=str_in(i:i)
				end if
			end do
			str_out=str_out(1:min(j,len(str_out)))
		end subroutine
