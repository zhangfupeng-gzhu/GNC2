

real(8) function get_clone_deep(en, log_bd_sep, en0)
    implicit none
    real(8) en, en0,log_bd_sep
    get_clone_deep=log(en/en0)/log_bd_sep+1
    if(get_clone_deep<0)then
        get_clone_deep=0
    end if
end function

logical function Istransit(en0,en1,i,nl)
	use com_main_gw
    implicit none
	real(8) en0, en1, en0p, en1p,get_clone_deep
	integer i,jk,leveln,nl

	i=-1
	Istransit=.false.
!	level0=log10(en0);level1=log10(en1)
    
	leveln=log10clone_emax
	!print*, leveln
	!stop
    en0p=get_clone_deep(en0, log_clone_bd_sep, ctl%clone_e0)
    en1p=get_clone_deep(en1, log_clone_bd_sep, ctl%clone_e0)
    !if(ctl%chattery.ge.3) then
    !    print*, "clone_emax=",clone_emax
    !    print*, "en0p, en1p=", en0p, en1p, leveln
    !end if
	do jk=1, leveln
    !    if(ctl%chattery.ge.3) print*, "jk=",jk 
		if(en0p<jk.and.en1p>jk)then
			Istransit=.true.
			i=jk
            nl=int(abs(en1p-en0p))+1
			return 
		end if
	end do
    !if(ctl%chattery.ge.3)then
    !    print*, "istransit=",istransit
    !endif
!	print*, "error in istransit",en0,en1
!	stop
end function

logical function Invtransit(en0,en1,i, nl)
	use com_main_gw
    implicit none
	real(8) en0, en1, en0p, en1p, get_clone_deep
	integer i,jk,leveln,nl

	Invtransit=.false.
	i=-1

	leveln=log10clone_emax
    en0p=get_clone_deep(en0, log_clone_bd_sep, ctl%clone_e0)
    en1p=get_clone_deep(en1, log_clone_bd_sep,  ctl%clone_e0)
    if(ctl%chattery.ge.5)then
        print*, "en0p, en1p, clone_e0, clone_emax=", en0p, en1p, ctl%clone_e0, clone_emax
    end if
	do jk=1, leveln
        !if(ctl%chattery.ge.3)then
        !    print*, "jk=", jk, en0p>10**jk, en1p<10**jk
        !end if
		if(en0p>jk.and.en1p<jk)then
			Invtransit=.true.
			i=jk
            nl=int(abs(en1p-en0p))+1
			return 
		end if
	end do
    
!	print*, "error in invtransit",en0,en1
!	stop
end function



subroutine clone_scheme(pt, en0, en1, amplifier,time, out_flag_clone)
    use com_main_gw
    implicit none
    type(chain_pointer_type)::pt
    integer amplifier
    real(8) en0, en1, time
    integer lvl, out_flag_clone,nl
    logical Istransit, Invtransit
    real(8) rnd, tmp
    out_flag_clone=0
    if(Istransit(en0, en1, lvl,nl))then
        if(ctl%chattery.ge.4)then
            !if(ctl%ntasks>1)then
            !    write(*,fmt=*) "time,en0,en1=",time/1d6/(2*pi), en0/ctl%energy0,&
            !        en1/ctl%energy0
            !    write(*,fmt=*) "crossing ",lvl,", create clone particle",  pt%idx
            !else
                print*, "time,en0,en1=",time/1d6/(2*pi), en0/ctl%energy0,en1/ctl%energy0
                print*, "crossing ",lvl,", create clone particle",  pt%idx
            !end if
        end if
        !print*, "ac_b:pt%ed%idx=",pt%ed%idx
        !print*, "transit, lvl=",lvl, en0, en1, en1/ctl%energy0, en0/ctl%energy0,log10clone_emax
        !read(*,*)
        call create_clone_particle(pt,lvl,amplifier**nl,time)
        !print*, "ac_e:pt%ed%idx=",pt%ed%idx
        !read(*,*)
        !Iscloned(lvl)=.true.
       ! print*, "clone"
        ctl%num_clone_created=ctl%num_clone_created+amplifier**nl-1
        out_flag_clone=1
        if(ctl%chattery.ge.4) write(*,*) "-------clone particle--",lvl,pt%idx, pt%ed%idx
    end if	
    !if(ctl%chattery.ge.3) print*, "???????", ctl%del_cross_clone
    if(ctl%del_cross_clone.ge.1)then
        !if(ctl%chattery.ge.3)then
        !    print*," nhiar, lvl=", sample%source,lvl
        !end if
        if(Invtransit(en0, en1,lvl,nl))then
            !print*, "inverse_transit, lvl=",lvl, en0, en1, en1/ctl%energy0, en0/ctl%energy0,log10clone_emax
            !read(*,*)
            tmp=rnd(0d0,dble(amplifier**nl))
            if(tmp>1d0)then
                !print*, "tmp=", tmp
                !read(*,*)
                if(ctl%chattery.ge.4)then
                    if(ctl%ntasks>1)then
                        write(unit=chattery_out_unit,fmt=*) "en0,en1=",en0/ctl%energy0,en1/ctl%energy0
                        write(unit=chattery_out_unit,fmt=*) &
                            "Invcrossing, deleting clone particle",  pt%idx
                    else
                        print*, "en0,en1=",en0/ctl%energy0,en1/ctl%energy0
                        print*, "Invcrossing, deleting clone particle",  pt%idx
                    end if
                end if
            !	read(*,*)
                ctl%num_clone_elim=ctl%num_clone_elim+1
                out_flag_clone=100
                return
            end if
            pt%ob%weight_real=pt%ob%weight_real*amplifier
            pt%ob%weight_clone=pt%ob%weight_clone*amplifier
            pt%ob%Lvl_clone=pt%ob%Lvl_clone-1
        end if		
    end if
end subroutine
