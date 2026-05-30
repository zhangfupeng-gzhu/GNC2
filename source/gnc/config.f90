module model_config
    type type_sub_para
        character*(50) name
        character*(200) str
    end type
    type type_para
        character*(50) name
        character*(200) str
        integer nsub
        type(type_sub_para),allocatable::sub_para(:)
        contains 
        procedure::init=>init_type_para
        procedure::get_sub_para
    end type
    

    type type_paras_all
        ! type(type_paras_real)::para_real
        ! type(type_paras_real_arr)::para_real_arr
        ! type(type_paras_int)::para_int
        ! type(type_paras_int_arr)::para_int_arr
        ! type(type_paras_log)::para_log
        type(type_para),allocatable::tp(:)
        integer n
    contains
        procedure::init=>init_type_paras_all
        procedure::print=>print_type_paras_all
        procedure::save_to_txt=>save_to_txt_type_paras_all
    end type
    
    type(type_paras_all)::pa_default, pa_usr_set, pa_now_used
    private::init_type_paras_all,get_sub_para
    private::print_type_paras_all,init_type_para,save_to_txt_type_paras_all
contains

subroutine get_sub_para(this,sub,str,ier)
    implicit none
    class(type_para)::this
    character*(50) str
    character*(*) sub
    integer i,ier
    ier=0
    do i=1, this%nsub
        if(trim(adjustl(sub))==trim(adjustl(this%sub_para(i)%name)))then
            str=trim(adjustl(this%sub_para(i)%str))
            return
        end if
    end do
    ier=-1
    print*, "error! no sub_para find:", trim(adjustl(sub))
    stop
end subroutine
subroutine init_type_paras_all(this,n)
    implicit none
    class(type_paras_all)::this
    integer n
    if(allocated(this%tp))deallocate(this%tp)
    allocate(this%tp(n))
    this%n=n
end subroutine
subroutine init_type_para(this,nsub)
    implicit none
    class(type_para)::this
    integer nsub
    if(allocated(this%sub_para))deallocate(this%sub_para)
    allocate(this%sub_para(nsub))
    this%nsub=nsub
end subroutine
subroutine save_to_txt_type_paras_all(this, fname)
    implicit none
    class(type_paras_all)::this
    character*(*) fname
    integer i,j
    integer,parameter::funit=1231123
    open(unit=funit,file=trim(adjustl(fname)))
    do i=1, this%n
        write(unit=funit,fmt="(A50,A80)") trim(adjustl(this%tp(i)%name))//"=",trim(adjustl(this%tp(i)%str))
        if(this%tp(i)%nsub.gt.0)then
            do j=1, this%tp(i)%nsub
                write(unit=funit,fmt="(A50,A80)") "sub:"//&
                trim(adjustl(this%tp(i)%sub_para(j)%name))//"=",trim(adjustl(this%tp(i)%sub_para(j)%str))
            end do
        end if
    end do
    close(unit=1231123)
end subroutine
subroutine print_type_paras_all(this, name)
    implicit none
    class(type_paras_all)::this
    character*(*) name
    integer i,j
    print*, "name=", trim(adjustl(name))
    do i=1, this%n
        !write(unit=*,fmt="(A50,A80)") "'"//trim(adjustl(this%tp(i)%name))//"'","'"//trim(adjustl(this%tp(i)%str))//"'"
        write(unit=*,fmt="(A50,A80)") trim(adjustl(this%tp(i)%name))//"=",trim(adjustl(this%tp(i)%str))
        if(this%tp(i)%nsub.gt.0)then
            !print*, "i,nsub=",i,this%tp(i)%nsub
            do j=1, this%tp(i)%nsub
                !write(unit=*,fmt="(A50,A80)") "sub:"//&
                !trim(adjustl(this%tp(i)%sub_para(j)%name))//"'","'"//trim(adjustl(this%tp(i)%sub_para(j)%str))//"'" 
                write(unit=*,fmt="(A50,A80)") "sub:"//&
                trim(adjustl(this%tp(i)%sub_para(j)%name))//"=",trim(adjustl(this%tp(i)%sub_para(j)%str))
            end do
        end if
    end do
    print*, "-----------------------------"
end subroutine
subroutine read_one_para(funit,pname,tp)
    implicit none
    type(type_para)::tp
    character*(200)  str_(2)
    integer::funit, ier
    !type(type_para)::tpa
    !open(unit=funit,file=trim(adjustl(fl)),status="old")
    integer nptot,i,j
!    logical flag_bg, flag_ed
    integer n_sub(1000)
    character*(2) sp
    character*(*) pname
    character*(200) str_sub(2,1000,10)


    call READPAR_STR_SPLIT_TWO(str_,funit,"#",ier)
    if(trim(adjustl(str_(1))).eq.trim(adjustl(pname)))then
        tp%name=trim(adjustl(str_(1)))
        tp%str=trim(adjustl(str_(2)))
    else
        print*, "error! expect parameter name:", trim(adjustl(pname))
        print*, "but we get:", trim(adjustl(str_(1))), "//"
    end if
            
end subroutine
subroutine read_paras(funit,pa)
    implicit none
    type(type_paras_all)::pa
    character*(200) str(2,1000),str_(2)!,str_trim(2)
    integer::funit, ier
    !type(type_para)::tpa
    !open(unit=funit,file=trim(adjustl(fl)),status="old")
    integer nptot,i,j
!    logical flag_bg, flag_ed
    integer n_sub(1000)
    character*(2) sp
    character*(200) str_sub(2,1000,10)
    nptot=0
    ier=0
!    flag_bg=.false.
!    flag_ed=.false.
    do while(ier.eq.0)        
        call READPAR_STR_SPLIT_TWO(str_,funit,"#",ier)
        if(trim(adjustl(str_(1))).ne.''.and.ier.eq.0)then
            sp=trim(adjustl(str_(1)))
            
            if(sp(1:2).eq.'--')then
                if(nptot.eq.0)then
                    print*, "error! sub para can not be the first one"
                    stop
                end if
                n_sub(nptot)=n_sub(nptot)+1
                !print*, "nptot=",nptot
                str_sub(1:2,nptot,n_sub(nptot))=str_(1:2)
                cycle
            end if
            nptot=nptot+1
            
            str(1:2,nptot)=str_(1:2)
            n_sub(nptot)=0
        end if
        !print*,"ier=",ier
    end do
    if(nptot>1000)then 
        print*, "error! too many parameters"
        stop
    end if
    !print*, "nptot=",nptot
    call pa%init(nptot)
    do i=1, nptot
        call remove_tab_in_string(trim(adjustl(str(1,i))),pa%tp(i)%name)
        call remove_tab_in_string(trim(adjustl(str(2,i))),pa%tp(i)%str)
        if(n_sub(i).gt.0)then
            call pa%tp(i)%init(n_sub(i))
           ! print*, "pa%tp(i)%nsub=",pa%tp(i)%nsub
            do j=1, n_sub(i)
                pa%tp(i)%sub_para(j)%name=trim(adjustl(str_sub(1,i,j)))
                pa%tp(i)%sub_para(j)%str=trim(adjustl(str_sub(2,i,j)))
            end do
        else
            pa%tp(i)%nsub=0
        end if
    end do
end subroutine
subroutine override_paras(pa_def, pa_set, pa_now)
    implicit none
    type(type_paras_all)::pa_def, pa_set, pa_now
    integer i,j,k,l

    pa_now=pa_def
loopi:    do i=1, pa_set%n
loopj:        do j=1, pa_now%n
                if(trim(adjustl(pa_set%tp(i)%name))==trim(adjustl(pa_now%tp(j)%name)))then
                    pa_now%tp(j)%str=pa_set%tp(i)%str
                    !print*, pa_now%tp(j)%nsub, pa_set%tp(i)%nsub
                    if(pa_now%tp(j)%nsub.ne.pa_set%tp(i)%nsub)then
                        call pa_now%tp(j)%init(pa_set%tp(i)%nsub)
                        !print*, "pa_now%n=",pa_now%tp(j)%nsub, size(pa_now%tp(j)%sub_para)
                    end if
                    pa_now%tp(j)%sub_para(1:pa_now%tp(j)%nsub)= &
                        pa_set%tp(i)%sub_para(1:pa_set%tp(i)%nsub)

                    cycle loopi
                end if
            end do loopj
        print*, "error! para not defined:", trim(adjustl(pa_set%tp(i)%name))
        stop
    end do loopi
end subroutine

end module