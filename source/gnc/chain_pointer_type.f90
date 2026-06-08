module md_chain_pointer
	use md_particle_sample
	implicit none
	type chain_pointer_type
		integer idx
!		class(*), allocatable, target:: ob
		class(particle_sample_type), allocatable:: ob
        type(chain_pointer_type),pointer:: Next =>null()
		type(chain_pointer_type),pointer:: Prev=>null()
		type(chain_pointer_type),pointer:: Ed=>null()
		type(chain_pointer_type),pointer:: Bg=>null() 
!	    integer flag_chain
    contains
        procedure::set_head=>set_item_chain_head_chain_type
        procedure::init_head=>init_chain_head_chain_type
        procedure::set_end=>set_item_chain_end_chain_type
        procedure::chain_to_arr_single=>chain_to_arr_chain_type_single
        procedure::arr_to_chain=>arr_to_chain_chain_type
		procedure::create_ob=>chain_create_ob_chain_type
		procedure::copy_to=>copy_chain_object_chain_type
		procedure::print_node=>print_node_chain_chain_type 
        procedure::attach_chain=>attach_chain_chain_type
        procedure::create_chain=>create_chain_chain_type
        procedure::get_sizeof=>get_sizeof_chain_chain_type
		!procedure::
	end type
!	integer,parameter::chain_pointer_bgn=1
!	integer,parameter::chain_pointer_end=2
!	integer,parameter::chain_pointer_cnt=0
    private::init_chain_head_chain_type
	private::chain_create_ob_chain_type,copy_chain_object_chain_type,&
         print_node_chain_chain_type, set_item_chain_end_chain_type
    private::chain_to_arr_chain_type_single, arr_to_chain_chain_type
    private:: attach_chain_chain_type
    private::get_sizeof_chain_chain_type,create_chain_chain_type
contains
	subroutine print_node_chain_chain_type(chain, str_)
		implicit none
		class(chain_pointer_type)::chain
		character*(*) str_
		print*, str_ 
        select type (ca=>chain%ob)
            class is (particle_sample_type)
                call ca%print(str_)
        end select 
		 
	end subroutine
    subroutine create_chain_chain_type(item,length) 
        class(chain_pointer_type),target::item
        type(chain_pointer_type),pointer::p=>null(), pc, pn
        integer length,i
        !print*, "create_chain_single:", chain%idx, chain%ob%id
        pc=>item
        if(.not.associated(pc))then
            print*, "pc not associated"
            stop
        end if	        
        if(.not.associated(item%next))then
            !item is the end of the chain
            p=>item
            do i=1, length
                allocate(pc%next)
                p=>pc%next
                p%prev=>pc
                p%bg=>pc%bg
				 
                p%idx=pc%idx+1
                pc=>pc%next
            end do    
            call set_item_chain_end_chain_type(p)
        else
            do i=1, length
                pn=>item%next
                allocate(p)
                pn%prev=>p
                item%next=>p
                p%prev=>item
                p%next=>pn
				 
                p%ed=>item%ed
                p%bg=>item%bg
                nullify(p)
            end do
        end if
       
        !print*, "p%idx=",p%idx, chain%ed%idx
        !call chain%output()
        !call chain%next%next%output()
        !read(*,*)
    end subroutine

	subroutine chain_to_arr_chain_type_single(chain, arr, n)
		implicit none
		class(chain_pointer_type),target::chain
		type(chain_pointer_type),pointer::p
		integer n,i
		type(particle_sample_type),dimension(n) ::arr

		p=>chain
		i=1
        do while(associated(p))
            if(.not.allocated(p%ob))then
                print*, "chain_to_arr:not allocated???i=", i, p%idx
                read(*,*)
            end if
            select type(ca=>p%ob)
                type is(particle_sample_type)
                    
                    arr(i)=ca
                    if(ca%write_down_track.gt.0) then
                        !print*, ca%length,ca%n_gene,ca%id
                        !print*, arr(i)%length
                    end if
                    arr(i)%idx=i
                    i=i+1
            end select
            p=>p%next
        end do
    endsubroutine
    
    subroutine arr_to_chain_chain_type(chain,  arr , n)
		implicit none
		class(chain_pointer_type),target::chain   
		type(chain_pointer_type),pointer::p
		integer n,i
		class(particle_sample_type) ::arr(n)
        !class(particle_sample_type),pointer::pa
		p=>chain
		i=1
		do while(associated(p))
            allocate(p%ob,source=arr(i))
            select type(pa=>arr(i))
            type is (particle_sample_type)
                select type(ca=>p%ob)
                type is(particle_sample_type)
                    ca=pa
                    ca%idx=i
                end select
            end select
			p=>p%next
			i=i+1
		end do
	end subroutine
	integer function get_sizeof_chain_chain_type(chain)
		implicit none
		class(chain_pointer_type),target::chain
		type(chain_pointer_type),pointer::p
		integer nsize
		p=>chain
		nsize=0
		do while(associated(p))
            select type(ca=>p%ob)
            class is (particle_sample_type)
			    nsize=nsize+sizeof(ca)
            end select
			p=>p%next
		end do
        get_sizeof_chain_chain_type=nsize
	end function
	
	recursive subroutine destroy_attach_pointer_chain_type(item)
		implicit none
		type(chain_pointer_type),target:: item
		type(chain_pointer_type),pointer:: ps
		!type(chain_particle_pointer_type),pointer:: left,right
		ps=>item
		!print*, "associated(left)?",associated(item%Append_left)
		 
		if(allocated(item%ob))then
			deallocate(item%ob)
		end if
		!print*,"deallocate",associated(item)
		deallocate(ps)
		!print*, "deallocate(item) sucess"
        !nullify(item)
	end subroutine

	recursive subroutine copy_chain_object_chain_type(item,cp)
		implicit none
		class(chain_pointer_type)::item
		type(chain_pointer_type):: cp
		
		 
        if(.not.allocated(item%ob)) then
            print*, "error! item%ob not allocated"
            stop
        end if
        select type(ca=>item%ob)
        type is(particle_sample_type)
            if(.not.allocated(cp%ob)) &
            allocate(particle_sample_type::cp%ob)
            select type(cb=>cp%ob)
            type is(particle_sample_type)
                cb=ca
            end select
        
        end select
    
		
		
	end subroutine

	recursive subroutine save_attach_points(sp,fileunit,version_number)
	implicit none
	integer fileunit,version_number
	type(chain_pointer_type)::sp
	integer,parameter::flag_left=0,flag_right=1, flag_write=2,flag_none=-1
    integer,parameter::flag_particle=3

	!if(sp%ob%N_gene.ge.2)then
	!	print*, "save_attach_points:sp%ob%N_gene", sp%ob%N_gene
		!call sp%print_node("save_attach_points")
	!end if 
    select type(ca=>sp%ob)
    type is(particle_sample_type) 
	    call ca%write_info(fileunit,version_number)
    end select
	!print*, "write, sizeof=", flag_write, sizeof(sp%ob%track), allocated(sp%ob%track)
end subroutine

recursive subroutine read_attach_points(sp,fileunit,version_number)
	implicit none
	integer fileunit,version_number
	type(chain_pointer_type)::sp
	!type(chain_particle_pointer_type),pointer:: left,right
	integer flag
	integer,parameter::flag_left=0,flag_right=1, flag_write=2,flag_none=-1
    integer,parameter::flag_particle=3 
  
	 
	if(.not.allocated(sp%ob)) then
		allocate(particle_sample_type::sp%ob)
	end if
	select type (ca=>sp%ob)
	type is (particle_sample_type)
		call ca%read_info(fileunit,version_number) 
	end select 
	
end subroutine

	subroutine attach_pointer_in_a_chain(item, pt, flag)
		implicit none
		type(chain_pointer_type),pointer::item, pt, prev, next
		integer flag
		if(.not.associated(pt))then
			print*, "attach_pointer:error, pt not allocated"
			stop
		end if
		if(associated(pt%prev))then 
			prev=>pt%prev
		else
			!The pt is the beginning of a chain
			write(*,*) "del:pt is the beginning of a chain"
			if(.not.associated(pt%next))then
				print*, "warnning, chain have only one item"
			end if
			call set_item_chain_head_chain_type(pt%next)
			pt%next%prev=>null()
			!print*, next%bg%idx,associated(item)
			goto 100
		end if

		if(associated(pt%next))then
			next=>pt%next
		else
			!The item is the end of a chain
			!	write(*,*) "del:item is the end of a chain"
			call set_item_chain_end_chain_type(pt%prev)
			prev%next=>null()
			goto 100
		end if

		prev%next=>next
		next%prev=>prev

100		pt%next=>Null()
		pt%prev=>Null()	
		pt%bg=>Null()
		pt%ed=>Null()
		 
	end subroutine
	 
	subroutine chain_pointer_delete_item_chain_type(item)
		implicit none
		type(chain_pointer_type),target:: item
		type(chain_pointer_type),pointer:: prev,next,ps
		!print*, "delete item single"
		ps=>item
		if(associated(ps))then
			!print*, associated(item%prev)
			if(associated(ps%prev))then
				prev=>ps%prev
			else
				!The item is the beginning of a chain
				write(*,*) "del:item is the beginning of a chain"
                stop
				if(.not.associated(ps%next))then
					print*, "warnning, chain have only one item"
				end if
				!call chain_output_list(item%bg)
				!call chain_pointer_delete_item(next)
                next=>ps%next
                call next%set_head()
				next%prev=>null()
				!print*, next%bg%idx,associated(item)
				!===================================
				call destroy_attach_pointer_chain_type(ps)
				!deallocate(item)
				!===================================
				return
			end if
			!print*, associated(item%next)
			if(associated(ps%next)) then
				next=>ps%next
			else
				!The item is the end of a chain
				!write(*,*) "del:item is the end of a chain"
				call prev%set_end()
				prev%next=>null()
				!===================================
				call destroy_attach_pointer_chain_type(ps)
				!deallocate(item)
				!===================================
				return
			endif

			prev%next=>next			
			next%prev=>prev
			!print*, "destroy_attach_pointer"
			!print*, "associated item?", associated(ps)
			!===================================
			!deallocate(item)
			call destroy_attach_pointer_chain_type(ps)
			!===================================
		else
			write(*,*) "error, item not allocated"
			stop
		end if
		
	end subroutine
	
	subroutine chain_create_ob_chain_type(sp)
		implicit none
		class(chain_pointer_type),target::sp
		type(chain_pointer_type),pointer::p
		p=>sp
		do while(associated(p))
			if(allocated(p%ob)) deallocate(p%ob)
			allocate(particle_sample_type::p%ob)
			p=>p%next
		end do
	end subroutine
	subroutine set_item_chain_end_chain_type(item)
		class(chain_pointer_type),target::item
		type(chain_pointer_type),pointer::p
		p=>item
		p%ed=>item
		do while(associated(p%prev))
			p=>p%prev
			p%ed=>item
		end do
	end subroutine
	subroutine set_item_chain_head_chain_type(item)
		class(chain_pointer_type),target::item
		type(chain_pointer_type),pointer::p
		!print*, "0"
		!p=>item
		!print*, associated(item)
        p=>item
        p%bg=>item
		!print*, "1"
		do while(associated(p%next))
			!print*, "2"
			p=>p%next
			!print*, "3"
			p%bg=>item
		end do
	end subroutine
	subroutine init_chain_head_chain_type(item)
		class(chain_pointer_type),target:: item
!			item%prev=>Null
!			item%next=>Null
			item%Bg=>item
			item%Ed=>item
			item%idx=0
!			read(*,*) 'error: item not created'
	end subroutine
	subroutine insert_after_item_chain_pointer(item) 
		type(chain_pointer_type),pointer:: item
		type(chain_pointer_type),pointer:: p, pn

		if(.not.associated(item))then
			! item is the begining of a chain
			!print*, "item is the beginning of a chain"
			allocate(item)
			call init_chain_head_chain_type(item)
			return
		end if

		if(.not.associated(item%next))then
			! item is the end of a chain
			!print*, "item is the end of a chain"
			allocate(item%next)
			p=>item%next
			p%prev=>item
			p%bg=>item%bg
			p%idx=item%idx+1
			call set_item_chain_end_chain_type(p)
			!print*, associated(p%ed)
			!call chain_output_list(item%bg)
			return
		end if
       ! print*, "1"
       ! print*, "item%idx=", item%idx
		pn=>item%next
		allocate(p)
		pn%prev=>p
		item%next=>p

		p%prev=>item
		p%next=>pn

		p%ed=>item%ed
		p%bg=>item%bg

		p%idx=item%idx+1
       ! print*, "p%idx=",p%idx
		p=>p%next	    
		do while(associated(p))
			p%idx=p%idx+1
            p=>p%next
		end do
       ! print*, "insert finished"
	end subroutine

	subroutine attach_chain_chain_type(this, c)
        class(chain_pointer_type):: this
		type(chain_pointer_type):: c
        !attach c to the head of this chain
        type(chain_pointer_type),pointer::phead, pend,p, chead, cend
		!print*, "1"
		!call this%output()
		!print*, "2"
        phead=>this%bg
        pend=>this%ed
       ! print*,"3"
        chead=>c%bg
		cend=>c%ed
       ! print*, "this%bg%bg%idx=",this%bg%bg%idx
        pend%next=>chead
        chead%prev=>pend
        !phead%idx=pend%idx+1
		
		!stop
		!print*,"4"
		!print*,"2:", associated(phead)
		!print*, "3:", associated(this%bg)
        p=>chead
        do while(associated(p))
            p%idx=p%prev%idx+1
            p=>p%next
        end do
		!print*,"5"
		!call this%output()
		!print*,"??"
		!print*, associated(phead)
        call phead%set_head()
        call cend%set_end()
        
		!print*, "7"
    end subroutine 

end module
