module md_by_particle
    use md_chain
    use md_particle_sample    
	
    type samples_type_pointer
		type(chain_pointer_type),pointer::sp=>null()
		integer rid, index
	end type
    type samples_type_pointer_arr
        integer n
        type(samples_type_pointer),allocatable::pt(:)
        contains
            procedure::init=>init_pointer_arr
    end type
contains
    subroutine init_pointer_arr(this,n)
        implicit none
        class(samples_type_pointer_arr)::this
        integer n
        if(allocated(this%pt))deallocate(this%pt)        
        allocate(this%pt(n))
        this%n=n
    end subroutine

    subroutine convert_sams_pointer_arr(sps, sps_pointer,type)
		!use com_main_gw
		implicit none	
		type(chain_type)::sps
		type(samples_type_pointer_arr)::sps_pointer
		type(chain_pointer_type),pointer::ps,psout
		integer nsel,i, n
        integer,optional::type
		integer typeI

		nsel=0
        typeI=0
		if(present(type))then
            typeI=type
        end if
        call sps%get_length(nsel,type=typeI)
		call sps_pointer%init(nsel)
		nsel=0
        ps=>sps%head
        do while (associated(ps))
            select case(typeI)
            case(0)
                nsel=nsel+1
                sps_pointer%pt(nsel)%sp=>ps        
            case(1)
                select type(ca=>ps%ob)
                type is(particle_sample_type)
                    nsel=nsel+1
                    sps_pointer%pt(nsel)%sp=>ps        
                end select
             
            end select
            ps=>ps%next
        end do
	
	end subroutine
subroutine chain_select(sps, sps_out, exitflag, obj_type)
	implicit none	
	type(chain_type)::sps, sps_out
	type(chain_pointer_type),pointer::ps,psout
	integer nsel,i, exitflag
	integer,optional:: obj_type
    integer objtype
    objtype=0
    if(present(obj_type)) objtype=obj_type

    call chain_select_by_condition(sps,sps_out, selection)

contains
	logical function selection(pt)
		implicit none
		type(chain_pointer_type)::pt
        logical cond
        cond=.false.
        select case(objtype)
        case(0)
            cond=.true.
        case(1)
            select type(ca=>pt%ob)
            type is(particle_sample_type)
               cond=.true.
            end select
         
        end select
        if(((pt%ob%exit_flag.eq.exitflag).or.exitflag.eq.-1)&
            .and.cond)then
            selection=.True.
        else
            selection=.False.
        end if
	end function
end subroutine



subroutine chain_select_replace(sps, exitflag, obj_type)
	implicit none	
	type(chain_type)::sps, sps_out
	type(chain_pointer_type),pointer::ps,psout
	integer nsel,i, exitflag
	integer,optional:: obj_type
    integer objtype
    objtype=0
    if(present(obj_type)) objtype=obj_type

    !call chain_select_replace_by_condition(sps, selection)

contains
	logical function selection(pt)
		implicit none
		type(chain_pointer_type)::pt
        logical cond
        cond=.false.
        select case(objtype)
        case(0)
            cond=.true.
        case(1)
            select type(ca=>pt%ob)
            type is(particle_sample_type)
               cond=.true.
            end select
        
        end select
        if(((pt%ob%exit_flag.eq.exitflag).or.exitflag.eq.-1)&
            .and.cond)then
            selection=.True.
        else
            selection=.False.
        end if
	end function
end subroutine

end module