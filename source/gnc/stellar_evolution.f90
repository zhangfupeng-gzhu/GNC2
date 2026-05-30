module md_stellar_evolution
    use md_particle_sample
    implicit none

contains        
    subroutine get_step_stellar_evl(sp,steps,ctime,debug)
        implicit none
        class(particle_sample_type)::sp
        real(8) ctime,steps,delta_next_evl
        logical debug
        character*(5) star_type_str
        real(8),parameter::tiny=1d-6

        if(sp%sh%cur_idx.eq.0)then
            print*, "error! sp%sh%cur_idx=0"
            stop
        end if
        if(sp%sh%cur_idx.eq.sp%sh%n)then
            return
        endif
        delta_next_evl=sp%sh%time(sp%sh%cur_idx+1)*2*pi*1d6-ctime
        if(delta_next_evl<=0)then
            return
        end if
        !delta_next_evl=max(delta_next_evl,tiny)
        steps=min(steps,delta_next_evl/sp%period*(1+tiny))
        !!!
        if(delta_next_evl.eq.0)then
            print*, "error! steps=0 in stellar evl"
            call sp%sh%print()
            print*, "ctime,time=",ctime/2d0/pi/1d6,sp%sh%time(sp%sh%cur_idx+1)
            print*, "type=", star_type_str
            print*, "cur_idx=",sp%sh%cur_idx
            print*, "delta_next_evl, evl_step, step=",delta_next_evl/2d6/pi, delta_next_evl/sp%period+tiny, steps
            stop
        end if
        if(debug)then
            print*, "===========stellar evolution step"
            !call sp%sh%print()
            call get_star_type(sp%obtype,star_type_str)
            print*, "ctime,time=",ctime/2d0/pi/1d6,sp%sh%time(sp%sh%cur_idx+1)
            print*, "type=", star_type_str
            print*, "cur_idx=",sp%sh%cur_idx
            print*, "delta_next_evl, evl_step, step=",delta_next_evl/2d6/pi, delta_next_evl/sp%period, steps
            print*, "============================"
            read(*,*)
        end if
    end subroutine
    subroutine prepare_bse_brown_data()
        implicit none
        integer i
        brown_data_age=(/0d0,1000d0,5000d0,10000d0/)
        do i=1, 4
            call prepare_table_of_brown_dwarf_radius(brown_data_age(i),brown_data_xs(1:14,i),&
            brown_data_ys(1:14,i),brown_data_y2(1:14,i))
        end do
    end subroutine
    subroutine get_brown_dwarf_radius_at_age(age, mass, radius)
        implicit none
        real(8) age, mass, radius
        integer i,idx
        if(age<=brown_data_age(1)) then
            idx=1
        elseif(age>=brown_data_age(4)) then
            idx=4
        else
            do i=1, 3
                if(age>brown_data_age(i).and.age<brown_data_age(i+1))then
                    idx=i
                end if
            end do
        end if
        call get_brown_dwarf_radius(mass,brown_data_xs(1:14,idx),brown_data_ys(1:14,idx),brown_data_y2(1:14,idx),radius)
    end subroutine
    subroutine add_mass_loss_to_gas_reservior(sp,mass_change)
        use md_mbh_evl_acc
        use com_main_gw
        implicit none 
        real(8) mass_change
        type(particle_sample_type)::sp
        mbh_mmg%mass_loss_stellar_evolution=mbh_mmg%mass_loss_stellar_evolution&
            +mass_change*ctl%stellar_evolution_fraction_mass_loss_to_reservior*sp%weight_real
        
    end subroutine
    
    subroutine get_current_particle_stellar_info(sp,ctime, flag, debug)
        use md_mobse_stellar_single,only:get_kstar_type!,get_kw_type
        use md_mbh_evl_acc,only:mbh_mmg
        implicit none
        class(particle_sample_type)::sp
        real(8) ctime
        integer ktype,kwtype, flag
        real(8) mass, radius
        real(8) pre_mass
        logical debug
        character*(10) ktype_str, star_type!, kwtype_str
        flag=0
        if(debug)then
            print*, "=========================update stellar info"
            call sp%sh%print()
        end if
        associate(sh=>sp%sh)
            if(sh%cur_idx==0)then
                call get_current_stellar_info(sh,ctime,ktype,kwtype,mass,radius)
            else
                if(sh%cur_idx<sh%n)then
                    if(ctime>=sh%time(sh%cur_idx+1))then
                        if(debug)then
                            print*, "ctime>=time_next, update", ctime, sh%time(sh%cur_idx+1)
                        end if
                        pre_mass=sh%mass(sh%cur_idx)
                        call get_current_stellar_info(sh,ctime,ktype,kwtype,mass, radius)
                        flag=1
                        if(pre_mass-mass<-1d-3)then
                            print*, "error! mass change<0,mass, pre_mass=",mass, pre_mass
                            call sh%print()
                            stop
                        endif
                        call add_mass_loss_to_gas_reservior(sp,pre_mass-mass)
                        if(debug)then
                            print*, "add gas to reservior=", (pre_mass-mass),sp%weight_real
                            print*, "current reservior=",mbh_mmg%mass_loss_stellar_evolution
                        end if
                    else
                        ! do nothing
                        if(debug)then
                            print*, "ctime<time_next, do nothing", ctime, sh%time(sh%cur_idx+1)
                        end if
                        return
                    end if
                else
                    if(debug)then
                        print*, "ctime>time_next, do nothing", ctime, sh%time(sh%cur_idx+1)
                    end if
                    return
                end if
            end if
    
            select case(ktype)
            case(0,1)
                if(sp%m<0.1d0)then
                    sp%obtype=star_type_bd
                else
                    sp%obtype=star_type_ms
                end if
                
            case(7,8,9)
                sp%obtype=star_type_nakedHe
            case(2,3,4,5,6)
                sp%obtype=star_type_rg
            case(10,11,12)
                sp%obtype=star_type_wd
            case(13)
                sp%obtype=star_type_ns
            case(14)
                sp%obtype=star_type_bh
            case(15)
                flag=-1 ! destroy
                if(debug)then
                    ktype_str=get_kstar_type(ktype)
                    print*, "ktype=",ktype_str
                    print*, "destroy!" 
                end if
                return
            end select
            !call get_type_from_ctl_obidx_sg(sp%obidx,sp%obtype)
            call get_obidx_from_type_sg(sp%obtype,sp%obidx)
            sp%m=mass
            
            sp%byot%ms%m=sp%m
            sp%byot%ms%radius=radius
            sp%byot%ms%obtype=sp%obtype
            sp%byot%ms%obidx=sp%obidx			
            if(debug)then
                ktype_str=get_kstar_type(ktype)
                call get_star_type(sp%obtype,star_type)
                print*, "update! current type, mass,radius=",ktype_str, mass, radius
                print*, "star type become=", star_type
                print*, "current idx=", sh%cur_idx 
            end if
        end associate
        if(debug)then
            print*, "============================================"
        end if
    end subroutine
end module
