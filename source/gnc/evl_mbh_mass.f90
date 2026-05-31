

subroutine get_mbh_increase(dmbh_tot)
    use com_main_gw
    use md_mbh_evl_acc
    implicit none
    real(8) dmbh_tot, mass_acc,dt
	if(rid.eq.0)then
    	print*, "update mbh masses..."
	end if
    call get_mass_mbh_growth(bksams,ctl%run_snap_time_inner_steps_i)
    if(rid.eq.0)then
        print*, "mbh mass increasement"
        write(unit=*,fmt="(20A12)") "TD", "LC(star)","LC(SBH)", "EMAX(SBH)", "EMRI(SBH)", "Std(MS)", "Std(RG)"
        write(unit=*,fmt="(1P20E12.3)") mbh_mmg%td_disc(1), mbh_mmg%lc_direct(1:2),mbh_mmg%emax_direct(2), &
        mbh_mmg%emri(2), mbh_mmg%stc_disc(1),mbh_mmg%stc_disc(6)
    end if
    if(rid.eq.0)then
        print*, "gas_resurvior_now=",  mbh_mmg%gas_reservior_left
        
    end if
    dt=ctl%run_snap_time_inner_steps_f-ctl%run_snap_time_inner_steps_i
    
    if(mbh_mmg%gas_reservior_left.gt.0)then
        call get_mass_acc(spp_new%mbh,dt, mass_acc,ctl%radiative_efficiency)

        if(mbh_mmg%gas_reservior_left>mass_acc)then
            mbh_mmg%gas_reservior_left=mbh_mmg%gas_reservior_left-mass_acc
        else
            mass_acc=mbh_mmg%gas_reservior_left
            mbh_mmg%gas_reservior_left=0            
        end if
        if(rid.eq.0)then
            print*, "af:dt,mass acc,gas_resurvior=", dt, mass_acc, mbh_mmg%gas_reservior_left
            print*, "direct_swallow=",mbh_mmg%mass_direct_swallow
        end if
    else
        mass_acc=0
    end if
    mbh_mmg%gas_reservior_consum=mass_acc
    dmbh_tot=mbh_mmg%mass_direct_swallow+mass_acc*(1-ctl%radiative_efficiency)
    if(rid.eq.0)then            
        print*, "mbh(or->new)=", spp_new%mbh, spp_new%mbh+dmbh_tot, spp_new%mbh_dmless, (spp_new%mbh+dmbh_tot)/m0_cl
    end if
    mbh_mmg%mass_tot=dmbh_tot
    !mbh=mbh+sum(ctl%evl_mbh_mass(1:ctl%enable_evl_mbh_mass_comp))
    !mbh_dmless=mbh/m0_cl
    !call update_arrays_single(.true.)
    !call get_sample_para(dms,bksams_arr_norm,.false.,spp_new)

end subroutine

subroutine get_mass_mbh_growth(bks,cur_time_i)
    use com_main_gw
    use md_mbh_evl_acc
     
    implicit none
    integer n, i,idx
    real(8) cur_time_i
	type(chain_type)::bks
    type(chain_pointer_type),pointer::pt

    
    pt=>bks%head
    !if(rid.eq.0)then
    !    print*, "enable mbh evl=", ctl%enable_evl_mbh_mass(1:3)
    !end if
    do while(associated(pt))
        select type (ca=>pt%ob)
        type is (particle_sample_type)
            call mbh_mmg%add_event_stype(ca,pt%ob%weight_real*pt%ob%m)
        end select
        pt=>pt%next
    end do
    !print*, "rid, evl_mass_star, bh, stellar_other=",rid, evl_mass_star, evl_mass_bh, evl_mass_stellar_other
    
    call mpi_collect_evl_masses(mbh_mmg)
    call update_acum_data(mbh_mmg)
    call update_tot_data(mbh_mmg)
    mbh_mmg%gas_reservior_left=mbh_mmg%gas_reservior_left+mbh_mmg%gas_reservior_add

end subroutine
subroutine get_mass_acc(mbh,dt,mass_acc,epsilon)
    implicit none
    real(8) mbh, dt, mass_acc
    real(8) medd
    real(8),external::mdot_edd_msun_yr
    real(8) epsilon
    medd=mdot_edd_msun_yr(epsilon,mbh)
    mass_acc=medd*dt*1e6
end subroutine
subroutine mpi_collect_evl_masses(mmg)
    use md_mbh_evl_acc
    use MPI_comu
    implicit none
    real(8) fx(mmg_real_length)
    type(mass_mbh_growth)::mmg
    integer ierr
    call mmg%pack(fx)
    call collection_and_avg_fx(fx,mmg_real_length)
   call mpi_barrier(MPI_comm_world, ierr)
   call mmg%unpack(fx)
end subroutine
 
 