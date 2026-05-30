

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


subroutine get_delta_e_mbh_evl(de,fphi, frho, fma, fjc,spd, sra, srp)
	use com_main_gw
	implicit none
	type(s1d_type)::  fma,fjc,fphi, frho
	!type(s1d_type)::aux
	type(s2d_type)::spd, sra, srp, de
	real(8) enx, ra, rp, pd,  jm,yout,jc
	integer idid, i,j
    do i=1, de%nx
		do j=1, de%ny
            enx=10**de%xcenter(i)
			jm=10**de%ycenter(j)
			ra=sra%fxy(i,j)
			rp=srp%fxy(i,j)
			pd=spd%fxy(i,j)
			jc=fjc%fx(i)

            call get_aux_function_for_period_pi2(common_aux,fphi, frho, fma,&
             enx,jm,jc,rp,ra,spp_new)
            !call aux%prepare_spline()
            !print*, "start"
            yout=0
            call my_integral_none(0d0,pi/2d0,yout,fcn,idid)
			de%fxy(i,j)=yout*2/pd

		end do
	end do	
	!call de%print("de")
contains
    subroutine fcn(n, x, y, f, par, ipar)
        use, intrinsic :: ieee_arithmetic
        implicit none
        integer n, ipar(100)
        real(8) x, y(n), f(n), par(100)
        real(8) aux_tmp,r

        !call splint_mylib(aux%xb,aux%fx,aux%y2a,aux%nbin, sin(x)**2, aux_tmp)
        call common_aux%get_value_s(x, aux_tmp)
        r=rp+(ra-rp)*sin(x)**2        
        f(1)=aux_tmp/r
        if(ieee_is_nan(f(1)).or.r.eq.0)then
            call common_aux%print("aux")
            call fma%print("fma")
            call fjc%print("jc")
            call srp%print("rp")
            call sra%print("ra")
            print*, "ex,jm=",enx,jm
            print*, "r,aux_tmp, rp, ra=",r,aux_tmp, rp, ra
            stop
        end if
    end subroutine
end subroutine

subroutine get_delta_e_mbh_evl_indvd_replace(sps,dmbh,fphi, frho, fma)
	use com_main_gw
	implicit none
    type(particle_samples_arr_type)::sps
	type(s1d_type)::  fma,fphi, frho
	real(8) enx, ra, rp, pd,  jm,yout,jc, de, dmbh
    integer idid, i,j

    do i=1, sps%n
        if(sps%sp(i)%x>emax_factor.and.spp_new%mbh_dmless.eq.0)then
            cycle
        end if
        !print*, "sps%sp(i)%x,emax_factor",sps%sp(i)%x,emax_factor,mbh_dmless
        call get_delta_e_mbh_evl_indvd_one(sps%sp(i),fphi,frho,fma,de)
        sps%sp(i)%x=sps%sp(i)%x+de*dmbh
        sps%sp(i)%en=sps%sp(i)%x*ctl%energy0
    end do
end subroutine

subroutine get_delta_e_mbh_evl_indvd_one(sp,fphi, frho, fma,de)
	use com_main_gw
	implicit none
    type(particle_sample_type)::sp
	type(s1d_type)::  fma,fphi, frho
	real(8) enx, ra, rp, pd,  jm,yout,jc, de
    integer idid, i,j

    enx=sp%x
    jm=sp%jm
    jc=sp%jc/(r0_cl*ctl%v0)
    rp=sp%rp/r0_cl
    ra=sp%ra/r0_cl
    pd=sp%period/(r0_cl/ctl%v0)

    call get_aux_function_for_period_pi2(common_aux,fphi, frho, fma,&
        enx,jm,jc,rp,ra,spp_new)
    if(common_aux%fx(common_aux%nbin)<-1d90)then
        print*, "sp%id=",sp%id
        stop
    end if
    !call aux%prepare_spline()
    !print*, "start"
    yout=0
    call my_integral_none(0d0,pi/2d0,yout,fcn,idid)
    de=yout*2/pd
contains
    subroutine fcn(n, x, y, f, par, ipar)
        use, intrinsic :: ieee_arithmetic
        implicit none
        integer n, ipar(100)
        real(8) x, y(n), f(n), par(100)
        real(8) aux_tmp,r

        !call splint_mylib(aux%xb,aux%fx,aux%y2a,aux%nbin, sin(x)**2, aux_tmp)
        call common_aux%get_value_s(x, aux_tmp)
        r=rp+(ra-rp)*sin(x)**2        
        f(1)=aux_tmp/r
        if(ieee_is_nan(f(1)).or.r.eq.0)then
            call common_aux%print("aux")
            !call fma%print("fma")
           ! call fjc%print("jc")
            !call srp%print("rp")
            !call sra%print("ra")
            print*, "ex,jm=",enx,jm
            print*, "r,aux_tmp, rp, ra=",r,aux_tmp, rp, ra, sp%id
            stop
        end if
    end subroutine
end subroutine