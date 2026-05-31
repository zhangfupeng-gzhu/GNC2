

subroutine get_system_dstr(bks,reset_weighting,cri_max)
    use com_main_gw
    implicit none
    integer n
    type(particle_samples_arr_type)::bks
    type(s1d_type)::fmden_new
	integer ierr,i,out_flag
    !integer,parameter::nmax=15
    integer niter
    !type(diffuse_mspec)::dms_tmp
    real(8) cri,cri_max
    logical::reset_weighting
    type(star_pot_para)::spp0
    real(8) t1, t2
    call mpi_barrier(mpi_comm_world,ierr)
    if(rid.eq.0)then
        call cpu_time(t1)
    end if 
    spp0=spp_new
    if(rid.eq.0.and.ctl%chattery.ge.2)then
        call dms%all%all%fmden%print("fmden old")
        !call dms%mb(1)%star%barj%print("barj old")
    end if
    call get_total_sample_number()
    call get_sample_erange(bks)
    !call set_edstr_bound()
     !print*, bks%sp(1:10)%x
     call get_ge(dms, bks)
    !call gen_nxj_fxj_gx_both(dms)
    niter=0
    
100 call get_nxgx_density(spp_new)
    
    fmden_new=dms%all%all%fmden     
    if(rid.eq.0)then
        print*, "sample_logrmin,logrmax=", sample_logrmin, sample_logrmax!,dms%logrmin,dms%logrmax
    end if 
    call get_rho_rmin(spp_new)    

    call init_diffuse_mspec_rtables(dms)
    call spp_new%init(dms%logrmin,dms%logrmax,ctl%dstr_bins_r)
    spp_new%frho_star=fmden_new
    
    spp_new%has_set_density=.true.
    call get_spp_from_rho(spp_new)   
    !fphi_new=spp_new%fphi_star

    
    call update_re_tables()
    call get_phi_diff(spp_new,spp0,cri)
    !call dms%mb(1)%all%barge%print("barge or")
    if(rid.eq.0.and.ctl%chattery.ge.2)then
        !call dms%fphi_star%print("phi_org")
        !if(ctl%chattery.ge.1)then
        call spp_new%frho_star%print("fmden new")
        call spp_new%fphi_star%print("phi new")
        !call dms%mb(1)%star%barj%print("barj new")
        !stop
    end if
    if(rid.eq.0)then
        print*, "get_self_consistent_pot:cri=",cri, niter
    end if
    
    if(cri>cri_max.and.niter<max_self_con_iter) then
        niter=niter+1 
        spp0=spp_new
        goto 100
    end if
    !call get_orbit_tables()
    call mpi_barrier(mpi_comm_world,ierr)
    if(rid.eq.0)then
        call cpu_time(t2)
        print*, "get_system_dstr used time:", t2-t1, " s"
    end if 
end subroutine
 
subroutine get_nxgx_density(spp)
    use com_main_gw
    implicit none
    integer ierr
    type(star_pot_para)::spp
 
     
	call gen_nxj_fxj_gx_iregular(dms,spp)  

    call mpi_barrier(mpi_comm_world,ierr) 
    call collect_dms_gx(dms)
    if(ctl%chattery.ge.3.and.rid.eq.0)then
        print*, "start get_cluster_density"
    end if
    call get_cluster_density(dms)
 
end subroutine
 
 
subroutine adb_cor_indvd_replace(bks_org, spp0,sppc,adb_method)
    use com_main_gw
    implicit none
    type(particle_samples_arr_type)::bks_org
    type(star_pot_para)::sppc,spp0
    real(8) de_tmp, t1,t2, jph_dmless, jc_xy, jmor
    integer i, ierr
    !integer,parameter::adb_method_acc=1,adb_method_fast=2,adb_method_nt=3
    !integer::adb_method=adb_method_fast
    integer adb_method

    if(rid.eq.0)then
        call cpu_time(t1)
    end if

    do i=1, bks_org%n
        if(spp0%mbh_dmless.eq.0)then
            if(bks_org%sp(i)%x>emax_factor)then
                print*, "bks_org%x,emax=",bks_org%sp(i)%x, emax_factor
                cycle
            end if
        end if
        !if(mod(i,1000).eq.0)print*, "i=",i,rid
        select case(adb_method)
        case(adb_est_method_acc)
            call get_dx_invariant_rad_action(bks_org%sp(i)%x,bks_org%sp(i)%jm,&
            bks_org%sp(i)%jc/(ctl%v0*r0_cl), spp0,sppc,dms%fr_phi,&
                bks_org%sp(i)%ra/r0_cl,bks_org%sp(i)%rp/r0_cl,bks_org%sp(i)%raq,de_tmp)
             
            
        case(adb_est_method_fast)
            call get_adb_cor_one(bks_org%sp(i),spp0,sppc,de_tmp)
 
        end select
       
        !    read(*,*)
        !read(*,*)
        bks_org%sp(i)%x=bks_org%sp(i)%x+de_tmp
        bks_org%sp(i)%en=bks_org%sp(i)%x*ctl%energy0
        jph_dmless=bks_org%sp(i)%jph/(ctl%v0*r0_cl)
        !bks_org%sp(i)%jm=bks_org%sp(i)%jph/bks_org%jc
        !print*, bks_org%sp(i)%x,bks_org%sp(i)%jm,bks_org%sp(i)%jc/(ctl%v0*r0_cl)
        if(bks_org%sp(i)%x>emin_factor.and.bks_org%sp(i)%x<emax_factor)then
            jmor=bks_org%sp(i)%jm
            
            call update_jm(dms,sppc,bks_org%sp(i)%x,jph_dmless,bks_org%sp(i)%jm,jc_xy)
            ! if(jmor>0.99)then
            !     print*, "bf:jm=", jmor, "af:jm=",bks_org%sp(i)%jm
            ! end if
        end if
        !print*, bks_org%sp(i)%x, bks_org%sp(i)%jm,jc_xy
        !read(*,*)
        !print*, "1:raq=",bks_org%sp(i)%raq
    end do
    !print*, "2"
    if(rid.eq.0)then
        call cpu_time(t2)
    end if
    call mpi_barrier(mpi_comm_world,ierr)
    if(rid.eq.0)then
        print*, "adb_cor_indvd, time=", t2-t1, " s"
    end if
end subroutine
 
subroutine get_adb_cor_one(sp,spp0,sppc,de)
    use com_main_gw
    type(particle_sample_type)::sp
    type(s1d_type)::fphi1, fphi0,frho0,fma0
    !type(s1d_type)::aux
    real(8) de, jc, rp, ra,pd, r, logr, phi_new,phi_old,yout,Dv, Jph_dm
    integer i,idid
    type(star_pot_para)::spp0,sppc


    !real(8) Ap, Aa, Bp, Ba
    !real(8) phic_rp, phic_ra, phi0_rp, phi0_ra
    !real(8) betac_ra,betac_rp, delta_phi_rp, delta_phi_ra
    !real(8) delta_m, Gv_rp, Gv_ra
    !real(8) delta_r_a, delta_r_p

    jc=sp%jc/(r0_cl*ctl%v0)
    rp=sp%rp/r0_cl
    ra=sp%ra/r0_cl
    pd=sp%period/(r0_cl/ctl%v0)

    if(rp.eq.ra)then
        call get_phi_star_full_range(sppc,log10(rp),phi_new)
        call get_phi_star_full_range(spp0,log10(rp),phi_old)
        de=10**phi_new-10**phi_old+sppc%mbh_dmless/rp-spp0%mbh_dmless/rp
        return
    end if
    call get_aux_function_for_period_pi2(common_aux,spp0,sp%x,sp%jm,jc,rp,ra)
    !call aux%prepare_spline()
    !de=0
    !do i=1, aux%nbin
    !    r=rp+(ra-rp)*sin(aux%xb(i))**2
    !    logr=log10(r)
    !    call get_phi_star_full_range(fphi_new,logr,phi_new)
    !    call get_phi_star_full_range(fphi_old,logr,phi_old)
    !    de=de+2d0/pd*aux%xstep*aux%fx(i)*(10**phi_new-10**phi_old)
    !end do
    yout=0
    if(spp0%mbh_dmless.eq.0.and.sppc%mbh_dmless.eq.0)then
        call my_integral_acc(0d0,pi/2d0,yout,pd_int_acc_a,pd_int_acc_r, fcn2,idid)
    else
        call my_integral_acc(0d0,pi/2d0,yout,pd_int_acc_a,pd_int_acc_r, fcn,idid)
    end if
    de=yout*2/pd
 
contains
	subroutine fcn(n, x, y, f, par, ipar)
		use, intrinsic :: ieee_arithmetic
		implicit none
		integer n, ipar(100)
		real(8) x, y(n), f(n), par(100)
		real(8) aux_tmp, phi_c, phi_0, r,logr

		!call splint_mylib(aux%xb,aux%fx,aux%y2a,aux%nbin, sin(x)**2, aux_tmp)
		call common_aux%get_value_s(x, aux_tmp)
		r=rp+(ra-rp)*sin(x)**2
        logr=log10(r)
		!call fphi_new%get_value_s(logr, phi_c)
        call get_phi_star_full_range(sppc,logr,phi_c)
        call get_phi_star_full_range(spp0,logr,phi_0)
		f(1)=aux_tmp*(10**phi_c-10**phi_0+sppc%mbh_dmless/r-spp0%mbh_dmless/r)
        !if(rid.eq.0)then
        !    print*, "r,phi_c,phi_0,phi_mbhc,phi_mbh0=",r,phi_c,phi_0,sppc%mbh_dmless/r,spp0%mbh_dmless/r
        !end if
	end subroutine   
    subroutine fcn2(n, x, y, f, par, ipar)
		use, intrinsic :: ieee_arithmetic
		implicit none
		integer n, ipar(100)
		real(8) x, y(n), f(n), par(100)
		real(8) aux_tmp, phi_c, phi_0, r,logr

		!call splint_mylib(aux%xb,aux%fx,aux%y2a,aux%nbin, sin(x)**2, aux_tmp)
		call common_aux%get_value_s(x, aux_tmp)
		r=rp+(ra-rp)*sin(x)**2
        logr=log10(r)
		!call fphi_new%get_value_s(logr, phi_c)
        call get_phi_star_full_range(sppc,logr,phi_c)
        call get_phi_star_full_range(spp0,logr,phi_0)
		f(1)=aux_tmp*(10**phi_c-10**phi_0)
        !if(rid.eq.0)then
        !    print*, "r,phi_c,phi_0,phi_mbhc,phi_mbh0=",r,phi_c,phi_0,sppc%mbh_dmless/r,spp0%mbh_dmless/r
        !end if
	end subroutine 
    !subroutine fcn3(n, x, y, f, par, ipar)
	!	use, intrinsic :: ieee_arithmetic
	!	implicit none
	!	integer n, ipar(100)
	!	real(8) x, y(n), f(n), par(100)
	!	real(8) aux_tmp, phi_c, phi_0, r,logr, beta_c, delta_r
!
	!	!call splint_mylib(aux%xb,aux%fx,aux%y2a,aux%nbin, sin(x)**2, aux_tmp)
	!	call common_aux%get_value_s(x, aux_tmp)
	!	r=rp+(ra-rp)*sin(x)**2
    !    logr=log10(r)
	!	delta_r=delta_r_p+Dv*(r-rp)
    !    call get_beta_full_range(sppc,logr, beta_c)
	!	f(1)=aux_tmp*(Jph_dm**2/r**3-beta_c/r**2)*delta_r
    !    !if(rid.eq.0)then
    !    !    print*, "r,phi_c,phi_0,phi_mbhc,phi_mbh0=",r,phi_c,phi_0,sppc%mbh_dmless/r,spp0%mbh_dmless/r
    !    !end if
	!end subroutine 
end subroutine

subroutine get_cluseter_with_mbh_adb_increase()
    use com_main_gw
    implicit none
    real(8) mbh_target,mbh_dm_target, dm, m0, m1
    integer nmod
    character*(4) tmpi
    integer i,ierr
    
    call get_system_dstr(bksams_arr_norm,.true.,ctl%gx_conv_cri)
    !call get_orbit_tables()
    call get_sample_para(dms,bksams_arr_norm,ctl%replace_sample_eceed_emax,spp_new) 
    call copy_spp_to_record(spp_new,spp_record(1))   
!    spp_record(1)=spp_new
    call update_sample_energy_indvd(ctl%replace_sample_eceed_emax)
    call update_arrays_single(.true.)

    if(rid.eq.0)then
        call output_dms_hdf5_pdf(dms, "./output/ini/hdf5/dms_0")
    end if

    n_record=1
    nmod=max(ctl%ini_adb_increase_nt/10,1)

    gx_func_max_step=-1d0; gx_func_min_step=-1d0
    do i=1, ctl%ini_adb_increase_nt
        if(ctl%adb_est_method.eq.adb_est_method_acc )then
            call get_sample_raq(bksams_arr_norm,spp_new)
        end if
        !print*, "bksams_arr_norm%sp(1)%rp:0=",bksams_arr_norm%sp(1)%rp
        spp_old=spp_new
        m1=10**((ctl%init_adb_mbh_log_factor_target-ctl%init_adb_mbh_log_factor)/dble(ctl%ini_adb_increase_nt-1)*dble(i-1)&
            +ctl%init_adb_mbh_log_factor)
        spp_new%mbh_dmless=m1
        spp_new%mbh=m1*m0_cl
        if(mod(i,nmod).eq.0) then
            n_record=n_record+1
            !spp_record(n_record)=spp_new
            call copy_spp_to_record(spp_new,spp_record(n_record)) 
        end if
        if(rid.eq.0)then
            print*, "i,m1=",i, spp_new%mbh_dmless
        end if

        call get_system_dstr(bksams_arr_norm,.true.,ctl%gx_conv_cri)

        call adb_cor_indvd_replace(bksams_arr_norm, spp_old,spp_new,ctl%adb_est_method)
                !print*, "bksams_arr_norm%sp(1)%rp:1=",bksams_arr_norm%sp(1)%rp
        !call get_orbit_tables()
        call get_sample_para(dms,bksams_arr_norm,ctl%replace_sample_eceed_emax,spp_new)   
        !print*, "bksams_arr_norm%sp(1)%rp:2=",bksams_arr_norm%sp(1)%rp
        !read(*,*)
        call mpi_barrier(mpi_comm_world,ierr) 
        call update_sample_energy_indvd(ctl%replace_sample_eceed_emax)
        call update_arrays_single(.true.)
        if(rid.eq.0)then
            write(unit=tmpi,fmt="(I4)") i
            call output_dms_hdf5_pdf(dms, "./output/ini/hdf5/dms_"//trim(adjustl(tmpi)))
        end if
    end do
        
    print*, "adb increase finished!"
end subroutine
subroutine get_system_dstr_adb_cor_one_time( )
    use com_main_gw
    implicit none
    integer n
    !type(s1d_type)::fmden_old, fphi_old, fmden_new, fphi_new, fmden_update,fphi_update
	integer ierr,i
    integer,parameter::nmax=20
    integer niter
    !type(diffuse_mspec)::dms_tmp
    real(8) cri,de_tmp, dmbh_tot, t1, t2
    type(s1d_type)::fma_new, fma_update
    if(rid.eq.0)then
        print*, "start of adb cor one time"
        call cpu_time(t1)
    end if 
    spp_old=spp_new
    call get_sample_para(dms,bksams_arr_norm,ctl%replace_sample_eceed_emax,spp_new)
    if(ctl%adb_est_method.eq.adb_est_method_acc)then
        call get_sample_raq(bksams_arr_norm,spp_new)
    end if
    call mpi_barrier(mpi_comm_world,ierr)
    if(ctl%enable_evl_mbh.ge.1)then
        !call get_sample_para(dms,bksams_arr_norm,.false.,spp_new)
        call get_mbh_increase(dmbh_tot)
        !dmbh_tot=600
        spp_new%mbh=spp_new%mbh+dmbh_tot
        spp_new%mbh_dmless=spp_new%mbh/m0_cl
        mbh_radius=spp_new%mbh/(my_unit_vel_c**2) 
    end if 
    call get_system_dstr(bksams_arr_norm,.true.,ctl%gx_conv_cri) 
    call adb_cor_indvd_replace(bksams_arr_norm, spp_old,spp_new,ctl%adb_est_method)
    
    if(rid.eq.0)then
        print*, "start get_orb_tables"
    end if
    call get_orbit_tables()

    call mpi_barrier(mpi_comm_world,ierr)
    if(rid.eq.0)then
        print*, "start get_sample_para_no_pd"
    end if
    call get_sample_para_no_pd(dms,bksams_arr_norm,ctl%replace_sample_eceed_emax,spp_new)
    
 
    if(rid.eq.0)then
        call cpu_time(t2)
        print*, "end of adb cor one time, used time=", t2-t1, " s"        
        
    end if
    !stop
end subroutine
  

subroutine get_sample_raq(sps,spp)
    use com_main_gw
    implicit none
    type(particle_samples_arr_type)::sps
    type(star_pot_para)::spp
    integer i
    real(8) jc_xy,ra_xy,rp_xy
    real(8) t1, t2
    if(rid.eq.0)then
        call cpu_time(t1)
    end if
    do i=1, sps%n
        associate(sp=>sps%sp(i))
            jc_xy=sp%jc/(r0_cl*ctl%v0)
            ra_xy=sp%ra/r0_cl
            rp_xy=sp%rp/r0_cl
            call get_radial_action(sp%x,sp%jm,jc_xy,spp,ra_xy,rp_xy,sp%raq)
        end associate
    end do
    if(rid.eq.0)then
        call cpu_time(t2)
        print*, "get_sample_raq time=", t2-t1, " s"
    end if
end subroutine


subroutine get_phi_diff(sppc, spp0, cri)
	use com_main_gw
	!type(s1d_type)::fc, f0
    type(star_pot_para)::spp0,sppc
	integer i
	real(8) fcy,f0y
    real(8) cri, xmin,xmax,xb(10)

	cri=0
    xmin=ctl%log10rmin_factor
    xmax=ctl%log10rmax_factor
    call set_range(xb,10,xmin,xmax,sts_type_dstr)
    !xmin=max(f0%xmin,fc%xmin)
    !xmax=min(f0%xmax,fc%xmax)
	do i=1, 10
    !    xb=(xmax-xmin)/dble(f0%nbin-1)*real(i-1)+xmin
!		call fc%get_value_l(xb,fcy)
!		call f0%get_value_l(xb,f0y)
        call get_phi_star_full_range(sppc,xb(i),fcy)
        call get_phi_star_full_range(spp0,xb(i),f0y)
		cri=cri+abs(fcy-f0y)
		!print*, "i, cri,dc=",i, cri, fcy, f0y
	end do
	cri=cri/10d0
	!print*, "cri=", cri
end subroutine
 
 
subroutine update_sample_energy_indvd(replace)
    use com_main_gw
    implicit none
    integer n,idx,out_flag_clone, n_chain_length
    integer nsam, i, j, i_dest,ierr 
    real(8),allocatable::en(:)
    integer,allocatable::id(:)
    type(chain_pointer_type),pointer::pt,ps
    logical,external:: Istransit,Invtransit
    integer nlvl,flag
    real(8),external:: rnd
    real(8) tmp, en0, en1
    integer amplifier, num_create_clone_bf
    logical condition,replace
    interface
        subroutine delete_some_samples(ps,pt,ch,flag)
            use com_main_gw
            implicit none
            type(chain_pointer_type),pointer::pt,ps
            type(chain_type)::ch
            integer flag
        end subroutine

	end interface


    if(allocated(en)) then
        deallocate(en, id)
    end if
    nsam=bksams_arr_norm%n
    allocate(en(nsam),id(nsam))

    do j=1, nsam
    !    en(j)=bksams_arr_norm%sp(j)%en
        id(j)=bksams_arr_norm%sp(j)%id
    end do

    ps=>bksams%head
    idx=1
    update_correction_emax=0
loop1: do while (associated(ps))
        pt=>ps%next
       if(ps%ob%id.eq.id(idx))then
            associate(sample=>ps%ob)
                !bksams_arr_norm%sp(idx)%en=ctl%energy_max-0.4d0
                if(spp_new%mbh_dmless.eq.0.and.replace)then
                    if(bksams_arr_norm%sp(idx)%en<ctl%energy_max)then
                    !    sample%exit_flag=exit_boundary_max
                    !    boundary_sts_emax_cros=boundary_sts_emax_cros+1
                        print*, "emax",bksams_arr_norm%sp(idx)%en,bksams_arr_norm%sp(idx)%x,bksams_arr_norm%sp(idx)%id
                        bksams_arr_norm%sp(idx)%en=ctl%energy_max*2-bksams_arr_norm%sp(idx)%en
                        bksams_arr_norm%sp(idx)%x=bksams_arr_norm%sp(idx)%en/ctl%energy0
                        call get_sample_para_one(dms,bksams_arr_norm%sp(idx),spp_new)
                        update_correction_emax=update_correction_emax+1
                    end if
                else
                    if(bksams_arr_norm%sp(idx)%en<ctl%energy_max)then
                        sample%exit_flag=exit_boundary_max
                        boundary_sts_emax_cros=boundary_sts_emax_cros+1 
                    end if
                end if

                if(bksams_arr_norm%sp(idx)%en>ctl%energy_boundary)then
                    sample%exit_flag=exit_boundary_min
                    print*, "out_boundary: rid,sp%e,emin",rid,bksams_arr_norm%sp(idx)%en,exit_boundary_min
                    ctl%mass_move_out_of_emin=ctl%mass_move_out_of_emin+sample%weight_real
                end if 
                en0=sample%en
                en1=bksams_arr_norm%sp(idx)%en
                sample%en=bksams_arr_norm%sp(idx)%en
                sample%x=bksams_arr_norm%sp(idx)%x
                sample%rp=bksams_arr_norm%sp(idx)%rp
                sample%ra=bksams_arr_norm%sp(idx)%ra
                sample%jc=bksams_arr_norm%sp(idx)%jc
                sample%period=bksams_arr_norm%sp(idx)%period
                sample%jm=bksams_arr_norm%sp(idx)%jm
                sample%jph=bksams_arr_norm%sp(idx)%jph
                sample%raq=bksams_arr_norm%sp(idx)%raq 
                    !
                condition=sample%exit_flag.eq.exit_boundary_min.or.sample%exit_flag.eq.exit_boundary_max
                if(ctl%clone_scheme.ge.1.and.(.not.condition))then
                    !ctl%chattery=4
                    call get_mass_idx(sample%m,sample_mass_idx)
                    if(sample_mass_idx.ne.-1)then
                        call clone_scheme(ps, en0, en1, ctl%clone_factor(sample_mass_idx),&
                            sample%exit_time, out_flag_clone)
                            !if(out_flag_clone.ge.1)then
                            !    particle_cloned=.true.
                            !end if
                        if(out_flag_clone.eq.100)then
                            sample%exit_flag=exit_invtransit
                            call delete_some_samples(ps,pt,bksams,flag)
                        end if
                    end if
                    !ctl%chattery=0
                end if
                !call sams_get_weight_clone_single_one(sample)
                !call get_sample_weight_real(sample)
                !call sams_get_weight_clone_single_one(bksams_arr_norm%sp(idx))
                !call get_sample_weight_real(bksams_arr_norm%sp(idx))
            end associate
            !print*, "idx,id=",idx, id(idx)
            if(idx<nsam)then
                idx=idx+1      
            else
                exit loop1
            end if

        end if
        !read(*,*)
        ps=>pt
    end do loop1     

    call collection_int(update_correction_emax)
    call collection_and_avg_real(ctl%mass_move_out_of_emin)
end subroutine
 
subroutine get_rc_given_jph(jph_xy,rc,spp)
    use md_star_pot
    use model_basic,only:ctl
    implicit none
    real(8) jph_xy,rc
    !real(8) mbh_dmless
    type(star_pot_para)::spp
    real(8) rtbis_yacc,par(50)
    integer ier,niter

    rc=10**rtbis_yacc(func,ctl%log10rmin_factor,ctl%log10rmax_factor,1d-9,&
     par,niter,1000,ier,.true.)
contains
    real(8) function func(x,par)
        implicit none
        real(8) x, par(50)
        real(8) beta_tmp,jc_xy
       ! print*, "x=",x
        call get_beta_full_range(spp,x,beta_tmp)
        jc_xy=((spp%mbh_dmless+beta_tmp)*10**x)**0.5
        !print*, "x,jc_xy,jph_xy=",x,jc_xy,jph_xy
        func=jc_xy-jph_xy
    end function    
end subroutine

subroutine get_maximum_ex_given_jph_fast(jph_xy ,ex,spp,mex)
    use com_sts_type
    use model_basic,only:emax_factor,emin_factor
    use md_star_pot
    implicit none
    real(8) jph_xy
    real(8) ex,  mex, rc
    type(star_pot_para)::spp
    real(8) beta_tmp,phi_tmp
    call get_rc_given_jph(jph_xy,rc,spp)
    call get_beta_full_range(spp,log10(rc),beta_tmp)
    call get_phi_star_full_range(spp,log10(rc),phi_tmp)
    mex=(spp%mbh_dmless+2*10**phi_tmp*rc-beta_tmp)/2d0/rc
end subroutine
 
subroutine get_dx_invariant_rad_action(ex,jm,jc,spp0,sppc,frphic,ra0,rp0,raq0,de)
    use com_sts_type
    use my_intgl
    use md_star_pot
    use model_basic,only:emin_factor,emax_factor
    implicit none
    !type(s1d_type)::fphi0,fphic
    type(s1d_type):: frphic!,fmac,frhoc
    type(star_pot_para)::spp0,sppc
    real(8) ex,jm,ra0,rp0,yout,jc,jlum,raq0,de,mex,exnew
    real(8) jph_xy
    real(8) rtbis_yacc,par(50), require_acc, xl, xh
    integer ier,niter

    jph_xy=jm*jc

    call get_maximum_ex_given_jph_fast(jph_xy,ex,sppc,mex)
    require_acc=1d-7; xl=max(ex-0.05,emin_factor); xh=min(ex+0.05,mex)
   ! print*, "xl,xh=",xl,xh
    exnew=rtbis_yacc(func,xl,xh,require_acc,par,niter,40,ier,.true.)
   ! print*, "niter=",niter
    !if(niter.ge.30)then
    !    print*, "exnew=",exnew
    !    require_acc=1d-6; xl=max(ex-0.05,emin_factor); xh=min(ex+0.05,mex)
    !    exnew=rtbis_yacc(func,xl,xh,require_acc,par,niter,1000,ier,.true.)
    !    print*, "exnew=", exnew
    !    read(*,*)
    !end if
    !print*, "raq1-raq0=",raq1,raq0,raq1-raq0
   ! print*, "exnew, ex=",exnew, ex
    de=exnew-ex
    if(ier<0)then
        print*, "exnew,ex,jm=",exnew,ex,jm
        stop
    end if
    !print*, "ex,exnew=",ex,exnew
    !stop
contains
    real(8) function func(x,par)
        implicit none
        real(8) x, par(50)
        real(8) jc_xy,rp_xy,ra_xy,pd_xy,raqc,jmnew
        !print*, "func:ex=",x
        call get_sample_para_one_xj_rpra(x,jm,jm*jc,sppc,frphic,jmnew,jc_xy,rp_xy,ra_xy)
        !print*, "jc_xy,rp_xy,ra_xy,pd_xy=",jc_xy,rp_xy,ra_xy,pd_xy
        call get_radial_action(x,jmnew,jc_xy,sppc,ra_xy,rp_xy, raqc)
        func=raqc-raq0
       ! print*, "x,raqc,raq0=",x,raqc,raq0,raqc-raq0
    end function
    subroutine test_func()
        implicit none
        type(s1d_type)::test
        integer i
        real(8) res,par(50)
        call test%init(xl,xh,20,sts_type_grid)
        call test%set_range()

        do i=1, 20
            test%fx(i)=func(test%xb(i),par)
        end do
        print*, "ex,jm,jc=",ex,jm,jc
        call test%print("test")
        stop
    end subroutine
end subroutine

subroutine get_radial_action(ex,jm,jc,spp,ra,rp,raq)
    use com_sts_type
    use my_intgl
    use md_star_pot
    use ieee_arithmetic
    use model_basic, only:pd_int_acc_a,pd_int_acc_r
    implicit none
    !type(s1d_type)::fphi 
    type(star_pot_para)::spp
    integer idid
    real(8) ex,jm,ra,rp,yout,raq,jc,jlum
    real(8) phi_tmp
    yout=0
    jlum=jm*jc
    !print*, "get_radial_action:ex,rp,ra=",ex,rp,ra
    call my_integral_acc(rp,ra,yout,pd_int_acc_a,pd_int_acc_r,fcn,idid)

    raq=yout*2
    if(idid<0)then
        print*, "rp,ra,ex,jm=",rp,ra,ex,jm
        print*, "spp%rho_min,phis,phis2=",spp%spt_rho_rmin,spp%phi_r1r2_s,spp%phi_r1r2_s2
        call get_phi_star_full_range(spp,log10(ra),phi_tmp)
        print*, "vr(ra)=",2*(10**phi_tmp+spp%mbh_dmless/ra-ex)-jlum**2/ra**2
        call get_phi_star_full_range(spp,log10(rp),phi_tmp)
        print*, "vr(rp)=",2*(10**phi_tmp+spp%mbh_dmless/rp-ex)-jlum**2/rp**2
        stop
    end if
contains
    subroutine fcn(n, x, y, f, par, ipar)
        use, intrinsic :: ieee_arithmetic
        implicit none
        integer n, ipar(100)
        real(8) x, y(n), f(n), par(100)
        real(8) logr,phi_out, ftmp
        logr=log10(x)
        call get_phi_star_full_range(spp,logr,phi_out)
        !print*,"rp", 2*(10**phi_out+mbh_dmless/x-ex)-jlum**2/x**2
        ftmp=2*(10**phi_out+spp%mbh_dmless/x-ex)-jlum**2/x**2
        if(ftmp<-1d-4)then
            print*, "error! in fcn of get_radial_action"
            print*, "x,f1**2=",x,2*(10**phi_out+spp%mbh_dmless/x-ex)-jlum**2/x**2
            print*, "rp,ra=",rp,ra
            stop
        end if
        f(1)=(abs(ftmp))**0.5
    end subroutine    
end subroutine