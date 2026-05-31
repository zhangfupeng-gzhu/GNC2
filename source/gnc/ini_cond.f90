subroutine get_ini_genx(spp)
    use com_main_gw
    implicit none
    integer i,j
    character*(4) tmpi
    real(8) xi,xisum,emax, emin
    type(star_pot_para)::spp
    integer nlvl,idx
    real(8) get_clone_deep,clone_factor

    ctl%ini_ge_tot%fx=0
    !xisum=0
    !do i=1, ctl%m_bins
    !    xisum=xisum+ctl%asymptot_ini(1,i)*ctl%bin_mass(i)
    !end do
    do i=1, ctl%m_bins
        !xi=ctl%asymptot_ini(1,i)*ctl%bin_mass(i)/xisum
        idx=ctl%ini_model_list_in(i)
        select case(ctl%ini_model_list(i))
        case(ini_den_model_dehnen)
           ! print*, "Dehnen,idx=",idx
            !print*, "dehnen%mtot=",ctl%dehnen(idx)%mtot
            call get_init_fe_gs_dehnen(ctl%bin_mass(i), dms%fr_phi, spp, &
                ctl%dehnen(idx),  emin_factor, spp%mbh_dmless, ctl%ini_ge(i))
        case(ini_den_model_plummer)
            call get_init_fe_gs_plummer(ctl%bin_mass(i), dms%fr_phi, spp, &
                ctl%plummer(idx), emin_factor, spp%mbh_dmless, ctl%ini_ge(i))
        case default 
            print*, "error! define ini gs model", ctl%ini_den_model
            stop
        end select
        ctl%ini_ge_tot%fx=ctl%ini_ge_tot%fx+ctl%ini_ge(i)%fx
       ! ctl%chattery=3
        if(ctl%chattery.ge.3.and.rid.eq.0)then
            write(unit=tmpi, fmt="(I4)") i
            call ctl%ini_ge(i)%print("ge"//trim(adjustl(tmpi)))
            !stop
        end if
    end do
    ! read(*,*)
    do i=1, ctl%m_bins
        
        ctl%ini_nx_log(i)=ctl%ini_nx_tot
        !ctl%ini_nx_log(i)%fx=log10(ctl%ini_nx_tot%fx*ctl%asymptot_ini(1,i))
        call get_nx(ctl%ini_nx_log(i),ctl%ini_ge(i), dms%fr_phi, spp%mbh_dmless, spp)
        
        do j=1, ctl%ini_nper_bin(i)%nbin
            nlvl=int(get_clone_deep(10**ctl%ini_nx_log(i)%xb(j),log_clone_bd_sep,  clone_e0_factor))
            clone_factor=ctl%clone_factor(i)
            ctl%ini_nper_bin(i)%fx(j)=10**ctl%ini_nx_log(i)%fx(j)*dble(clone_factor)**nlvl*dms%barp_ir%xsteps(j)
        end do
        call ctl%ini_nx_log(i)%prepare_spline()
    end do
    do i=1, ctl%m_bins
        ctl%ini_nx_tot%fx=ctl%ini_nx_tot%fx+10**ctl%ini_nx_log(i)%fx
    end do
    !ctl%ini_nx_tot%fx=log10(ctl%ini_nx_tot%fx)
    call ctl%ini_nx_tot%prepare_spline()
        
    if(ctl%chattery.ge.3.and.rid.eq.0)then
        call ctl%ini_nx_tot%print("nx_tot")
    end if
    ! print*, "xxxxxxxxxxxxxx"
    !ctl%debug=0
    !stop
end subroutine
subroutine get_spp_from_rho(spp)
    use md_star_pot
    use model_basic,only:ctl
    use MPI_comu,only:rid,mpi_comm_world
    implicit none
    type(s1d_type)::frho
    type(star_pot_para)::spp
    integer ier
    real(8) t1, t2
    if(rid.eq.0)then
        call cpu_time(t1)
    end if
    call check_spp_data(spp,ier)
    if(ier<0)then
        print*, "error in get_spp_from_rho"
        stop
    end if
    call get_spp_starpt(spp)
    call get_spt_phi_constants(spp)
    !if(rid.eq.0.and.ctl%chattery.ge.2)then
    if(rid.eq.0)then
        print*, "phi_r1r2_s,phi_r1r2_s2=", spp%phi_r1r2_s,spp%phi_r1r2_s2
        !if(rid.eq.0)then
        !    call spp_new%frho_star%print("0 rho")    
        !end if
        !!call mpi_barrier(mpi_comm_world,ier)
        !if(rid.eq.15)then
        !    call spp_new%frho_star%print("15 rho")    
        !end if
    end if
    !call mpi_barrier(mpi_comm_world,ier)

    call get_spp_fma(spp)
    !call set_fma_y2(spp)
    !call get_spp_beta(spp)
    call set_spt_y2(spp)
    if(rid.eq.0)then
        call cpu_time(t2)
        print*, "get_spp_from_rho used time:", t2-t1, " s"
        !call spp%fphi_star%print("fphi_star")
        !stop
    end if
end subroutine
subroutine get_all_star_fmden_init()
    use com_main_gw
    implicit none
    integer i
    dms%all%star%fmden%xb=ctl%ini_frho(1)%xb
    dms%all%star%fmden%fx=0
    do i=1, ctl%m_bins
        if(ctl%asymptot_ini(2,i)>0)then
            dms%all%star%fmden%fx=dms%all%star%fmden%fx+ctl%ini_frho(i)%fx*ctl%asymptot_ini(2,i)
            !print*, "i=",i,ctl%asymptot_ini(2,i)
        end if
    end do
    !if(rid.eq.0)then
        !print*, "warnning: Assuming the first mass bin is all star, if it isn't, change this line"
    !end if

endsubroutine
subroutine prepare_tables()
    use com_main_gw
    implicit none
    integer i,j
    real(8) rho_rmax,n_rmax
    if(rid.eq.0.and.ctl%init_adb_mbh_inc.ge.1)then
        print*, "init adb mbh inc, mbh_dmless=", spp_new%mbh_dmless
    end if
    !print*, "0"
    sample_logrmax=ctl%log10rmax_factor
    sample_logrmin=ctl%log10rmin_factor
    call init_diffuse_mspec_rtables(dms)
    call spp_new%init(ctl%log10rmin_factor,ctl%log10rmax_factor,ctl%dstr_bins_r)
 
	call init_stellar_obj_rtables(dms)
    if(rid.eq.0.and.ctl%chattery.ge.2)then
        print*, "start prepare tables"
    end if
    call init_ctl_rtables()
    !call ctl%ini_frho_tot%print("ini frho")
    !print*, "1"
    !call get_ini_dehnen_pot(fphi_theory)
    call get_all_star_fmden_init()
    
    dms%ALL%ALL%fmden=ctl%ini_frho_tot
    do i=1, dms%n
        dms%mb(i)%all%fna=ctl%ini_fna(i)
        do j=1, n_tot_comp_sg
            dms%mb(i)%dsp(j)%p%fden%fx=ctl%ini_frho(i)%fx/ctl%bin_mass(i)*ctl%asymptot_ini(j+1,i)
        end do
    end do

    spp_new%frho_star=ctl%ini_frho_tot
    spp_new%has_set_density=.true.
    call ctl%ini_frho_tot%get_value_s(sample_logrmin,spp_new%spt_rho_rmin)
    !spp_new%spt_rho_rmin=10**spp_new%spt_rho_rmin
    spp_new%has_set_rhomin=.true.
    
    !do i=1, dms%n
    !    call ctl%ini_frho(i)%get_value_s(sample_logrmax,rho_rmax)
    !    n_rmax=rho_rmax/ctl%bin_mass(i)
    !    do j=1, n_tot_comp_sg
    !        dms%mb(i)%dsp(j)%p%spt_rho_rmax=n_rmax*ctl%asymptot_ini(j+1,i)
    !        !print*, "i,j,rho_rmax=",i,j,rho_rmax,n_rmax,dms%mb(i)%dsp(j)%p%spt_rho_rmax
    !       ! dms%mb(i)%dsp(j)%p%spt_rho_rmax=0
    !    end do
    !end do
    !read(*,*)
    call get_spp_from_rho(spp_new)
    

    if(ctl%chattery.ge.2.and.rid.eq.0)then
        call ctl%ini_frho_tot%print("ctl%ini_frho_tot")
        call spp_new%fphi_star%print("fphi_star")
    endif  
    if(ctl%chattery.ge.2.and.rid.eq.0)then
		print*, "rho_min,M_within_rmax,phi_r1r2_s,phi_r1r2_s2=", &
            spp_new%spt_rho_rmin, spp_new%M_r_within_max,spp_new%phi_r1r2_s,spp_new%phi_r1r2_s2
	end if

    call get_diffuse_mspec_rtables(dms,spp_new)
    
    call get_ini_ebounds()
    call get_ini_sample_ebound()
    call get_diffuse_mspec_ebound()

	call init_diffuse_mspec_etables(dms)

	call set_mass_bin_mass_given(dms, ctl%bin_mass, ctl%bin_mass_m1,&
            ctl%bin_mass_m2, ctl%asymptot_ini, ctl%m_bins)
    if(rid.eq.0)then
        print*, "start etables"
    end if
    !call get_dms_dlnx()
    call get_frphi(spp_new%fphi_star,dms%fr_phi)
    !print*, "chattery=",ctl%chattery

    sample_nxgx_logemin=max(sample_logemin,log10(dms%emin))
    sample_nxgx_logemax=min(log10(dms%emax),sample_logemax)

    !ctl%barge_grid_type=barge_grid_type_iregular_phi
    
    !call set_barp_step_size(dms,dms%barp_ir,spp_new)

    call dms%barp_ir%init(sample_nxgx_logemin,sample_nxgx_logemax,&
    dms%dstr_bins_e,sts_type_dstr)

    !call set_nx_ranges_xb_ir(dms%barp_ir%xb,dms%barp_ir%xsteps, dms%barp_ir%xmin, &
    !dms%barp_ir%xmax,dms%barp_ir%nbin,dms%barp_ir%bin_type,dms,spp_new,ctl%barge_grid_type)
    call set_nx_ranges_xb_ir(dms%barp_ir%xb,dms%barp_ir%xsteps, dms%barp_ir%xmin, &
    dms%barp_ir%xmax,dms%barp_ir%nbin,dms%barp_ir%bin_type,dms,spp_new,&
        ctl%barge_grid_type)

    call get_barp_ir(dms,dms%barp_ir,spp_new)
    if(rid.eq.0.and.ctl%chattery.ge.2)then
        call dms%barp_ir%print("ini cond:barp_ir")
    end if
    !if(rid.eq.0.and.ctl%chattery.ge.2)then
    !    if(ctl%barge_grid_type.eq.barge_grid_type_iregular)then
    !        call dms%barp_ir%print("barp_ir")
    !    else
    !        call dms%barp%print("barp")
    !    end if
    !end if
    !read(*,*)
    !stop
    !if(ctl%barge_grid_type.ne.barge_grid_type_regular)then
    !    call set_etable_ir_bins()
    !end if

    !call get_init_etables()
    !print*, "1"
    call get_orbit_tables()
    
    call set_ini_nx_ranges()

    call get_ini_genx(spp_new)
    !if(rid.eq.0)then
    !    call dms%all%all%fmden%print("fmden0.9")
    !end if
    !print*, "end"
    !stop
    !if(ctl%chattery.ge.2)then
    !ctl%debug=0
    !end if
endsubroutine
subroutine get_ini_sample_ebound()
    use com_main_gw
    implicit none
    real(8) emax,emin
    !if(mbh_dmless.eq.0)then
        sample_logemin=log10emin_factor
        sample_logemax=log10emax_factor!-abs(log10emax_factor)*0.01
        sample_emin=10**sample_logemin
        sample_emax=10**sample_logemax
    !else
    !    call get_phi_star_full_range(dms%fphi_star,ctl%log10rmin_factor,emax,spp_new)
    !    call get_phi_star_full_range(dms%fphi_star,ctl%log10rmax_factor,emin,spp_new)
    !    sample_logemin=max(log10(10**emin+mbh_dmless/10**ctl%log10rmax_factor),ctl%loge_min_factor)
    !    sample_logemax=log10(10**emax+mbh_dmless/10**ctl%log10rmin_factor)
    !    !print*, "emin,emax,rmin,rmax=",emin,emax,ctl%log10rmin_factor,ctl%log10rmax_factor
    !end if
    if(rid.eq.0)then
        print*, "init: sample_logemin,emax=",sample_logemin,sample_logemax
    end if
    !call set_edstr_bound()
end subroutine
!subroutine set_edstr_bound()
!    use com_main_gw
!    implicit none
!
!    emax_dstr_factor=min(10**sample_logemax,emax_factor-abs(emax_factor)*0.01)
!    emin_dstr_factor=10**sample_logemin
!    if(rid.eq.0)then
!        print*, "dstemin,dstremax=",log10(emin_dstr_factor),log10(emax_dstr_factor)
!    end if
!end subroutine
subroutine set_ini_nx_ranges()
    use com_main_gw
    implicit none
    real(8) xsteps(dms%dstr_bins_e)
    integer i
    real(8) emax, emin
    type(s1d_type)::barp
    select case(ctl%ebin_type)
    case(ebin_type_log)
        emin=log10emin_factor
        emax=log10emax_factor 
    end select
    
    call ctl%ini_nx_tot%init(emin,emax,ctl%dstr_bins_e,sts_type_dstr)
    call ctl%ini_ge_tot%init(emin,emax,ctl%dstr_bins_e,sts_type_dstr)
    do i=1, ctl%m_bins
        call ctl%ini_ge(i)%init(emin,emax,ctl%dstr_bins_e,sts_type_dstr)
        !print*, "ctl%ini_ge(i)%nbin=", ctl%ini_ge(i)%nbin, size(ctl%ini_ge(i)%xb)
    end do
    !select case(ctl%barge_grid_type)
    !case(barge_grid_type_iregular_phi,barge_grid_type_iregular_barp,barge_grid_type_regular)

        associate(b=>ctl%ini_nx_tot,d=>ctl%ini_ge_tot)
            if(dms%barp_ir%nbin.ne.b%nbin)then
                print*, "error! nbin size not match"
                stop
            end if
            b%xb=dms%barp_ir%xb
            d%xb=b%xb
        end associate

        do i=1, ctl%m_bins
            associate(b=>ctl%ini_ge(i))
                b%xb=dms%barp_ir%xb
            end associate
        end do
    !case(barge_grid_type_regular)
    !    call ctl%ini_nx_tot%set_range()
    !    call ctl%ini_ge_tot%set_range()
    !    do i=1, ctl%m_bins
    !        call ctl%ini_ge(i)%set_range()
    !    end do
    !end select
    do i=1, ctl%m_bins
        ctl%ini_nper_bin(i)=ctl%ini_nx_tot
    end do
end subroutine


real(8) function gen_ran_from_dstr_consider_clone(x,y,y2,n, xmin,xmax,ymax, clone_e0, amplifier, int_type)
	use constant
    use model_basic,only:ctl,log_clone_bd_sep
	implicit none
	real(8) ranout
	integer n,int_type
	real(8) x(n),y(n), rnd, ymax, get_clone_deep
	real(8) xmin,xmax, xtmp,ytmp, fy,yp1,ypn, clone_e0
	real(8) y2(n)
    integer nlvl, amplifier,ier

	!allocate(y2(n))
	if(xmin>xmax)then
        print*, "error in gen_ran_from_data"
        read(*,*)
    end if
	!print*, "xmin,xmax=", xmin, xmax
	!print*, "ymax=", ymax
100	xtmp=rnd(xmin,xmax)
	ytmp=rnd(0d0,ymax)
	select case(int_type)
	case(1)
		!yp1=1d30;ypn=1d30
		!call spline_mylib(x,y,n,yp1,ypn,y2)
		call splint_mylib(x,log10(y),y2,n,xtmp,fy,ier)
	case(2)
		call linear_int(x,log10(y),n,xtmp,fy)
    case default
		print*, 'error, define int_type'
        read(*,*)
	end select
        fy=10**fy
        nlvl=int(get_clone_deep(10**xtmp,log_clone_bd_sep,  clone_e0))
        !print*, 'clone_e0=',clone_e0
        !read(*,*)
        !print*, "10**xtmp, clone_e0, nlvl, fy=", 10**xtmp, clone_e0, nlvl, fy
        if(nlvl>=1)then
            fy=fy*dble(amplifier)**dble(nlvl)
        end if
        !print*, "amplifier,fy=", amplifier, fy
        !read(*,*)
	    !print*, "xtmp, fy=",xtmp,fy
        !print*, "xtmp, ytmp, fy=",xtmp, ytmp, fy,amplifier, clone_e0
	if(ytmp>fy) then
		goto 100
	end if
	gen_ran_from_dstr_consider_clone=xtmp
    !print*, "gen_ran_from_dstr_consider_clone:xtmp=",xtmp
    !if(ctl%debug.ge.1)then
    !    if(xtmp>0.92)then
    !        print*, "xtmp=",xtmp
    !        print*, "fy=",fy
    !        print*, "ytmp=",ytmp
    !        print*, "error!"
    !        stop
    !    end if
    !end if
    !print*, "gen_ran_from_data=",gen_ran_from_data
!	pause	
end function

real(8) function gen_ran_from_dstr_consider_clone_cum(x,y,n, xmin,xmax)
	use constant
    use model_basic,only:ctl,log_clone_bd_sep
	implicit none
	real(8) ranout
	integer n,i
	real(8) x(n),y(n), rnd, ymax, get_clone_deep
	real(8) xmin,xmax, xtmp,ytmp, fx,yp1,ypn, clone_e0
    integer nlvl, amplifier,ier

	!allocate(y2(n))
	if(xmin>xmax)then
        print*, "error in gen_ran_from_data"
        read(*,*)
    end if
	

	ytmp=rnd(0d0,1d0)

    call linear_int(y,x,n,ytmp,fx)
    !print*, ytmp, fx
	gen_ran_from_dstr_consider_clone_cum=fx
    !read(*,*)
end function

subroutine get_orbit_tables()
    !use md_dms
    use model_basic!,only:ctl
    use md_star_pot
    use MPI_comu 
    implicit none
	real(8) t1, t2
	integer i,ierr
    if(rid.eq.0)then
        call cpu_time(t1)
    end if
    !ctl%debug=ctl%chattery
    !if(ctl%barge_grid_type.ne.barge_grid_type_regular)then
        !print*, "1"
        !if(rid.eq.0)then
        !    print*, "start set_etable_ir_bins"
        !end if
        call set_etable_ir_bins()
    !end if
    call mpi_barrier(mpi_comm_world,ierr)
    !if(rid.eq.0)then
   !     print*, "start get_xrc"
    !end if

    call get_xrc(dms%frc_x,spp_new)
	!call cpu_time(t1)
	!do i=1, 100
    !call mpi_barrier(mpi_comm_world,ierr)
    !if(rid.eq.0)then
    !    print*, "start get_rc"
    !end if
    call get_rc(dms%rc,spp_new)
	!end do
	!call cpu_time(t2)
	!print*, "rc,time=",t2-t1, (t2-t1)/100d0/dms%rc%nbin

    if(ctl%chattery.ge.2.and.rid.eq.0)then
        call dms%rc%print("rc")
    end if
	!call cpu_time(t1)
	!do i=1, 100
    !if(rid.eq.0)then
    !    print*, "start get_jc"
    !end if
    !call mpi_barrier(mpi_comm_world,ierr)

    call get_jc_dmless(dms%jc,dms%rc, spp_new)
	!end do
	!call cpu_time(t2)
	!print*, "jc,time=",t2-t1, (t2-t1)/100d0/dms%jc%nbin

    if(ctl%fden_ana_est_method.eq.fden_ana_est_method_2d)then
        call get_jc_dmless_sample_erange(dms%jc_sample_erange, spp_new)
        !call dms%jc_sample_erange%print("jc_sample_erange")
    end if
    if(ctl%chattery.ge.2.and.rid.eq.0)then
        call dms%jc%print("jcdm")
        call dms%fr_phi%print("fr_phi")
    end if 
    !end if
    call mpi_barrier(mpi_comm_world,ierr) 
    call get_rp_ra_dm_mpi(spp_new, dms%jc, dms%rc, dms%fr_phi, dms%rp, dms%ra)
    !stop
  	!call cpu_time(t2)
	!print*, "rpra,time=",t2-t1, (t2-t1)/dble(dms%rp%nx)/dble(dms%rp%ny)
	!read(*,*)
    if(ctl%chattery.ge.3.and.rid.eq.0)then
        call dms%rp%print("rp")
    end if
    !if(rid.eq.0)then
    !    print*, "start get_pd"
    !end if
    call mpi_barrier(mpi_comm_world,ierr)

    !call get_pd_dmless(spp_new, dms%jc, dms%rp, dms%ra, dms%pd)
    !print*, "1:dms%pd%fx(10,10:80)=", dms%pd%fxy(10,10:80)
    !print*, "1:dms%pd%fx(60,10:80)=", dms%pd%fxy(60,10:80)
    call get_pd_dmless_mpi(spp_new, dms%jc, dms%rp, dms%ra, dms%pd)
    !print*, "2:rid=",rid, dms%pd%fxy(10,10:80)
    !print*, "3:rid=",rid, dms%pd%fxy(60,10:80)
    !stop
    if(ctl%chattery.ge.3.and.rid.eq.0)then
        call dms%pd%print("pd")
        print*, "pd finished"
    end if     
    !ctl%debug=0
    !stop 
    !if(rid.eq.0)then
    !    print*, "start get surface density"
    !    call dms%all%star%fmden%print("fmden")
    !end if

    !call mpi_barrier(mpi_comm_world,ierr)

    call get_surface_density(dms%all%star%fmden,dms%surface_den)
    !if(rid.eq.0)then
        !call dms%surface_den%print("surface_den")
    !    print*, "start get_cumulative_surface_density"
    !!end if
    !call mpi_barrier(mpi_comm_world,ierr)

    call get_cumulative_surface_density(dms%surface_den,dms%cum_sur_den)
    !if(rid.eq.0)then
        !call dms%cum_sur_den%print("cum_sur_den")
    !    print*, "end of orbit table"
    !end if
    !call dms%all%star%fmden%print("fmden")
    !call dms%surface_den%print("sur")
    !call dms%cum_sur_den%print("cum_sur")
    !read(*,*)
    call get_reff_now(nsc_radius_eff)
    if(rid.eq.0)then
        call cpu_time(t2)
        print*, "orbit table time=",t2-t1
    end if


end subroutine
 
subroutine get_surface_density(frho,sur_den)
    use com_main_gw
    implicit none
    type(s1d_type),intent(in)::frho
    type(s1d_type)::sur_den
    real(8) rev, rmax,rmin
    integer i,idid
    real(8) int_out

    call sur_den%init(dms%logrmin,dms%logrmax,dms%dstr_bins_r,sts_type_dstr)
    call sur_den%set_range()
    rmin=10**dms%logrmin
    rmax=10**dms%logrmax
    do i=1, sur_den%nbin
        int_out=0
        rev=10**sur_den%xb(i)
        call my_integral_acc(0d0,pi/2d0, int_out,&
            fr_funcs_int_acc_a,fr_funcs_int_acc_r, fcn,idid)
        sur_den%fx(i)=int_out*rev*2d0
        if(idid<0)then
            print*, "in get_surf_density"
            print*, "rmin,rmax=",rmin,rmax
            call frho%print("frho")
        end if
    end do

contains 
    subroutine fcn(n, x, y, f, par, ipar)
        implicit none
		integer n, ipar(100)
		real(8) x, y(n), f(n), par(100)
        real(8) r, rhov
        r=rev/cos(x)
        if(r>rmin.and.r<rmax)then
            !call get_rho_full_range_spp(spp_new,log10(r),rhov)
           ! call frho%get_value_s(log10(r),rhov)
           ! f(1)=10**rhov/(cos(x)**2)
            call frho%get_value_s(log10(r),rhov)
            f(1)=rhov/(cos(x)**2)
            if(rhov<0)f(1)=0
        else
            f(1)=0d0
        end if
        
    end subroutine
end subroutine


subroutine get_cumulative_surface_density(sur_den,cum_sur_den)
    use com_main_gw
    implicit none
    type(s1d_type),intent(in)::sur_den
    type(s1d_type)::cum_sur_den
    real(8) rev,rmin,rmax
    integer i,idid
    real(8) int_out

    call cum_sur_den%init(dms%logrmin,dms%logrmax,dms%dstr_bins_r,sts_type_dstr)
    call cum_sur_den%set_range()
    !rmin=dms%logrmin
    !rmax=dms%logrmax
    do i=1, cum_sur_den%nbin
        int_out=0
        call my_integral_acc(dms%logrmin-2d0,cum_sur_den%xb(i), int_out,&
            fr_funcs_int_acc_a,fr_funcs_int_acc_r, fcn,idid)
            cum_sur_den%fx(i)=int_out*2*pi*log(10d0)
        if(idid<0)then
            print*, "in get_cumulative_surface_density"
            call sur_den%print("sur_den")
        end if
    end do

contains 
    subroutine fcn(n, x, y, f, par, ipar)
        implicit none
		integer n, ipar(100)
		real(8) x, y(n), f(n), par(100)
        real(8) r, rhov
        !print*, "x,rmin,rmax=",x,rmin,rmax
        if(x>dms%logrmin)then
            call sur_den%get_value_s(x,rhov)
        else
            rhov=sur_den%fx(1)
        end if
        f(1)=rhov*10**(x*2)
        !else
        !    f(1)=0d0
        !end if
        
    end subroutine
end subroutine

subroutine get_jc_dmless_sample_erange(jc,spp)
    use com_main_gw
    use md_star_pot
    implicit none
    type(s1d_type)::jc
    type(star_pot_para)::spp
    real(8) rc, jc_dmless, r_c_iter
    integer i,ier
    real(8) ex

    do i=1, jc%nbin
        ex=10**jc%xb(i)
        !print*, "ex, 10emax=",ex,log10(dms%emax)
        if(ex<dms%emax)then
            rc=r_c_iter(spp, ex,ier)
            jc%fx(i)=jc_dmless(rc,spp)
        else
            rc=0
            jc%fx(i)=0
        end if
    end do
    !if(dms%mb(1)%all%gxj_ir%nx>0)then
    !    print*, jc%xb
    !    print*, dms%mb(1)%all%gxj_ir%xcenter
    !    read(*,*)
    !end if
end subroutine
subroutine get_ini_self_consistent_solution()
    use com_main_gw
    implicit none
    call update_arrays_single(.true.)
	if(ctl%barge_evl_method.eq.barge_evl_method_direct)then
		call get_sample_para(dms,bksams_arr_norm,ctl%replace_sample_eceed_emax,spp_new)
	end if
	!print*, "what"
	!stop
   ! print*, "XXXXXXXXXXXXXX",RID
	if(ctl%init_adb_mbh_inc.ge.1)then
		call get_cluseter_with_mbh_adb_increase()
		call get_orbit_tables()
	else
		call get_system_dstr(bksams_arr_norm,.true.,ctl%gx_conv_cri)
		call get_orbit_tables()
	end if
	!if(rid.eq.0)then
	!	call dms%fphi_star%print("fphi_star1")
	!	call dms%all%all%fmden%print("fmden")
	!end if
    if(rid.eq.0)then
        print*, "start get sample_para"
    end if
    call get_sample_para(dms,bksams_arr_norm,ctl%replace_sample_eceed_emax,spp_new) 
	call update_sample_energy_indvd(ctl%replace_sample_eceed_emax)
end subroutine
 
subroutine get_diffuse_mspec_ebound()
    !use md_dms
    use model_basic
    use MPI_comu,only:rid
    use md_star_pot
    implicit none

    dms%emin=emin_factor; 
    if(spp_new%mbh_dmless.eq.0)then
        dms%emax=emax_factor-abs(emax_factor)*1d-9
    else
        dms%emax=emax_factor
    end if
    if(rid.eq.0)then
        print*, "dmsemin,dmsemax=",log10(dms%emin),log10(dms%emax)
    end if



    !dms%emin=max(10**sample_logemin,dms%emin)
    !dms%emax=min(dms%emax,10**sample_logemax)
end subroutine

subroutine get_ini_ebounds()
    use model_basic
    use md_star_pot
    use MPI_comu,only:rid
    implicit none
    real(8) emax,emin, logmbh_hz_emax,mbh_hz_emax,mbh_radius_dmless
    !call spp_new%fphi_star%print("fphi_star")
    call get_phi_star_full_range(spp_new,ctl%log10rmin_factor,emax)
    true_log10emax_factor=10**emax+spp_new%mbh_dmless/10**ctl%log10rmin_factor

    if(spp_new%mbh_dmless.ne.0)then
        mbh_radius=spp_new%mbh/(my_unit_vel_c**2)
        mbh_radius_dmless=mbh_radius/r0_cl
        call get_phi_star_full_range(spp_new,log10(mbh_radius_dmless),logmbh_hz_emax)
        mbh_hz_emax=10**logmbh_hz_emax+spp_new%mbh_dmless/mbh_radius_dmless
        emax_factor=min(true_log10emax_factor,ctl%emax_boundary,mbh_hz_emax)
    else
        emax_factor=min(true_log10emax_factor,ctl%emax_boundary)
    end if
    !if(ctl%emax_boundary.ge.1)then
    !    emax_factor=min(true_log10emax_factor,1d5)
    !else
    !    emax_factor=true_log10emax_factor
    !end if

    log10emax_factor=log10(emax_factor)
	true_log10emax_factor=log10(true_log10emax_factor)

    call get_phi_star_full_range(spp_new,ctl%log10rmax_factor,emin)
    emin_factor=10**emin+spp_new%mbh_dmless/10**ctl%log10rmax_factor
    log10emin_factor=log10(emin_factor)

    if(log10emin_factor.ge.1.or.log10emin_factor.ge.log10emax_factor)then
        print*, "phiemin, mbhe=",10**emin, spp_new%mbh_dmless/10**ctl%log10rmax_factor
        print*, "log10emin_factor=",log10emin_factor,log10emin_factor.ge.1
        print*, "log10emax_factor=",log10emax_factor, log10emin_factor.ge.log10emax_factor
        !call spp_new%fphi_star%print("fphi_star")
        print*, "rho_min,phi_r1r2_s,phi_r1r2_s2=",spp_new%spt_rho_rmin,&
            spp_new%phi_r1r2_s,spp_new%phi_r1r2_s2
        print*, "log10emin_factor, true_log10emax_factor=", log10emax_factor, true_log10emax_factor,mbh_hz_emax

    end if

    !log10emin_factor=log10(dms%emin)
    !log10emax_factor=log10(dms%emax)

    ctl%energy_min=ctl%energy0*emin_factor; 
    ctl%energy_max=ctl%energy0*emax_factor
    ctl%x_boundary=emin_factor
    !ctl%x_boundary=0d0
    ctl%energy_boundary=ctl%x_boundary*ctl%energy0
    !print*, "logemin,logemax=",dms%logemin, dms%logemax
    if(ctl%clone_scheme.ge.1)then
        clone_emax=ctl%energy_max/ctl%clone_e0
        log10clone_emax=log(clone_emax)/log(clone_bd_sep)+1
    end if
    if(rid.eq.0)then
        select case(ctl%ebin_type)
        case(ebin_type_log)
            if(spp_new%mbh_dmless.ne.0)then
                print*, "logemin,logemax,logmbh_hz_emax=",log10emin_factor,log10emax_factor,log10(mbh_hz_emax)
            else
                print*, "logemin,logemax=",log10emin_factor,log10emax_factor
            end if 
        end select
        print*, "logrmin,logrmax=",dms%logrmin, dms%logrmax 
    end if
end subroutine
subroutine get_init_fe_gs_dehnen(bin_mass, rmax,spp,dehnen,  xmin, mbhin, fe)
	use com_sts_type
	use my_intgl
	use constant
	use model_basic,only:ctl,type_mdehnen
	use md_coeff
    use md_star_pot
	implicit none
	type(s1d_type)::rmax,fe!, frho
    type(type_mdehnen)::dehnen
    type(star_pot_para)::spp
	real(8) xmin,xmax, logrmin,logrmax, bin_mass,rho_rmin
	real(8) gamma, ra, enx,phi_min,phi_max,mtot
	real(8) drhodr,drhodr2,dphidr,dphidr2, phi_tmp
	real(8) rho_tmp, rho_tmp0,beta_tmp,  gamma_tmp, gamma_tmp2,fout
	real(8) rho_tmp_max0,beta_tmp_max,phi_tmp_max
	real(8) drhodr_max,dphidr_max, mbhin,radius_max
	integer i,ier, n

	!if(fe%xb(1)>phi_tot%fx(1))then
	!	print*, "error! fe%xb(1)>phi_tot%fx(1),check!"
	!	stop
	!endif
	n=spp%frho_star%nbin
	!logrmax=rho%xb(n)
	!rho_tmp_max=rho%fx(n)
	!beta_tmp_max=fma%fx(n)
    gamma=dehnen%gamma
    ra=dehnen%ra_crit
    mtot=dehnen%mtot

    logrmax=spp%frho_star%xmax
    radius_max=10**logrmax 
    call get_ini_dehnen_dens_func(rho_tmp_max0,radius_max,mtot,ra,gamma) 
    call get_beta_full_range(spp,logrmax,beta_tmp_max)

	!phi_tmp_max=phi_tot%fx(n)
	dphidr_max=-(mbhin+beta_tmp_max)/(radius_max)**2
	drhodr_max=-rho_tmp_max0*(gamma/radius_max+(4-gamma)/(radius_max+ra))
	!print*, "max:rho,beta,phi=",rho_tmp_max,beta_tmp_max, phi_tmp_max
	!print*, "max:dphi_dr,drho_dr=", dphidr_max,drhodr_max
    !print*, fe%nbin, size(fe%xb)

	do i=1, fe%nbin
		fout=0
		select case(ctl%ebin_type)
		case(ebin_type_log)
			enx=10**fe%xb(i) 
		end select 
		if(xmin<enx)then
            ier=-99
			call my_gs_integral(xmin,enx, fout, 0d0,-0.5d0,1d-10,1d-10,1,fx,ier)
        !call my_gs_integral(1d-5,enx, fout, 0d0,-0.5d0,1d-10,1d-10,1,fx,ier)
		else
			fout=0
            return
		end if
		
		fe%fx(i)=(fout +(enx-xmin)**(-0.5)*drhodr_max/dphidr_max)/pi**0.5
        !fe%fx(i)=fout! +(enx-xmin)**(-0.5)*drhodr_max/dphidr_max)/pi**0.5
        !print*, fout, drhodr_max,dphidr_max
        !fe%fx(i)=fout/pi**0.5
		if(fe%fx(i)<0d0)then
			print*, "get_init_fe_gs_dehnen:xmin, en, fe%fx(i)=", xmin, enx, fe%fx(i)
			fe%fx(i)=0
			!read(*,*)
		end if
	end do
	fe%fx=fe%fx/bin_mass
    !read(*,*)
contains 
	real(8) function fx(x)
		use, intrinsic :: ieee_arithmetic
		implicit none
		real(8) x, r,phi,logr

		phi=x
		!call phi_tot_sorted%print("phi_tot_sorted")
		!call phi_tot_sorted%get_value_s_y(log10(phi),logr)
		call get_rmax_accurate(spp,rmax,log10(phi),logr)
		r=10**logr
		!call rmax%print("rmax")
		!print*, "logphi,logr, r=", log10(phi), logr, r

		!call rho%get_value_s(logr,rho_tmp)
        call get_rho_full_range_spp(spp,logr,rho_tmp) 
        call get_ini_dehnen_dens_func(rho_tmp0,r,mtot,ra,gamma)

		call get_beta_full_range(spp,logr, beta_tmp)
		!call beta%get_value_s(logr,beta_tmp)
		
		dphidr=-(mbhin+beta_tmp)/r**2
		dphidr2=(2*(mbhin+beta_tmp)-4*pi*r**3*rho_tmp)/r**3
		gamma_tmp=(gamma/r+(4-gamma)/(r+ra))
		gamma_tmp2=(gamma/r**2+(4-gamma)/(r+ra)**2)

		drhodr=-rho_tmp0*gamma_tmp
		drhodr2=rho_tmp0*gamma_tmp**2+rho_tmp0*gamma_tmp2

		fx=(drhodr2/dphidr**2-drhodr*dphidr2/dphidr**3)
       ! if(fx<0) fx=0
		!print*, "logphi,logr, r=", log10(phi), logr, r, fx
		!print*, "fx=", fx
		!stop
		!if(fx<0)then
		!	print*, "x,logx,rho,beta,r=",x,log10(x),rho_tmp,beta_tmp, logr
		!	print*, "drhodr2,dphidr,drhodr,dphidr2=",drhodr2,dphidr,drhodr,dphidr2
		!	print*, "drhodr2*dphidr-drhodr*dphidr2=",drhodr2*dphidr-drhodr*dphidr2
		!	print*, "fx=",(drhodr2*dphidr-drhodr*dphidr2)/dphidr**3
		!	print*, "fx=",fx
		!	read(*,*)
		!end if
	end function
end subroutine

subroutine get_init_fe_gs_plummer(bin_mass, rmax,spp,plummer, xmin, mbhin, fe)
    use ieee_arithmetic
	use com_sts_type
	use my_intgl
	use constant
	use model_basic,only:ctl,type_mplummer
	use md_coeff
    use md_star_pot
	implicit none
	type(s1d_type)::rmax,fe
    type(star_pot_para)::spp
    type(type_mplummer)::plummer
	real(8) xmin,xmax, logrmin,logrmax, bin_mass
	real(8) gamma, ra, enx,phi_min,phi_max,mtot
	real(8) drhodr,drhodr2,dphidr,dphidr2, phi_tmp
	real(8) rho_tmp, beta_tmp,  gamma_tmp, gamma_tmp2,fout
	real(8) beta_tmp_max,phi_tmp_max,rho_tmp0
	real(8) drhodr_max,dphidr_max, mbhin, radius_max,rho_tmp_max0
	integer i,ier, n

	ra=plummer%ra_crit
    mtot=plummer%mtot

	n=spp%frho_star%nbin
	logrmax=spp%frho_star%xmax
	radius_max=10**logrmax
	!call get_rho_plummer_outside(logrmax,rho_tmp_max)
	!rho_tmp_max=rho%fx(n)
    !call rho%get_value_s(logrmax,rho_tmp_max)
    !call get_rho_full_range_spp(spp,logrmax,rho_tmp_max)
    call get_plummer_den(mtot,ra,logrmax,rho_tmp_max0)
    call get_beta_full_range(spp,logrmax,beta_tmp_max)
	!beta_tmp_max=fma%fx(n)
	!print*, "beta_tmp_max=",beta_tmp_max
	!stop
	!phi_tmp_max=phi_tot%fx(n)
	dphidr_max=-(mbhin+beta_tmp_max)/radius_max**2
	drhodr_max=-rho_tmp_max0*5*radius_max/(radius_max**2+ra**2)
    !print*, "logrmax,radius_max=", logrmax, radius_max
    !print*, "rho_tmp_max,beta_tmp_max,dphidr_max, drhodr_max=", &
    !    rho_tmp_max,beta_tmp_max,dphidr_max, drhodr_max

	!print*, drhodr_max/dphidr_max, 3*ra**2/4d0/pi/ctl%plummer_model_mtot**4*&
	!	(ctl%plummer_model_mtot/(radius_max**2+ra**2)**0.5)**4*5
	!read(*,*)
	!print*, "beta_tmp_max=",beta_tmp_max

	do i=1, fe%nbin
		fout=0
		select case(ctl%ebin_type)
		case(ebin_type_log)
			enx=10**fe%xb(i) 
		end select 
		if(xmin<enx)then
            ier=-99
			call my_gs_integral(xmin,enx, fout, 0d0,-0.5d0,1d-8,1d-7,1,fx,ier)
		else
			fout=0
            return
		end if
		
		fe%fx(i)=(fout +(enx-xmin)**(-0.5)*drhodr_max/dphidr_max)/pi**0.5
        !fe%fx(i)=(fout)/pi**0.5
		!print*, "fout,enx,dp0=",fout,fe%xb(i),enx**(-0.5)*drhodr_max/dphidr_max
        !read(*,*)
		if(fe%fx(i)<0d0.or.ieee_is_nan(fe%fx(i)))then
			print*, "get_init_fe_gs_plummer:i,xmin, en, fe%fx(i)=", i, xmin, enx, fe%fx(i),fout
			fe%fx(i)=0
			!read(*,*)
		end if
	end do
	!fe%fx=fe%fx
    !call fe%print("fe")
    fe%fx=fe%fx/bin_mass
	!read(*,*)
contains 
	real(8) function fx(x)
		use, intrinsic :: ieee_arithmetic
		implicit none
		real(8) x, r,phi,logr

		phi=x 
		call get_rmax_accurate(spp,rmax,log10(phi),logr)
		 
		r=10**logr
        call get_rho_full_range_spp(spp,logr,rho_tmp)
		!call rho%get_value_s(logr,rho_tmp)
        call get_plummer_den(mtot,ra,logr,rho_tmp0)
		call get_beta_full_range(spp,logr, beta_tmp)
		
		dphidr=-(mbhin+beta_tmp)/r**2
		dphidr2=(2*(mbhin+beta_tmp)-4*pi*r**3*rho_tmp)/r**3
		
		gamma_tmp=5*r/(r**2+ra**2)
		gamma_tmp2=7*r/(r**2+ra**2)-1d0/r

		drhodr=-rho_tmp0*gamma_tmp
		drhodr2=rho_tmp0*gamma_tmp*gamma_tmp2

		fx=(drhodr2/dphidr**2-drhodr*dphidr2/dphidr**3)
        !if(fx<0) fx=0
		!print*, "logphi,logr, r,rho,beta,fx=", log10(phi), logr, r,rho_tmp,beta_tmp, fx
		!print*, "fx=", fx
		!stop
		!if(fx<0)then
		!	print*, "x,logx,rho,beta,r=",x,log10(x),rho_tmp,beta_tmp, logr
		!	print*, "drhodr2,dphidr,drhodr,dphidr2=",drhodr2,dphidr,drhodr,dphidr2
		!	print*, "drhodr2*dphidr-drhodr*dphidr2=",drhodr2*dphidr-drhodr*dphidr2
		!	print*, "fx=",(drhodr2*dphidr-drhodr*dphidr2)/dphidr**3
		!	print*, "fx=",fx
		!	read(*,*)
		!end if
	end function
end subroutine

subroutine solving_ini_para_dehnen(mbh_mass,gamma,r0_cl )
    use MPI_comu,only:rid
    use constant
    implicit none
    real(8) mbh_mass, ra, mcl,rh,nh,gamma,r0_cl
    real(8),external::rtbis_yacc
    integer niter, ier
    real(8),parameter::n0=3e4
    real(8) par(50),rib
    print*, "mbh_mass=",mbh_mass,gamma
    rh=3.1*(mbh_mass/4d6)**0.55  !pc
    nh=n0*(mbh_mass/4d6)**(-0.65)  !pc^{-3}
    rib=rd_sun*(mbh_mass)**(1d0/3d0)/206264 ! pc
    ra=10**rtbis_yacc(func,-5d0,3d0,1d-4,par,niter,1000,ier, .false.) !pc
    mcl=nh*4/(3-gamma)*pi*rh**gamma*(rh+ra)**(4-gamma)/ra
    if(rid.eq.0)then
        print*, "dehen model, solving ini para dehnen: mbh_mass,rh, ra,mcl=",mbh_mass,rh, ra/rh,mcl/mbh_mass
    end if
    !print*, "nh=", mcl/2d0/pi*
   ! stop
contains 
    real(8) function func(x,par)
    implicit none
        real(8) x, par(50),mcl_tmp,alpha, masswithinrib
        mcl_tmp=nh*4*pi/(3-gamma)*rh**gamma*(rh+10**x)**(4-gamma)/10**x
        masswithinrib=mcl_tmp*(rib/(rib+10**x))**(3-gamma)
        alpha=((mbh_mass*2d0-masswithinrib)/mcl_tmp)**(1d0/(3-gamma))
        func=alpha*10**x/(1-alpha)-rh
        !print*, "x,mcl_tmp,alpha, func, M(<r0)=",x,mcl_tmp,alpha, func, masswithinrib
    end function
end subroutine