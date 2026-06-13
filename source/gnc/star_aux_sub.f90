subroutine get_sample_para_one_grids(dm,sp,spp)
     use com_main_gw
     implicit none
     type(particle_sample_type)::sp
     type(diffuse_mspec)::dm
     type(star_pot_para)::spp
     real(8) ex, logex, jc, jc_dm
     real(8) rmax, rc,jm,rp_dm,ra_dm
     real(8) jc_dmless,pd_dm,p_EJ_dmless,r_c_iter
     integer ier
     logical condition1, condition2
     if(ctl%chattery.ge.3)then
        print*, "==get_sample_para_grid======================"
                !start=======================================
    end if
     ex=sp%x
     logex=log10(ex)
     sp%en=sp%x*ctl%energy0
     !x=sP%x

     !logx=log10(ex)
    condition1=logex>dm%jc%xb(1).and.logex<dm%jc%xb(dm%jc%nbin)
    condition2=sp%jm<0.99  ! to make sure the result is accurate when jm is very close to one.
    if(condition1.and.condition2)then
        select case(ctl%ebin_type)
        case(ebin_type_log)
            call dm%jc%get_value_l(logex,jc_dm) 
        end select
        if(jc_dm<0)then
            print* ,"jc_dm=",logex, jc_dm
            call dm%jc%print("jc_dm")
            stop
        end if
    else
         
        rc=r_c_iter(spp, ex,ier)
       ! print*, "rc=",rc
        !read(*,*)
        jc_dm=jc_dmless(rc,spp)
    end if

    sp%jc=jc_dm*ctl%v0*r0_cl

     !===========================
    sp%jm=sp%jph/sp%jc
    call set_jm_bound(sp%jm)
     !===========================
    jm=sp%jm
    call get_jm_idx(sp%jm, sample_table_idy,sample_table_rdy,sample_evjum)
     
    sp%jph=jm*sp%jc
    select case(ctl%method_interpolate)
    case(method_int_nearst)
        rp_dm=dm%rp%fxy(sample_table_idx,sample_table_idy)
        sp%rp=rp_dm*r0_cl
        ra_dm=dm%ra%fxy(sample_table_idx,sample_table_idy)

        sp%ra=ra_dm*r0_cl
        pd_dm=dm%pd%fxy(sample_table_idx,sample_table_idy)
        sp%period=pd_dm*r0_cl/ctl%v0
    case(method_int_linear)    
        !print*, "sizeof(rp)=",sizeof(dm%rp%fxy), sample_table_idx,sample_table_idy
        call linear_int_2d_xy(sample_table_idx,sample_table_idy,sample_table_rdx,sample_table_rdy,&
        dm%rp%fxy,dms%df_coe_bins,dms%df_coe_bins,rp_dm)
        !print*, "rp_dm=",rp_dm
        !rp_dm=dm%rp%fxy(sample_table_idx,sample_table_idy)
        !print*, "rp_dm=",rp_dm
        
        sp%rp=rp_dm*r0_cl

        call linear_int_2d_xy(sample_table_idx,sample_table_idy,sample_table_rdx,sample_table_rdy,&
        dm%ra%fxy,dms%df_coe_bins,dms%df_coe_bins,ra_dm)
        !print*, "ra_dm=",ra_dm
        !ra_dm=dm%ra%fxy(sample_table_idx,sample_table_idy)
        sp%ra=ra_dm*r0_cl 

        call linear_int_2d_xy(sample_table_idx,sample_table_idy,sample_table_rdx,sample_table_rdy,&
        common_pd_log%fxy,dms%df_coe_bins,dms%df_coe_bins,pd_dm) 
        
        sp%period=10**pd_dm*r0_cl/ctl%v0 
        !=================================
    end select
    if(ctl%chattery.ge.5)then
        print*, "==end of get_sample_para_grid===================" 
     end if
 end subroutine


subroutine get_barp_xy(logxb,fx,n,rmax,mbhin,spp)
	use com_sts_type
	use constant
	use my_intgl
	use md_star_pot
	implicit none
	!type(s1d_type):: fphi_star
	type(s1d_type)::rmax
	integer n
	real(8) logxb(n), fx(n)
	integer i, idid
	real(8) logx,fout, logrmin, logr, mbhin, rmax_tmp
    type(star_pot_para)::spp

	do i=1, n
       ! print*, "logxb(i)=",logxb(i)
        call get_barp_one(logxb(i),fx(i),rmax,mbhin,spp)
	end do

end subroutine
subroutine get_rho_rmax()
    use com_main_gw
    implicit none
    integer i,j
    do i=1, dms%N
        do j=1,n_tot_comp_sg
    !        !call dms%mb(i)%dsp(j)%p%fden%get_value_s(sample_logrmax,dms%mb(i)%dsp(j)%p%spt_rho_rmax)
            dms%mb(i)%dsp(j)%p%spt_rho_rmax=0
    !        !if(dms%mb(i)%dsp(j)%p%spt_rho_rmax<0) dms%mb(i)%dsp(j)%p%spt_rho_rmax=0
    !        !if(dms%mb(i)%dsp(j)%p%spt_rho_rmax.ne.0)then
    !            !print*, "rid,i,j,rho=",rid,i,j,dms%mb(i)%dsp(j)%p%spt_rho_rmax
    !        !end if
    !        
        end do
    end do
    !read(*,*)
    !type(star_pot_para)::spp
    !call get_rho_full_range_spp(spp,sample_logrmax,spp%spt_rho_rmax)
end subroutine
subroutine get_rho_rmin(spp)
    use com_main_gw
    implicit none
    integer i,j,flag_out
    real(8) rho_rmin(dms%n),rho_rmax
    type(star_pot_para)::spp
    type(s1d_ird_type)::common_gx_ir
    !type(particle_samples_arr_type)::bks
    
    rho_rmin=0
    select case(ctl%fden_ana_est_method)
    case(fden_ana_est_method_1d_iso)

        do i=1, dms%n
            common_gx_ir=dms%mb(i)%all%barge_ir
            !call get_none_zero_s1d_ir(dms%mb(i)%all%barge_ir,common_gx_ir)
        ! print*, "111"
            !do j=1, common_gx_ir%nbin
            !    if(common_gx_ir%fx(j)>0)then
            !        common_gx_ir%fx(j)=log10(dms%mb(i)%all%barge_ir%fx(j))
            !    else
            !        print*, "error! =0!???"
            !        stop
            !!        common_gx_ir%fx(i)=-100
            !    end if
            !end do
            rho_rmax=0
            do j=1,n_tot_comp_sg
                rho_rmax=rho_rmax+dms%mb(i)%dsp(j)%p%spt_rho_rmax
            end do 
            call get_fden_ird_one(sample_logrmin,rho_rmin(i),common_gx_ir ,dms,spp,rho_rmax,flag_out)
            if(flag_out.eq.1)then
                if(rid.eq.0)then
                    print*, "in get_rho_rmin, sample_logrmin=", sample_logrmin
                end if
            end if 
        end do
    case(fden_ana_est_method_2d)
         
        rho_rmin=0
    end select

    spp%spt_rho_rmin=0
    do i=1, dms%n
        spp%spt_rho_rmin=spp%spt_rho_rmin+rho_rmin(i)*dms%mb(i)%mc
    end do 
    spp%has_set_rhomin=.true. 
    if(ctl%fden_ana_est_method.ne.fden_ana_est_method_2d)then
        if(rid.eq.0)then
            print*, "spt_rho_rmin=",spp%spt_rho_rmin
            print*, "rho_rmin(:)=", rho_rmin(:)
        end if
    end if
end subroutine
subroutine get_sample_jlc(ex,mbhin, rt,jc,spp,jlc,ier) 
    use model_basic,only:ctl,r0_cl
    use md_star_pot
    implicit none
    real(8) rt, jc, jlc, jphlc2,ex,mbhin, phi0
    integer ier
    type(star_pot_para)::spp
    call get_phi_star_full_range(spp,log10(rt),phi0)
    jphlc2=2*(mbhin/rt+10**phi0-ex)
    ier=0
    if(jphlc2<=0.and.mbhin.ne.0d0)then
        ier=1
        return 
    end if    
    if(jc.eq.0d0)then
        print*, "get_sample_jlc:error! sp%jc=0"
        stop
    end if
    jlc=rt*(jphlc2)**0.5/jc
    if(ctl%chattery.ge.4)then
        print*, "===get_r_td================================="
        print*, "sample%r_lc,rtd_dm=",rt*r0_cl,rt
        print*, "sample_jlc_dimless=",jlc
     end if
end subroutine
 

subroutine get_sample_rrange_mpi()
    use com_main_gw
    implicit none
    integer i,ierr
    real(8) sendbuffer_r(2,ctl%ntasks), reivbuffer_r(2,ctl%ntasks)
    real(8) ra_dmless, rp_dmless
    
    sample_logrmax=-1d99
    sample_logrmin=1d99
    !sample_logemin=sample_logemax
    do i=1, bksams_arr_norm%n
        ra_dmless=bksams_arr_norm%sp(i)%ra/r0_cl
        if(sample_logrmax<log10(ra_dmless)) sample_logrmax=log10(ra_dmless)
        rp_dmless=bksams_arr_norm%sp(i)%rp/r0_cl
        if(sample_logrmin>log10(rp_dmless)) sample_logrmin=log10(rp_dmless)
    end do

    sendbuffer_r(1,:)=sample_logrmax+0.000001*abs(sample_logrmax)
    sendbuffer_r(2,:)=sample_logrmin-0.000001*abs(sample_logrmin)!+log10(2d0)

	call mpi_alltoall(sendbuffer_r,2, MPI_DOUBLE, reivbuffer_r, 2, &
			MPI_DOUBLE, mpi_comm_world, ierr)    

    sample_logrmax=min(maxval(reivbuffer_r(1,:)),ctl%log10rmax_factor)
    sample_logrmin=max(minval(reivbuffer_r(2,:)),ctl%log10rmin_factor)

   call mpi_barrier(MPI_comm_world, ierr)

   if(rid.eq.0)then
        print*, "log10sample_rmax, log10sample_rmin=",sample_logrmax, sample_logrmin
   end if
end subroutine
subroutine get_sample_erange()
    use com_main_gw
    implicit none
    integer i,ierr
    real(8) sendbuffer_r(2,ctl%ntasks), reivbuffer_r(2,ctl%ntasks)
    
    sample_logemax=-1d99
    sample_logemin=1d99
    !sample_logemin=sample_logemax
    do i=1, bksams_arr_norm%n
        if(sample_logemax<log10(bksams_arr_norm%sp(i)%x)) sample_logemax=log10(bksams_arr_norm%sp(i)%x)
        if(sample_logemin>log10(bksams_arr_norm%sp(i)%x)) sample_logemin=log10(bksams_arr_norm%sp(i)%x)
    end do
    !print*, "bf:rid,sample_logemin,logemax=",rid,sample_logemin,sample_logemax
    sendbuffer_r(1,:)=sample_logemax
    sendbuffer_r(2,:)=sample_logemin

	call mpi_alltoall(sendbuffer_r,2, MPI_DOUBLE, reivbuffer_r, 2, &
			MPI_DOUBLE, mpi_comm_world, ierr)    

    sample_logemax=maxval(reivbuffer_r(1,:))
    sample_logemax=sample_logemax+0.000001*abs(sample_logemax)
    sample_logemin=minval(reivbuffer_r(2,:))
    sample_logemin=sample_logemin-0.000001*abs(sample_logemin)
    !print*, "reivbuffer_r=", reivbuffer_r
   ! print*, "sample_emax,rid=", sample_emax, rid
    sample_emin=10**sample_logemin
    sample_emax=10**sample_logemax
    !print*, "af:rid,sample_logemin,logemax=",rid,sample_logemin,sample_logemax
   call mpi_barrier(MPI_comm_world, ierr)


   if(rid.eq.0)then
        print*, "log10sample_emax, log10sample_emin=",sample_logemax, sample_logemin
   end if

end subroutine
subroutine get_total_sample_number()
    use com_main_gw
    implicit none
    integer ns
    ns=bksams_arr_norm%n
    call collection_int(ns)
    ctl%n_tot_samples=ns
end subroutine
subroutine get_sample_rrange(logrmin,logrmax)
    use com_main_gw
    implicit none
    real(8) logrmin, logrmax

    select case(ctl%fden_ana_est_method)
    case(fden_ana_est_method_1d_iso)
        call get_rmax_accurate(spp_new ,dms%fr_phi,sample_nxgx_logemin,logrmax)
        call get_rmax_accurate(spp_new ,dms%fr_phi,sample_nxgx_logemax,logrmin)
        dms%logrmin=logrmin 
        dms%logrmax=logrmax 
        sample_logrmin=dms%logrmin; sample_logrmax=dms%logrmax
    case(fden_ana_est_method_2d)
        call get_sample_rrange_mpi()
        dms%logrmin=sample_logrmin
        dms%logrmax=sample_logrmax
    end select
 
end subroutine
subroutine get_barp_one(logenx,fx,rmax,mbhin,spp)
	use com_sts_type
	use constant
	use my_intgl
	use md_star_pot
    use model_basic,only:ctl,dms!,log10rmin_intrinisc
    use md_coeff
    use MPI_comu,only:rid
	implicit none
	!type(s1d_type):: fphi_star
	type(s1d_type)::rmax
    type(star_pot_para)::spp
	real(8) logenx, fx, enx
	integer i, idid
	real(8) fout, logrmin, logr, mbhin, rmax_tmp
    
	logrmin=ctl%log10rmin_factor  !log10rmin_intrinisc;
    fout=0
    enx=10**logenx 

    call get_rmax_accurate(spp,rmax,logenx,rmax_tmp) 
    call my_integral_acc(logrmin,rmax_tmp,fout,1d-24,1d-14, FCN, idid) 
    fx=2**1.5d0*abs(fout)*log(10d0)
contains 
	subroutine FCN(N,X,Y,F,RPAR,IPAR)
		implicit none
		integer n, ipar(100)
		real(8) x, y(n), f(n), rpar(100),ysp, phi_tmp, radius
		!call fphi_tot%get_value_s(x, phi_tmp)
		call get_phi_star_full_range(spp,x,phi_tmp)
        select case(ctl%ebin_type)
        case(ebin_type_log)
            radius=10**x 
            F(1)=(radius**2.5d0)*(max((10**phi_tmp-enx)*radius+mbhin,0d0))**0.5 
        end select 
	end subroutine
end subroutine 
subroutine dms_so_get_fslope(so, source)
    use md_stellar_object
    use model_basic,only:ctl
    use MPI_comu
    implicit none
    type(dms_stellar_object)::so
    type(s1d_type)::fden_temp, fden
    integer i, ilast, source,ier
    real(8) xmin
 
        select case(source)
        case(1)
            fden=so%fden
        case(2)
            fden=so%fden_simu
        end select
        
        ilast=0
        do i=1,fden%nbin
            if(fden%fx(i).gt.0)then
                ilast=ilast+1
            endif
        end do
        
        if(ilast>2)then
            call fden_temp%init(fden%xmin,fden%xmax,ilast,sts_type_dstr)
            ilast=0
            do i=1, fden%nbin
                if(fden%fx(i)>0)then
                    ilast=ilast+1
                    fden_temp%xb(ilast)=fden%xb(i)
                    fden_temp%fx(ilast)=log10(fden%fx(i))
                end if
            end do 
            call fden_temp%get_dfdx_spline(so%fslope) 
        else
            so%fslope%fx=-999d0
        end if 
end subroutine

subroutine get_slope0(dm)
    use com_main_gw
    implicit none
    type(diffuse_mspec)::dm
    integer i,j,k,ier
     
    call get_fslope(dm,source_ana) 
end subroutine
  
subroutine set_gx_nx_ranges_ir(dm)
    use md_dms
    use model_basic,only:ctl
    use MPI_comu,only:rid 
    implicit none
    type(diffuse_mspec)::dm
    integer i
 
    do i=1, dm%n
       call set_gx_nx_ranges_mb(dm%mb(i))
    end do
    call set_gx_nx_ranges_mb(dm%all)
    
end subroutine
subroutine set_nx_1d_ranges(dm)
    use md_dms
    implicit none
    type(diffuse_mspec)::dm
    integer i
    do i=1, dm%n
        call set_nx1d_ranges_mb(dm%mb(i))
    end do
    call set_nx1d_ranges_mb(dm%all)
end subroutine

subroutine set_nx_1d_ranges_ir(dm)
    use md_dms
    implicit none
    type(diffuse_mspec)::dm
    integer i
    do i=1, dm%n
        call set_nx1d_ird_ranges_mb(dm%mb(i),dm%barp_ir%xb,dm%barp_ir%xsteps,dm%barp_ir%nbin,&
        dm%barp_ir%xmin,dm%barp_ir%xmax)
    end do
    call set_nx1d_ird_ranges_mb(dm%all,dm%barp_ir%xb,dm%barp_ir%xsteps,dm%barp_ir%nbin,&
    dm%barp_ir%xmin,dm%barp_ir%xmax)
end subroutine
subroutine set_gx_nx2d_ranges_ir(dm)
    use md_dms
    use model_basic,only:ctl
    use MPI_comu,only:rid 
    implicit none
    type(diffuse_mspec)::dm
    integer i
 
    do i=1, dm%n
       call set_gx_nx2d_ranges_mb(dm%mb(i))
    end do
    call set_gx_nx2d_ranges_mb(dm%all)
    
end subroutine

subroutine update_re_tables()
    use com_main_gw
    implicit none
    real(8) t1, t2 
    call get_diffuse_mspec_rtables(dms,spp_new)
     
    call get_ini_ebounds()
    call get_diffuse_mspec_ebound() 
    call init_diffuse_mspec_etables(dms)
    
    call get_frphi(spp_new%fphi_star,dms%fr_phi)
    call get_xrc(dms%frc_x,spp_new)
     
    if(ctl%chattery.ge.2.and.rid.eq.0)then
        print*, "end of update_re_tables"
    end if
     
end subroutine
 
subroutine set_barp(dm,spp)
    use com_main_gw
    implicit none
    type(diffuse_mspec)::dm
    type(star_pot_para)::spp
    real(8) t1, t2
    !if(rid.eq.0)then
    !    call cpu_time(t1)
    !end if
    
    if(dm%e_iregular)then
        !print*, "set_barp:sample_logemin,logemax=", sample_logemin,sample_logemax
        sample_nxgx_logemin=max(sample_logemin,log10(dms%emin))
        sample_nxgx_logemax=min(log10(dms%emax),sample_logemax)


        !if(rid.eq.0.and.ctl%chattery.ge.1)then
        !    print*, "setting barp_ir:min,max,dms%logemax,sample_emax=",&
        !        sample_logemin,min(log10(dms%emax),10**sample_logemax),log10(dms%emax),10**sample_logemax
        !end if
        !if(rid.eq.0)then
        !    print*, "semin,semax=",sample_nxgx_logemin,sample_nxgx_logemax
        !end if
        call set_barp_ir(dm,dm%barp_ir,spp)

        !stop
    else
        call dm%barp%init(10**sample_logemin,10**sample_logemax,dm%dstr_bins_e,sts_type_dstr)
        call dm%barp%set_range()
        call set_barp_re(dm%barp,spp)
        if(rid.eq.0.and.ctl%chattery.ge.2)then
            call dm%barp%print("barp")
        end if
    end if
    !if(rid.eq.0)then
       ! call cpu_time(t2)
        !print*, "set_barp used time=", t2-t1, " s"
        !call dm%barp_ir%print("barp")
        !read(*,*)
    !end if
end subroutine
subroutine set_barp_re(barp,spp)
    use com_main_gw
    implicit none
    integer i
    type(s1d_type)::barp
    type(star_pot_para)::spp
    !common_barp=dm%mb(1)%star%barge
    select case(ctl%ebin_type)
    case(ebin_type_log) 
        call get_barp_xy(barp%xb,barp%fx,barp%nbin, &
        dms%fr_phi,spp%mbh_dmless,spp) 
    end select

end subroutine
 
subroutine set_barp_step_size(dm, barp_ir,spp)
    use com_main_gw
    implicit none
    integer i
    type(s1d_ird_type)::barp_ir
    type(s1d_type)::barp
    type(diffuse_mspec)::dm
    type(star_pot_para)::spp
    integer nt
 
    call dms%barp_ir%init(sample_nxgx_logemin,sample_nxgx_logemax,&
    dms%dstr_bins_e,sts_type_dstr)

    call set_nx_ranges_xb_ir(barp_ir%xb,barp_ir%xsteps, barp_ir%xmin, &
    barp_ir%xmax,barp_ir%nbin,barp_ir%bin_type,dm,spp,ctl%barge_grid_type)

end subroutine
subroutine set_barp_ir(dm, barp_ir,spp)
    use com_main_gw
    implicit none
    integer i
    type(s1d_ird_type)::barp_ir
    type(s1d_type)::barp
    type(diffuse_mspec)::dm
    type(star_pot_para)::spp
    integer nt
     
    call set_barp_step_size(dm, barp_ir,spp)
     
    call get_barp_ir(dm,barp_ir,spp)
end subroutine
subroutine get_barp_ir(dm,barp_ir,spp)
    use com_main_gw
    implicit none
    integer i
    type(s1d_ird_type)::barp_ir
    type(s1d_type)::barp
    type(diffuse_mspec)::dm
    type(star_pot_para)::spp

    select case(ctl%ebin_type)
    case(ebin_type_log)
        call get_barp_xy(barp_ir%xb,barp_ir%fx,barp_ir%nbin, &
        dm%fr_phi,spp%mbh_dmless,spp) 
    end select
end subroutine
 
subroutine set_gx_nx_ranges_mb(mb)
    use com_main_gw
    implicit none
    integer i
    type(mass_bins)::mb 
    do i=1, n_tot_comp 
        associate(b=>mb%dsp(i)%p%nx_ir, d=>mb%dsp(i)%p%barge_ir) 

            call b%init(sample_logemin,sample_logemax,dms%dstr_bins_e,use_weight=.true.)
 
            call d%init(sample_logemin,sample_logemax,dms%dstr_bins_e,sts_type_dstr)
             b%xb=dms%barp_ir%xb
             b%xsteps=dms%barp_ir%xsteps 
             d%xb=b%xb
             d%xsteps=b%xsteps

        end associate
    end do
    associate(b=>mb%all%nx_ir,d=>mb%all%barge_ir)

        call b%init(sample_logemin,sample_logemax,dms%dstr_bins_e,use_weight=.true.)
 
        call d%init(sample_logemin,sample_logemax,dms%dstr_bins_e,sts_type_dstr)
        b%xb=dms%barp_ir%xb
        b%xsteps=dms%barp_ir%xsteps 
        d%xb=b%xb
        d%xsteps=b%xsteps
    end associate
end subroutine

subroutine set_nx1d_ranges_mb(mb)
    use com_main_gw
    implicit none
    integer i
    type(mass_bins)::mb
    real(8) ysteps 
    do i=1, n_tot_comp
        associate(nx=>mb%dsp(i)%p%nx)
            call nx%init(sample_nxgx_logemin,sample_nxgx_logemax,30 ,use_weight=.true.)
            call nx%set_range()
        end associate
    end do
    associate(nx=>mb%all%nx)
        call nx%init(sample_nxgx_logemin,sample_nxgx_logemax,30 ,use_weight=.true.)
        call nx%set_range()        
    end associate

end subroutine

subroutine set_nx1d_ird_ranges_mb(mb, xb,xsteps,n,xmin,xmax)
    use com_main_gw
    implicit none
    integer i
    type(mass_bins)::mb
    real(8) ysteps
    integer n
    real(8) xb(n),xsteps(n),xmin,xmax

    
    do i=1, n_tot_comp
        associate(nx=>mb%dsp(i)%p%nx_ir)
            call nx%init(xmin,xmax,n ,use_weight=.true.)
             nx%xb=xb
             nx%xsteps=xsteps !*ctl%dstr_bins_e_amplifier
        end associate
    end do
    associate(nx=>mb%all%nx_ir)
        
        call nx%init(xmin,xmax,n ,use_weight=.true.)

        nx%xb=xb
        nx%xsteps=xsteps !*ctl%dstr_bins_e_amplifier
        
    end associate

end subroutine


subroutine set_gx_nx2d_ranges_mb(mb)
    use com_main_gw
    implicit none
    integer i
    type(mass_bins)::mb
    real(8) ysteps
    integer n
    !type(s1d_type)::barp

    !barp=dms%jc
    !call set_barp_re(barp)
    !call barp%print("barp")
    n=dms%barp_ir%nbin 
    do i=1, n_tot_comp
        !call set_gx_ranges_obj(mb%dsp(i)%p%barge_ir)
        associate(nx2=>mb%dsp(i)%p%nxj_ir, &
                gxj=>mb%dsp(i)%p%gxj_ir, barge=>mb%dsp(i)%p%barge_ir)
            
            call nx2%init(n,ctl%dstr_bins_j,sample_logemin,sample_logemax,&
                log10(dms%jmin),log10(dms%jmax),&
                use_weight=.true.)

            call gxj%init(n,ctl%dstr_bins_j,sample_logemin,sample_logemax,&
                log10(dms%jmin),log10(dms%jmax),sts_type_dstr)

            call barge%init(sample_logemin,sample_logemax,n,sts_type_dstr)
            barge%xb=dms%barp_ir%xb
            barge%xsteps=dms%barp_ir%xsteps!*2!*ctl%dstr_bins_e_amplifier

             nx2%xcenter=barge%xb
             nx2%xsteps=barge%xsteps 
             call set_range(nx2%ycenter,nx2%ny,nx2%ymin,nx2%ymax,sts_type_dstr)
             ysteps=nx2%ycenter(2)-nx2%ycenter(1)
             nx2%ysteps=ysteps
             gxj%xcenter=nx2%xcenter
             gxj%ycenter=nx2%ycenter
             gxj%ystep=nx2%ysteps(1)
        end associate
    end do
    associate(nx2=>mb%all%nxj_ir, &
            gxj=>mb%all%gxj_ir, barge=>mb%all%barge_ir)
        
        call nx2%init(n,ctl%dstr_bins_j,sample_logemin,sample_logemax,&
            log10(dms%jmin),log10(dms%jmax),&
            use_weight=.true.)

        call gxj%init(n,ctl%dstr_bins_j,sample_logemin,sample_logemax,&
            log10(dms%jmin),log10(dms%jmax),sts_type_dstr) 
        call barge%init(sample_logemin,sample_logemax,n,sts_type_dstr)
        barge%xb=dms%barp_ir%xb
        barge%xsteps=dms%barp_ir%xsteps!*2

        nx2%xcenter=barge%xb
        nx2%xsteps=barge%xsteps
        call set_range(nx2%ycenter,nx2%ny,nx2%ymin,nx2%ymax,sts_type_dstr)
        ysteps=nx2%ycenter(2)-nx2%ycenter(1)
        nx2%ysteps=ysteps
        gxj%xcenter=nx2%xcenter
        gxj%ycenter=nx2%ycenter
        gxj%ystep=nx2%ysteps(1)
        if(ctl%fden_ana_est_method.eq.fden_ana_est_method_2d)then
            call dms%jc_sample_erange%init(sample_logemin,sample_logemax,dms%dstr_bins_e,sts_type_dstr)
            dms%jc_sample_erange%xb=barge%xb
        end if
    end associate

end subroutine


subroutine set_etable_ir_bins()
    use com_main_gw
    implicit none
    integer ierr, n 
    call dms%dlxb_ir%init(log10(dms%emin),log10(dms%emax),&
        dms%df_coe_bins,sts_type_dstr)
    n=dms%dlxb_ir%nbin
    call set_nx_ranges_xb_ir(dms%dlxb_ir%xb(1:n),dms%dlxb_ir%xsteps(1:n), dms%dlxb_ir%xmin, &
    dms%dlxb_ir%xmax,n,dms%dlxb_ir%bin_type,dms,spp_new,ctl%barge_grid_type)
 
    dms%rc%xb=dms%dlxb_ir%xb
 
    dms%jc%xb=dms%dlxb_ir%xb 
    dms%pd%xcenter=dms%dlxb_ir%xb
    dms%rp%xcenter=dms%dlxb_ir%xb
    dms%ra%xcenter=dms%dlxb_ir%xb        
 
    dms%rc%xmax=dms%dlxb_ir%xmax
    dms%jc%xmax=dms%dlxb_ir%xmax 
    dms%pd%xmax=dms%dlxb_ir%xmax
    dms%rp%xmax=dms%dlxb_ir%xmax
    dms%ra%xmax=dms%dlxb_ir%xmax  
 
    dms%rc%xmin=dms%dlxb_ir%xmin
    dms%jc%xmin=dms%dlxb_ir%xmin
    dms%pd%xmin=dms%dlxb_ir%xmin
    dms%rp%xmin=dms%dlxb_ir%xmin
    dms%ra%xmin=dms%dlxb_ir%xmin      
end subroutine
subroutine set_gx_ranges_xb(xb,xsteps, xmin,xmax,n,bin_type,jc, pd,&
    gx_func_max_step,gx_func_min_step)
    !use com_main_gw
    use com_sts_type
    implicit none
    integer n, i, bin_type
    real(8) xb(n),xsteps(n),xmin,xmax,ex
    real(8) xstep_max, xstep_min,gx_func_max_step,gx_func_min_step
    real(8) jc_xy,pd_xy,xsum
    type(s1d_type)::jc
    type(s2d_type)::pd
    !print*, "??"
    !stop
    if(n<1) return
    call set_range(xb,n,xmin,xmax,bin_type)
    xstep_max=(xmax-xmin)/real(n)*gx_func_max_step
    xstep_min=(xmax-xmin)/real(n)*gx_func_min_step

    do i=1, n
        ex=10**xb(i)
        call jc%get_value_s(xb(i),jc_xy)
        call pd%get_value_d(xb(i),-0.5d0,pd_xy)
        xsteps(i)=(ex*jc_xy**2*pd_xy)**0.5
    end do
    xsum=sum(xsteps)
    do i=1, n
        xsteps(i)=max(min(xstep_max,xsteps(i)/xsum*(xmax-xmin)),xstep_min)
    end do
    xsum=sum(xsteps)
    do i=1, n
        xsteps(i)=xsteps(i)/xsum*(xmax-xmin)
    end do
    
    xb(1)=xmin+xsteps(1)/2d0
    do i=2, n
        xb(i)=xb(i-1)+(xsteps(i)+xsteps(i-1))/2d0
    end do
    
end subroutine 

subroutine set_nx_ranges_xb_ir(xb,xsteps, xmin,xmax,n,bin_type,dm,spp,type)
    use com_main_gw
    use com_sts_type
    implicit none
    integer n, i, bin_type,type
    real(8) xb(n),xsteps(n),xmin,xmax,ex, rmin,rmax
    real(8) xstep_max, xstep_min!,gx_func_max_step,gx_func_min_step
    real(8) phi_out,phi_tot
    type(diffuse_mspec)::dm
    type(star_pot_para)::spp
    integer ierr
    !real(8)

    select case(type)
    case(barge_grid_type_iregular_phi)
        call set_nx_ranges_xb_ir_phi(xb(1:n),xsteps(1:n), xmin,xmax,n,bin_type,dm,spp)
        !call set_nx_ranges_xb_ir_phi_n(xb(1:n),xsteps(1:n), xmin,xmax,n,bin_type,spp)
    case(barge_grid_type_iregular_jc)

        call set_nx_ranges_xb_ir_jc(xb(1:n),xsteps(1:n), xmin,xmax,n,bin_type,dm,spp)
    case default
        print*, "error, set_nx_ranges_xb_ir:define barge_grid type", type
        stop
    end select
 
end subroutine  
 
subroutine set_nx_ranges_xb_ir_phi(xb,xsteps, xmin,xmax,n,bin_type,dm,spp)
    use md_star_pot
    use com_sts_type
    use md_dms
    use MPI_comu,only:rid
    use model_basic,only:gx_func_max_step,gx_func_min_step
    implicit none
    integer n, i, bin_type,nt
    real(8) xb(n),xsteps(n),xmin,xmax,ex, rmin,rmax
    real(8) xstep_max, xstep_min!,gx_func_max_step,gx_func_min_step
    real(8) phi_out,phi_tot
    type(diffuse_mspec)::dm
    type(star_pot_para)::spp
    real(8) rb(n+1), xg(n+1)
    !real(8) xstep_min,xstep_max
    if(n<1) return
    !print*, "set_nx_ranges_xb_ir: xmin,xmax=",xmin,xmax
    !print*, "dm%fr_phi%xmin,xmax=",dm%fr_phi%xmin,dm%fr_phi%xmax
    !read(*,*)
    !print*, "bg",rid
    
    call get_rmax_accurate(spp,dm%fr_phi,xmin,rmax)
    call get_rmax_accurate(spp,dm%fr_phi,xmax,rmin)
    !print*, "xmin,xmax,rmin,rmax=",xmin,xmax,rmin,rmax,rid
    call set_range(rb,n+1,rmin,rmax,sts_type_grid)
    
    do i=n+1,1,-1
        call get_phi_star_full_range(spp,rb(i),phi_out)
        if(spp%mbh_dmless.ne.0)then
            xg(n-i+2)=log10(10**phi_out+spp%mbh_dmless/10**rb(i))
        else
            xg(n-i+2)=phi_out
        end if
    ! print*, "r, xb=",rb(i),xg(n-i+1)
    end do
    do i=1, n
        
        xb(i)=(xg(i+1)+xg(i))/2d0
        if(xg(i+1).ne.xg(i))then
            xsteps(i)=(xg(i+1)-xg(i))
        else
            xsteps(i)=1d-20
        end if
    end do
    !print*, "xb=",xb
    !print*, "xsteps=",xsteps

    xstep_min=(xmax-xmin)/dble(n)*gx_func_min_step
    xstep_max=(xmax-xmin)/dble(n)*gx_func_max_step

   ! print*, "xstep_min,max=",xstep_min,xstep_max
    if(xstep_min<0.or.xstep_max<0) return
    do i=1, n
        if(xsteps(i)<xstep_min)then
            xsteps(i)=xstep_min
        end if
        if(xsteps(i)>xstep_max)then
            xsteps(i)=xstep_max
        end if
    end do
   ! print*, "xsteps=",xsteps
    xsteps=xsteps/sum(xsteps)*(xmax-xmin)
    !print*, "xsteps=",xsteps
    xb(1)=xmin+xsteps(1)/2d0
    do i=2, n-1
        xb(i)=xb(i-1)+(xsteps(i)+xsteps(i-1))/2d0
    end do
    xb(n)=xmax-xsteps(n)/2d0
    !print*, "xb=",xb
    !stop
    !print*, "ed=",rid
end subroutine

subroutine set_nx_ranges_xb_ir_jc( xb,xsteps, xmin,xmax,n,bin_type,dm,spp)
    use md_star_pot
    use com_sts_type
    use md_dms
    use MPI_comu,only:rid
    use model_basic,only:gx_func_max_step,gx_func_min_step
    implicit none
    integer n, i, bin_type,nt
    real(8) xb(n),xsteps(n),xmin,xmax,ex,jcmin,jcmax
    real(8) xstep_max, xstep_min!,gx_func_max_step,gx_func_min_step
    real(8) phi_out,phi_tot, rcmin,rcmax
    type(diffuse_mspec)::dm
    type(s1d_type)::jc_inv,s1d_jc,rc
    type(star_pot_para)::spp
    real(8) jc(n+1), xg(n+1)
    integer ier 
    if(n<1) return 
    call get_xrc(dm%frc_x,spp)
    call get_jc_minmax(xmin,xmax,rcmin,rcmax,jcmin,jcmax,spp) 
    call set_range(jc,n+1,log10(jcmin),log10(jcmax),sts_type_grid)

    xg(n+1)=xmax
    xg(1)=xmin
    do i=2,n
        !call get_phi_star_full_range(spp,rb(i),phi_out)
        !call jc_inv%get_value_s(jc(i),xg(n-i+2))
        call get_x_given_jc(rcmin,rcmax,xg(n-i+2),10**jc(i),spp,dm)
        !print*, "xg,jc=",xg(n-i+2),10**jc(i)
    end do
    
    do i=1, n
        
        xb(i)=(xg(i+1)+xg(i))/2d0
        if(xg(i+1).ne.xg(i))then
            xsteps(i)=(xg(i+1)-xg(i))
        else
            xsteps(i)=1d-20
        end if
    end do
    !print*, "xb=",xb
    !print*, "xsteps=",xsteps
   ! read(*,*)

    xstep_min=gx_func_min_step*(xmax-xmin)
    xstep_max=gx_func_max_step*(xmax-xmin)

   ! print*, "xstep_min,max=",xstep_min,xstep_max
    if(xstep_min<0.or.xstep_max<0) return
    do i=1, n
        if(xsteps(i)<xstep_min)then
            xsteps(i)=xstep_min
        end if
        if(xsteps(i)>xstep_max)then
            xsteps(i)=xstep_max
        end if
    end do
   ! print*, "xsteps=",xsteps
    xsteps=xsteps/sum(xsteps)*(xmax-xmin)
    !print*, "xsteps=",xsteps
    xb(1)=xmin+xsteps(1)/2d0
    do i=2, n-1
        xb(i)=xb(i-1)+(xsteps(i)+xsteps(i-1))/2d0
    end do
    xb(n)=xmax-xsteps(n)/2d0
    !print*, "xb=",xb
    !stop
    !print*, "ed=",rid
end subroutine



subroutine get_none_zero_s1d_ir_one_zero(barge, barge_nonzero)
    use com_sts_type
    implicit none
    type(s1d_ird_type)::barge, barge_nonzero
    !real(8) xb(barge%nbin), fx(barge%nbin)
    integer i, nnz, nf, ib, ie, ibn, ien
    logical first_zero_ib, first_zero_ie
        ! remove zeros
    ib=1; ie=barge%nbin
    first_zero_ib=.true.
    first_zero_ie=.true.
    nnz=0
    do while (ib<ie)
        if(first_zero_ib)then
            if(barge%fx(ib).ne.0)then
                first_zero_ib=.false.
            end if
            nnz=nnz+1
            !print*, "first_zero_ib,barge%fx(ib),nnz=",first_zero_ib,barge%fx(ib),nnz
        else
            
            if(barge%fx(ib).ne.0)then
                nnz=nnz+1
            end if
            !print*, "first_zero_ib,barge%fx(ib),nnz=",first_zero_ib,barge%fx(ib),nnz
        endif
        
        if(barge%fx(ie).ne.0)then
            first_zero_ie=.false.
            nnz=nnz+1
        end if
        ib=ib+1
        ie=ie-1
    end do
    !print*, "nnz,ib,ie=",nnz,ib,ie
    !end if
    if(nnz.ne.barge%nbin) then
        
        call barge_nonzero%init(barge%xmin,barge%xmax,nnz,sts_type_dstr)

        ib=1; ie=barge%nbin
        first_zero_ib=.true.
        first_zero_ie=.true.
        ibn=0; ien=nnz+1
        do while (ib<ie)
            if(first_zero_ib)then
                
                if(barge%fx(ib).ne.0)then
                    first_zero_ib=.false.
                end if
                ibn=ibn+1
                barge_nonzero%xb(ibn)=barge%xb(ib)
                barge_nonzero%fx(ibn)=barge%fx(ib)
                !print*, "first_zero_ib,barge%fx(ib),nnz=",first_zero_ib,barge%fx(ib),ibn
            else
                
                if(barge%fx(ib).ne.0)then
                    ibn=ibn+1
                    barge_nonzero%xb(ibn)=barge%xb(ib)
                    barge_nonzero%fx(ibn)=barge%fx(ib)
                end if
                !print*, "first_zero_ib,barge%fx(ib),nnz=",first_zero_ib,barge%fx(ib),ibn
            endif
            
            if(barge%fx(ie).ne.0)then
                ien=ien-1
                barge_nonzero%xb(ien)=barge%xb(ie)
                barge_nonzero%fx(ien)=barge%fx(ie)    
            end if            
            ib=ib+1
            ie=ie-1
        end do    
        do i=barge_nonzero%nbin, 1, -1
            if(barge_nonzero%fx(i).eq.0)then
                barge_nonzero%xmax=barge_nonzero%xb(i)
            else
                return
            end if
        end do
    end if
    call barge_nonzero%prepare_spline()
end subroutine
subroutine get_none_zero_s1d_ir_two_zero(barge, barge_nonzero)
    use com_sts_type
    implicit none
    type(s1d_ird_type)::barge, barge_nonzero
    !real(8) xb(barge%nbin), fx(barge%nbin)
    integer i, nnz, nf, ib, ie, ibn, ien
    logical first_zero_ib, first_zero_ie
        ! remove zeros
    ib=1; ie=barge%nbin
    first_zero_ib=.true.
    first_zero_ie=.true.
    nnz=0
    do while (ib<ie)
        if(first_zero_ib)then
            if(barge%fx(ib).ne.0)then
                first_zero_ib=.false.
            end if
            nnz=nnz+1
            !print*, "first_zero_ib,barge%fx(ib),nnz=",first_zero_ib,barge%fx(ib),nnz
        else
            
            if(barge%fx(ib).ne.0)then
                nnz=nnz+1
            end if
            !print*, "first_zero_ib,barge%fx(ib),nnz=",first_zero_ib,barge%fx(ib),nnz
        endif
        if(first_zero_ie)then
            if(barge%fx(ie).ne.0)then
                first_zero_ie=.false.
            end if
            nnz=nnz+1
            !print*, "first_zero_ie,barge%fx(ie),nnz=",first_zero_ie,barge%fx(ie),nnz
        else
            if(barge%fx(ie).ne.0)then
                nnz=nnz+1
            end if
            !print*, "first_zero_ie,barge%fx(ie),nnz=",first_zero_ie,barge%fx(ie),nnz
        endif        
        ib=ib+1
        ie=ie-1
    end do
    !print*, "nnz,ib,ie=",nnz,ib,ie
    !end if
    if(nnz.ne.barge%nbin) then
        call barge_nonzero%init(barge%xmin,barge%xmax,nnz,sts_type_dstr)

        ib=1; ie=barge%nbin
        first_zero_ib=.true.
        first_zero_ie=.true.
        ibn=0; ien=nnz+1
        do while (ib<ie)
            if(first_zero_ib)then
                
                if(barge%fx(ib).ne.0)then
                    first_zero_ib=.false.
                end if
                ibn=ibn+1
                barge_nonzero%xb(ibn)=barge%xb(ib)
                barge_nonzero%fx(ibn)=barge%fx(ib)
                !print*, "first_zero_ib,barge%fx(ib),nnz=",first_zero_ib,barge%fx(ib),ibn
            else
                
                if(barge%fx(ib).ne.0)then
                    ibn=ibn+1
                    barge_nonzero%xb(ibn)=barge%xb(ib)
                    barge_nonzero%fx(ibn)=barge%fx(ib)
                end if
                !print*, "first_zero_ib,barge%fx(ib),nnz=",first_zero_ib,barge%fx(ib),ibn
            endif
            if(first_zero_ie)then
                if(barge%fx(ie).ne.0)then
                    first_zero_ie=.false.
                end if
                ien=ien-1
                barge_nonzero%xb(ien)=barge%xb(ie)
                barge_nonzero%fx(ien)=barge%fx(ie)
                !print*, "first_zero_ie,barge%fx(ie),nnz=",first_zero_ie,barge%fx(ie),ien
            else
                if(barge%fx(ie).ne.0)then
                    ien=ien-1
                    barge_nonzero%xb(ien)=barge%xb(ie)
                    barge_nonzero%fx(ien)=barge%fx(ie)
                end if
                !print*, "first_zero_ie,barge%fx(ie),nnz=",first_zero_ie,barge%fx(ie),ien
            endif        
            ib=ib+1
            ie=ie-1
        end do    
    end if
    call barge_nonzero%prepare_spline()
end subroutine

subroutine get_none_zero_two_side_s1d_ir(barge, barge_nonzero)
    use com_sts_type
    implicit none
    type(s1d_ird_type)::barge, barge_nonzero
    !real(8) xb(barge%nbin), fx(barge%nbin)
    integer i, nnz, nf, ib, ie, ibn, ien
    real(8) xmin, xmax
    logical::first_i,first_e
        ! remove zeros
    
    nnz=0
    xmin=barge%xmin; xmax=barge%xmax
    first_i=.true.
    first_e=.true.
    do i=1, barge%nbin
        if(barge%fx(i).ne.0)then
            if(i>=2.and.first_i)then
                if(barge%fx(i-1).eq.0)then
                    xmin=barge%xb(i)-barge%xsteps(i)/2d0
                end if
                first_i=.false.
            end if
            nnz=nnz+1
        end if
    end do

    do i=barge%nbin,1,-1
        if(barge%fx(i).ne.0)then
            if(i<=barge%nbin-1.and.first_e)then
                if(barge%fx(i+1).eq.0)then
                    xmax=barge%xb(i)+barge%xsteps(i)/2d0
                end if
                first_e=.false.
            end if
        end if
    end do

    !print*, "nnz,ib,ie=",nnz,ib,ie
    !end if
    if(nnz.ne.barge%nbin) then  
        call barge_nonzero%init(xmin,xmax,nnz,sts_type_dstr)
        nnz=0
        do i=1, barge%nbin
            if(barge%fx(i).ne.0)then
                nnz=nnz+1
                barge_nonzero%xb(nnz)=barge%xb(i)
                barge_nonzero%fx(nnz)=barge%fx(i)
                barge_nonzero%xsteps(nnz)=barge%xsteps(i)
            end if
        end do
    end if
    call barge_nonzero%prepare_spline()
end subroutine

subroutine get_none_zero_s1d_ir(barge, barge_nonzero)
    use com_sts_type
    implicit none
    type(s1d_ird_type)::barge, barge_nonzero
    !real(8) xb(barge%nbin), fx(barge%nbin)
    integer i, nnz, nf, ib, ie, ibn, ien
    real(8) xmin, xmax
    logical::first_i,first_e
        ! remove zeros
    
    nnz=0
    
    do i=1, barge%nbin
        if(barge%fx(i).ne.0)then
            nnz=nnz+1
        end if
    end do

    if(nnz.ne.barge%nbin) then  
        call barge_nonzero%init(barge%xmin,barge%xmax,nnz,sts_type_dstr)
        nnz=0
        do i=1, barge%nbin
            if(barge%fx(i).ne.0)then
                nnz=nnz+1
                barge_nonzero%xb(nnz)=barge%xb(i)
                barge_nonzero%fx(nnz)=barge%fx(i)
                barge_nonzero%xsteps(nnz)=barge%xsteps(i)
            end if
        end do
    end if
    barge_nonzero%xmin=barge_nonzero%xb(1)-barge_nonzero%xsteps(1)/2d0
    barge_nonzero%xmax=barge_nonzero%xb(nnz)+barge_nonzero%xsteps(nnz)/2d0
    call barge_nonzero%prepare_spline()
end subroutine

subroutine get_none_zero_s1d(barge, barge_nonzero)
    use com_sts_type
    implicit none
    type(s1d_type)::barge, barge_nonzero
    !real(8) xb(barge%nbin), fx(barge%nbin)
    integer i, nnz, nf, ib, ie, ibn, ien
    real(8) xmin, xmax
    logical::first_i,first_e
        ! remove zeros
    
    nnz=0
    xmin=barge%xmin; xmax=barge%xmax
    first_i=.true.
    first_e=.true.
    do i=1, barge%nbin
        if(barge%fx(i).ne.0)then
            if(i>=2.and.first_i)then
                if(barge%fx(i-1).eq.0)then
                    xmin=barge%xb(i)-barge%xstep/2d0
                    first_i=.false.
                end if
            end if
            if(i<=barge%nbin-1.and.first_e)then
                if(barge%fx(i+1).eq.0)then
                    xmax=barge%xb(i)+barge%xstep/2d0
                    first_e=.false.
                end if
            end if
            nnz=nnz+1
        end if
    end do
    !print*, "nnz,ib,ie=",nnz,ib,ie
    !end if
    if(nnz.eq.barge%nbin) return
    
    call barge_nonzero%init(xmin,xmax,nnz,sts_type_dstr)
    nnz=0
    do i=1, barge%nbin
        if(barge%fx(i).ne.0)then
            nnz=nnz+1
            barge_nonzero%xb(nnz)=barge%xb(i)
            barge_nonzero%fx(nnz)=barge%fx(i)
            !barge_nonzero%xsteps(nnz)=barge%xsteps(i)
        end if
    end do
end subroutine

subroutine get_r_for_sample_stpt(bkarrs, spp,r)
    use com_main_gw
    use md_star_pot
    implicit none
    type(star_pot_para)::spp
    type(particle_samples_arr_type)::bkarrs
    real(8) r(bkarrs%n)
    integer i
    do i=1, bkarrs%n
        call get_r_stpt(spp, bkarrs%sp(i)%x,bkarrs%sp(i)%jm,bkarrs%sp(i)%jc/(ctl%v0*r0_cl),&
        bkarrs%sp(i)%rp/r0_cl,bkarrs%sp(i)%ra/r0_cl, r(i))
    end do
end subroutine
subroutine get_v_from_r(spp,r,x,jm, jc, v,vr,vt)
    use com_main_gw
    implicit none
    type(star_pot_para)::spp
    real(8) r, v,phi_tmp,x,vr,vt, jm, jc
    real(8) v2tmp,vr2tmp
    call get_phi_star_full_range(spp,log10(r),phi_tmp)
    v2tmp=(10**phi_tmp+spp%mbh_dmless/r-x)*2
    if(v2tmp<0)then
        print*, "warnning: v2tmp<0, v2tmp=",v2tmp, r, x, 10**phi_tmp, spp%mbh_dmless/r
        v2tmp=abs(v2tmp)
    end if
    v=(v2tmp)**0.5
    vt=jm*jc/r
    vr2tmp=v**2-vt**2
    if(vr2tmp<0)then
        if(vr2tmp<-0.001)then
            print*, "warnning: v2rtmp<0, v2rtmp=", vr2tmp, r, x, v, vt, jm, jc
        end if
        vr2tmp=abs(vr2tmp)
    end if
    vr=vr2tmp**0.5
    
    !print*, "phi_tmp,mbh,r,x,v=",10**phi_tmp,spp%mbh_dmless,r,x,v
end subroutine

subroutine get_rvm_for_sample_stpt(bkarrs, spp,r,vm)
    use com_main_gw
    use md_star_pot
    implicit none
    type(star_pot_para)::spp
    type(particle_samples_arr_type)::bkarrs
    real(8) r(bkarrs%n),vm(bkarrs%n), phi_tmp,vr,vt
    integer i
    do i=1, bkarrs%n
        call get_r_stpt(spp, bkarrs%sp(i)%x,bkarrs%sp(i)%jm,bkarrs%sp(i)%jc/(ctl%v0*r0_cl),&
        bkarrs%sp(i)%rp/r0_cl,bkarrs%sp(i)%ra/r0_cl, r(i))
        call get_v_from_r(spp,r(i),bkarrs%sp(i)%x,bkarrs%sp(i)%jm,bkarrs%sp(i)%jc/(ctl%v0*r0_cl), vm(i),vr,vt)
        !call get_phi_star_full_range(spp,log10(r(i)),phi_tmp)
        !vm(i)=((10**phi_tmp+spp%mbh_dmless/r(i)-)*2)**0.5
    end do
    
end subroutine

subroutine get_r_stpt(spp, e,jm,jc,rp,ra, r)
    use com_sts_type
    use constant
    use ieee_arithmetic
    use model_basic,only:common_aux
    use md_star_pot
    implicit none
    !! ra rp in unit of r0_cl
    type(star_pot_para)::spp
    real(8) phi_tmp, beta_tmp_rp, beta_tmp_ra
    real(8)   e,jm,r
    real(8) s,  rnd, tmpy, vr, rp, ra,phi_out, jc
    real(8) g1, g_1, gr,  vt
    real(8) theta, phi, ag(3), vrout(3), vtout(3), vecvin(3), vecrin(3)
    real(8) grange, sign_of_vr, sign_of_vt, grs,sins2
    !type(s1d_type)::aux

    if(jm.ge.0.99999d0) then
        r=rp
        return
    end if

    call get_beta_full_range(spp, log10(rp), beta_tmp_rp)
    call get_beta_full_range(spp, log10(ra), beta_tmp_ra)
    g_1=2d0**0.5d0*(ra-rp)**0.5*(-(spp%mbh_dmless+beta_tmp_rp)/rp**2+jm**2*jc**2/rp**3)**(-0.5)
    g1=2d0**0.5d0*(ra-rp)**0.5*((spp%mbh_dmless+beta_tmp_ra)/ra**2-jm**2*jc**2/ra**3)**(-0.5)
    grange=max(g1,g_1)

    !print*, "e,jm=",e,jm

    call get_aux_function_for_period_pi2(common_aux,spp, e,jm,jc,rp,ra)
100 s=rnd(0d0,pi/2d0)
    sins2=sin(s)**2
    call common_aux%get_value_s(s, gr)
    
    if(gr>grange+10.0d0)then
        print*, "gr, grange, g_1, g1=",gr,  grange, g_1, g1
        print*, "ex,jm,jc=",e,jm,jc
        print*, "s,ra, rp,beta_rp,beta_ra,sin(s)**2=", s,ra, rp,beta_tmp_ra,beta_tmp_rp,sin(s)**2
        call common_aux%print("aux")
        call get_phi_star_full_range(spp,log10(rp),phi_tmp)
        print*, "v_r(rp)**2=", 2*(10**phi_tmp+spp%mbh_dmless/rp-e)-(jc*jm)**2/rp**2
        call get_phi_star_full_range(spp,log10(ra),phi_tmp)
        print*, "v_r(ra)**2=", 2*(10**phi_tmp+spp%mbh_dmless/ra-e)-(jc*jm)**2/ra**2
        !read(*,*)
    end if
    tmpy=rnd(0d0,grange)
    if(tmpy>gr) goto 100
    r=rp+(ra-rp)*sins2
    
end subroutine
!subroutine set_x2_for_plummer()
!    use com_main_gw
!    implicit none
!    real(8) ml, r0, r1,r2, rho1
!    ml=ctl%plummer_model_mtot
!    r0=ctl%plummer_model_ra_critical
!    r1=10**dms%logrmin
!    r2=10**dms%logrmax
!   ! print*, "r1,r2=",r1,r2
!    call get_plummer_den(ml, r0, dms%logrmin,rho1)
!    sample_emax=ml*r0**2*(1d0/(r1**2+r0**2)**1.5-1d0/(r2**2+r0**2)**1.5)+ &
!                r1**2*rho1*2*pi
!    sample_logemin=log10(sample_emax)
!end subroutine

subroutine init_ctl_rtables()
    use com_main_gw
    implicit none
    integer i
    real(8) weight_tot
    real(8) xb(dms%dstr_bins_r),fx(dms%dstr_bins_r)
    integer i_idx,output_flag
    type(s1d_type)::fden
    
    !real(8) mbh
    !print*,":", ctl%ini_den_model
    !print*, "0"
    call ctl%ini_frho_tot%init(dms%logrmin,dms%logrmax,dms%dstr_bins_r,coeff_sts_type_dc)  ! mass density
    call ctl%ini_frho_tot%set_range()
    call ctl%ini_fphi_tot%init(dms%logrmin,dms%logrmax,dms%dstr_bins_r,coeff_sts_type_dc)
    call ctl%ini_fphi_tot%set_range()
    call ctl%ini_fna_tot%init(dms%logrmin,dms%logrmax,dms%dstr_bins_r,coeff_sts_type_dc)  ! number 
    call ctl%ini_fna_tot%set_range()
    call ctl%ini_fma_tot%init(dms%logrmin,dms%logrmax,dms%dstr_bins_r,coeff_sts_type_dc)
    call ctl%ini_fma_tot%set_range()
    
    ctl%ini_frho_tot%fx=0
    xb=ctl%ini_fphi_tot%xb
    do i=1, ctl%m_bins
        call ctl%ini_fphi(i)%init(dms%logrmin,dms%logrmax,dms%dstr_bins_r,coeff_sts_type_dc)
        call ctl%ini_fphi(i)%set_range()
        i_idx=ctl%ini_model_list_in(i)
        !print*, "i,i_idx=",i,i_idx
        select case(ctl%ini_model_list(i))
            case(ini_den_model_dehnen)
                call get_ini_dehnen_pot(ctl%ini_fphi(i),ctl%dehnen(i_idx))
            case(ini_den_model_plummer)
                call get_ini_plummer_pot(ctl%ini_fphi(i),ctl%plummer(i_idx))
            case default
            !    print*, "error, define ini_den_model", ctl%ini_den_model
            !    stop
        end select
        ctl%ini_fphi_tot%fx=ctl%ini_fphi_tot%fx+ctl%ini_fphi(i)%fx
    end do
    !call ctl%ini_fphi_tot%Print("fphi_tot")
    !print*, "1",ctl%m_bins
    !weight_tot=0
    !do i=1, ctl%m_bins
        !weight_tot=weight_tot+ctl%asymptot_ini(1,i)*ctl%bin_mass(i)
    !end do
    do i=1, ctl%m_bins
        call ctl%ini_frho(i)%init(dms%logrmin,dms%logrmax,dms%dstr_bins_r,coeff_sts_type_dc)
        call ctl%ini_frho(i)%set_range()
        i_idx=ctl%ini_model_list_in(i)
        select case(ctl%ini_model_list(i))
        case(ini_den_model_dehnen)
            call get_ini_dehnen_dens(ctl%ini_frho(i), ctl%dehnen(i_idx)%mtot, &
        ctl%dehnen(i_idx)%ra_crit, ctl%dehnen(i_idx)%gamma)
        case(ini_den_model_plummer)
            call get_ini_plummer_dens(ctl%ini_frho(i), ctl%plummer(i_idx)%mtot, &
        ctl%plummer(i_idx)%ra_crit)
         
        end select
        !call ctl%ini_frho(i)%init_intrp()
        !ctl%ini_frho(i)%fx=log10(ctl%ini_frho(i)%fx)
        ctl%ini_frho_tot%fx=ctl%ini_frho_tot%fx+ctl%ini_frho(i)%fx!*ctl%bin_mass(i)

        !call ctl%ini_frho
        !print*, "ctl%chattery=",ctl%chattery
        !ctl%chattery=2
        if(ctl%chattery.ge.2)then
            !print*, rid
            if(rid.eq.0)then
                print*, "i=",i
                call ctl%ini_frho(i)%print("ini_frho")
            !    stop
            end if
        end if
    end do    
    !ctl%ini_frho_tot%fx=log10(ctl%ini_frho_tot%fx)
    
    do i=1, ctl%m_bins
        call ctl%ini_fna(i)%init(dms%logrmin,dms%logrmax,dms%dstr_bins_r,coeff_sts_type_dc)
        call ctl%ini_fna(i)%set_range()
        i_idx=ctl%ini_model_list_in(i)
        select case(ctl%ini_model_list(i))
        case(ini_den_model_dehnen)
            call get_ini_dehnen_fna(ctl%ini_fna(i), ctl%dehnen(i_idx)%mtot/ctl%bin_mass(i), &
            ctl%dehnen(i_idx)%ra_crit, ctl%dehnen(i_idx)%gamma)
        case(ini_den_model_plummer)
            call get_ini_plummer_fna(ctl%ini_fna(i), ctl%plummer(i_idx)%mtot/ctl%bin_mass(i), &
            ctl%plummer(i_idx)%ra_crit)
 
        end select
        !call ctl%ini_frho(i)%init_intrp()
        ctl%ini_fna_tot%fx=ctl%ini_fna_tot%fx+ctl%ini_fna(i)%fx!*ctl%bin_mass(i)
        !print*, "ctl%chattery=",ctl%chattery
        if(ctl%chattery.ge.2)then
            !print*, rid
            if(rid.eq.0)then
                print*, "i=",i
                call ctl%ini_fna(i)%print("ini_fna")
            end if
        end if
    end do    
    ctl%ini_fma_tot%fx=0
    do i=1, ctl%m_bins
        !call ctl%ini_fma(i)%init(dms%logrmin,dms%logrmax,dms%dstr_bins_r,coeff_sts_type_dc)
        !call ctl%ini_fma(i)%set_range()
        
        !call ctl%ini_frho(i)%init_intrp()
        ctl%ini_fma_tot%fx=ctl%ini_fma_tot%fx+ctl%ini_fna(i)%fx*ctl%bin_mass(i)
        !print*, "ctl%chattery=",ctl%chattery
    end do    
    !call ctl%ini_fma_tot%print("fma_tot")
    !read(*,*)
    !print*, "3" 
    !call ctl%ini_frho_tot%init_intrp()
    !call ctl%ini_frho_tot%print("frho_tot")
    !read(*,*)
end subroutine