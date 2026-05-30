subroutine get_sample_para_one_grids(dm,sp,spp)
     use com_main_gw
     implicit none
     type(particle_sample_type)::sp
     type(diffuse_mspec)::dm
     type(star_pot_para)::spp
     real(8) ex, logex, jc, jc_dm
     real(8) rmax,r_c, rc,jm,rp_dm,ra_dm
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
        case(ebin_type_lin)
            call dm%jc%get_value_l(ex,jc_dm)
        end select
        if(jc_dm<0)then
            print* ,"jc_dm=",logex, jc_dm
            call dm%jc%print("jc_dm")
            stop
        end if
    else
        !jc_dm=dms%jc%fx(sample_table_idx)
        !rc=r_c(spp,ex,ier)
        !print*, "rc=",rc
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
        !print*, "ra_dm=",ra_dm
        !rc=r_c(spp,x,ier)
        !jc_dm=jc_dmless(rc,spp)
        !call get_rmax_accurate(spp ,  dm%fr_phi, logex,rmax)
        !call get_rpra_dmless(spp, ex, jm, jc_dm, &
        !log10(rc), rmax, rp_dm,ra_dm)
        !print*, "rp_dm,ra_dm=",rp_dm,ra_dm

        call linear_int_2d_xy(sample_table_idx,sample_table_idy,sample_table_rdx,sample_table_rdy,&
        common_pd_log%fxy,dms%df_coe_bins,dms%df_coe_bins,pd_dm)
        !print*, "pd_dm=",pd_dm
        !pd_dm=dm%pd%fxy(sample_table_idx,sample_table_idy)
        
        sp%period=10**pd_dm*r0_cl/ctl%v0
        !if(sp%period<0)then
        !    print*, "error! period<0, pd_dm=",sp%period, pd_dm
        !    print*, "idx,idy=",sample_table_idx, sample_table_idy, sample_table_rdx, sample_table_rdy
        !    print*, "pd(idx-1,idy)=", dm%pd%fxy(sample_table_idx-1,sample_table_idy-1),&
        !         dm%pd%fxy(sample_table_idx-1,sample_table_idy),&
        !        dm%pd%fxy(sample_table_idx-1,sample_table_idy+1)
        !    print*, "pd(idx,idy)=", dm%pd%fxy(sample_table_idx,sample_table_idy-1),&
        !         dm%pd%fxy(sample_table_idx,sample_table_idy),&
        !        dm%pd%fxy(sample_table_idx,sample_table_idy+1)
        !        print*, "pd(idx+1,idy)=", dm%pd%fxy(sample_table_idx+1,sample_table_idy-1), &
        !            dm%pd%fxy(sample_table_idx+1,sample_table_idy),&
        !        dm%pd%fxy(sample_table_idx+1,sample_table_idy+1)
        !    print*, "r0_cl,v0=",r0_cl,ctl%v0
        !    stop
        !end if
        !print*, "pd_dm=",pd_dm
        !=================================
        !call get_rmax_accurate(spp ,  dm%fr_phi, logex,rmax)
    
        !print*, "rmax=",rmax
        !call get_rmax_accurate(dm%fphi_star,  dm%fr_phi, logex,rmax)
        !print*, "rmax=",rmax
        !read(*,*)
        !rc=r_c(spp,ex,ier)
        !rc=r_c_iter(spp, ex,ier)

        !jc_dm=jc_dmless(rc,spp)
        
        !call get_rpra_dmless(spp, ex, jm, jc_dm, &
        !            log10(rc), rmax, rp_dm,ra_dm)
        !pd_dm=p_EJ_dmless(spp, ex,sp%jm,  jc_dm, rp_dm,ra_dm)
        !print*, "pd_dm=",pd_dm
        !read(*,*)
        !read(*,*)
        !=================================
    end select
    if(ctl%chattery.ge.5)then
        print*, "==end of get_sample_para_grid==================="
        !block 
        !    use md_star_pot
        !    real(8) r_c_iter,r_c, jc_dmless, rc_dm
        !    rc_dm=r_c_iter(spp_new,ex,ier)
        !    !print*, "rc=", rc, r_c_iter(spp_new,ex,ier), r_c(spp_new,ex,ier)
        !    print*, "ex,logx,jc=", ex,logex, jc_dm, jc_dmless(rc_dm, spp_new)
        !    call dms%jc%print("jc")
        !    read(*,*)
        !    !read(*,*)
        !end block
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
            !print*, "rid,i,rho_rmax=",rid,i,rho_rmax
            !call get_fden_ird_one_log(sample_logrmin,rho_rmin(i),common_gx_ir ,dms,spp)
            call get_fden_ird_one(sample_logrmin,rho_rmin(i),common_gx_ir ,dms,spp,rho_rmax,flag_out)
            if(flag_out.eq.1)then
                if(rid.eq.0)then
                    print*, "in get_rho_rmin, sample_logrmin=", sample_logrmin
                end if
            end if
            !call get_fden_ird_one_fast(sample_logrmin,rho_rmin(i),common_gx_ir ,dms,spp)
        ! print*, "222"
            !if(rid.eq.0)then
            !    print*, "i=",i
            !    call common_gx_ir%print("common_gx_ir")
            !end if
           ! print*, "2:rid,rho_rmin=",rid, rho_rmin(i)
        end do
    case(fden_ana_est_method_2d)
        !do i=1, dms%n
            !common_gx_ir=dms%mb(i)%all%barge_ir
            !call get_none_zero_s1d_ir(dms%mb(i)%all%barge_ir,common_gx_ir)
        ! print*, "111"
            !call get_fden_ird_one(sample_logrmin,rho_rmin(i),common_gx_ir ,dms,spp)
            
           ! call get_fden_ird_2d_one(sample_logrmin,rho_rmin(i),dms%mb(i)%all%gxj_ir,dms, spp)
        ! print*, "222"
            !if(rid.eq.0)then
            !    print*, "i=",i
            !    call common_gx_ir%print("common_gx_ir")
            !end if
       ! end do
        rho_rmin=0
    end select

    spp%spt_rho_rmin=0
    do i=1, dms%n
        spp%spt_rho_rmin=spp%spt_rho_rmin+rho_rmin(i)*dms%mb(i)%mc
    end do
    !call get_fden_ird_one(sample_logrmin,rho_rmin(1), dms%all%all%barge_ir,1d0,dms,spp)
    !spp%spt_rho_rmin=rho_rmin
    !m_r_within_min=0d0
    !do i=1, bks%n
    !    if(bks%sp(i)%x>emax_factor.and.mbh_dmless.eq.0)then
    !        m_r_within_min=m_r_within_min+bks%sp(i)%weight_real/m0_cl
    !    end if
    !end do
    !print*, "m_r_within_min=", m_r_within_min,m_r_within_min/(4*pi/3d0*10**(sample_logrmin*3))
    !spt_rho_rmin=rho_rmin+m_r_within_min/(4*pi/3d0*10**(sample_logrmin*3))
    
    spp%has_set_rhomin=.true.
    !m_r_within_min=m_r_within_min+rho_rmin*4*pi/3d0*10**(sample_logrmin*3)
    if(ctl%fden_ana_est_method.ne.fden_ana_est_method_2d)then
        if(rid.eq.0)then
            print*, "spt_rho_rmin=",spp%spt_rho_rmin
            print*, "rho_rmin(:)=", rho_rmin(:)
        end if
    end if
end subroutine
subroutine get_sample_jlc(ex,mbhin, rt,jc,spp,jlc,ier)
    !use com_main_gw
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
    !    print*, "?????? jphlc2<=0"
    !    print*, "rtd_dm, r0_cl, sp%x, jphlc2, sp%r_lc, star_type_str(sp%obtype)"
    !    print*, rt, r0_cl, x, jphlc2, sp%r_lc, star_type_str(sp%obtype), sp%byot%ms%radius
    !    print*, "mbh_dmless,mbh=",mbh_dmless,mbh
    !    stop
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
    !call rmax%get_value_s(logx,rmax_tmp)
    !call get_rmax_add_accurate(fphi_star,dms%all%all%fmden,rmax,logx,rmax_tmp)

    call get_rmax_accurate(spp,rmax,logenx,rmax_tmp)
   ! print*, "logx, rmax_tmp=",logx, rmax_tmp
    !call get_rmax_add_accurate(fphi_star,dms%all%all%fmden, rmax,logx,rmax_tmp)
    !print*, "rmax_tmp=",rmax_tmp
    !stop
    call my_integral_acc(logrmin,rmax_tmp,fout,1d-24,1d-14, FCN, idid)
    !if(rid.eq.0)then
    !    print*, "logrmin,rmax,fout=",logrmin,rmax_tmp,fout
    !end if
    fx=2**1.5d0*abs(fout)*log(10d0)
contains 
	subroutine FCN(N,X,Y,F,IPAR,RPAR)
		implicit none
		integer n, ipar(100)
		real(8) x, y(n), f(n), rpar(100),ysp, phi_tmp, radius
		!call fphi_tot%get_value_s(x, phi_tmp)
		call get_phi_star_full_range(spp,x,phi_tmp)
        select case(ctl%ebin_type)
        case(ebin_type_log)
            radius=10**x
		    !F(1)=(10**(x*3))*(abs(10**logx-10**phi_tmp-mbh/10**x))**0.5
            F(1)=(radius**2.5d0)*(max((10**phi_tmp-enx)*radius+mbhin,0d0))**0.5
        case(ebin_type_lin)
            radius=x
            F(1)=(radius**2.5d0)*(max((phi_tmp-enx)*radius+mbhin,0d0))**0.5
        end select
        !print*, "x,f(1)=",x,f(1),logx,phi_tmp,mbh
		!if((phi_tmp-10**logx)<-0.1)then
		!	call fphi_star%print("fphi_star")
		!	print*, "logx, x,phi_tmp=",x, 10**logx,phi_tmp
		!	print*, "F(1)=", x, F(1), 10**(x*3), (abs(10**logx-10**phi_tmp-1d0/10**x))**0.5
		!	read(*,*)
		!end if
	end subroutine
end subroutine


subroutine get_rv(  fphi_star, frho, fma, mbh, e,jm,jc,rp,ra, rout,vout)
    use com_sts_type
    use constant
    use ieee_arithmetic
    use model_basic,only:common_aux
    use md_star_pot
    implicit none
    !! ra rp in unit of r0_cl
    type(s1d_type):: fphi_star, fma, frho
    real(8) beta_tmp, beta_tmp_rp, beta_tmp_ra
    real(8)   e,jm,r,mbh
    real(8) s,  rnd, tmpy, vr, rp, ra,phi_out, jc
    real(8) g1, g_1, gr, rout(3), vout(3), vt
    real(8) theta, phi, ag(3), vrout(3), vtout(3), vecvin(3), vecrin(3)
    real(8) grange, sign_of_vr, sign_of_vt, grs,sins2
    !type(s1d_type)::aux

    rout=0; vrout=0
    !print*, "x, j, jc=", e, jm, jc
    !call beta%print("beta")
    !call beta%get_value_s(log10(rp), beta_tmp_rp)
    !call beta%get_value_s(log10(ra), beta_tmp_ra)
    call get_beta_full_range(spp_new, log10(rp), beta_tmp_rp)
    call get_beta_full_range(spp_new, log10(ra), beta_tmp_ra)
    g_1=2d0**0.5d0*(ra-rp)**0.5*(-(mbh+beta_tmp_rp)/rp**2+jm**2*jc**2/rp**3)**(-0.5)
    g1=2d0**0.5d0*(ra-rp)**0.5*((mbh+beta_tmp_ra)/ra**2-jm**2*jc**2/ra**3)**(-0.5)
    !print*, (1+beta_tmp_ra)/ra**2-jm**2*jc**2/ra**3
    !print*, "ra,rp, g_1, g1, beta_rp, beta_ra=", ra, rp, g_1, g1, beta_tmp_rp, beta_tmp_ra
    !read(*,*)
    !call fphi_star%print("fphi_star")
    grange=max(g1,g_1)
    !print*, "grange=", grange

    !call aux%init(0d0,1d0,nbins,sts_type_grid)
	!call aux%set_range()
    call get_aux_function_for_period_pi2(common_aux,spp_new,e,jm,jc,rp,ra)
100 s=rnd(0d0,pi/2d0)
    sins2=sin(s)**2
    call common_aux%get_value_s(s, gr)
    !if(ieee_is_nan(gr)) then
    !    print*, "error!"
    !    stop
    !end if
    !print*, "gr=",gr
    !r=rp+(ra-rp)*sins2
    !call get_phi_star_full_range(fphi_star, log10(r), phi_out)
    !!print*, "phi,e,1/r=", 10**phi_out, e, 1d0/r
    !!print*, "r, ra,rp, vr**2=", r, ra, rp, 2*(10**phi_out+1d0/r-e)-(jm*jc)**2/r**2
    !vr=(2*(10**phi_out+1d0/r-e)-(jm*jc)**2/r**2)**0.5d0
    !!if(s>1d0-1d-3)then
    !!    grs=g1*(ra-rp)**0.5/(ra-r)**0.5
    !!end if
    !!if(s<-1d0+1d-3)then
    !!    grs=g_1*(ra-rp)**0.5/(r-rp)**0.5
    !!end if
    !gr=2*(ra-rp)*(sins2*(1-sins2))**0.5d0/vr
    !print*, "gr=",gr
    !read(*,*)
    if(gr>grange+0.1d0)then
        print*, "gr, grange, g_1, g1=",gr,  grange, g_1, g1
        print*, "r,s,ra, rp,vr,sintheta^2=",r,s,ra, rp,vr,sin(s)**2
        !block
        !    type(s1d_type)::aux
        !    real(8) tmp
        !    call get_aux_function_for_period(aux,fphi_star,beta,e,jm,jc,rp,ra)
        !    !call aux%print("aux")
        !    call aux%get_value_s(sin(s)**2, tmp)
        !    print*, "tmp=",tmp
        !end block
        read(*,*)
    end if
    tmpy=rnd(0d0,grange)
    if(tmpy>gr) goto 100
    r=rp+(ra-rp)*sins2
    vt=jm*jc/r
    call get_phi_star_full_range(spp_new, log10(r), phi_out)
    vr=(2*(10**phi_out+mbh/r-e)-(jm*jc)**2/r**2)**0.5d0

    !print*, "r, vr, vt=", r, vr, vt
    !print*, "vr*r=",vr*r
    !read(*,*)
    vecrin=(/r,0d0,0d0/)
    if((rnd(0d0,1d0)-0.5d0)>0)then
        sign_of_vr=1
    else
        sign_of_vr=-1
    endif
    if((rnd(0d0,1d0)-0.5d0)>0)then
        sign_of_vt=1
    else
        sign_of_vt=-1
    endif
    vecvin=(/sign_of_vr*vr, sign_of_vt*vt, 0d0/)

    theta=rnd(0d0,2*pi)
    call rotation_z(vecrin, rout,theta)
    vecrin=rout
    
    call rotation_z(vecvin, vout,theta)
    vecvin=vout

    phi=rnd(0d0,2*pi)
    call rotation_x(vecrin, rout,phi)
    call rotation_x(vecvin, vout,phi)
    !block
    !    real(8) vecmag
    !    call vector_mag(vecrin, vecmag)
    !    print*, "after vecr=", vecmag
    !end block
    !print*, "vecvin, vout,theta=", vecvin, vout, phi
    !print*, "vec\dot r=",(vecrin(1)*vecvin(1)+vecrin(2)*vecvin(2)+vecrin(3)*vecvin(3))/&
    !    r/(vr**2+vt**2)**0.5
    !read(*,*)
    !read(*,*)
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

    ! if(so%n>ctl%min_sample_in_mass_bin)then
        !print*, "dms_so_get_fslope"
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
        !call so%fslope%print("fslope")
        !read(*,*)
    ! end if
end subroutine

subroutine get_slope0(dm)
    use com_main_gw
    implicit none
    type(diffuse_mspec)::dm
    integer i,j,k,ier
    
    !call dms%mb(1)%star%fden%print("fden before")
    !read(*,*)
    ! do i=1, ctl%ntasks
    !     if(rid.eq.i-1)then
    !         print*, "rid=",rid
    !         do j=1, dm%n
    !             print*, "mass i=",j
    !             call dm%mb(j)%all%fden%print("all fden")
    !             print*, "all%n=",dm%mb(j)%all%n
    !             do k=1, n_tot_comp_sg
    !                 print*, 'type=',k
    !                 print*, "so%n=",dm%mb(j)%dsp(k)%p%n
    !                 if(dm%mb(j)%dsp(k)%p%n>0)then
    !                     call dm%mb(j)%dsp(k)%p%fden%print("fden")
    !                 end if
    !             end do
    !         end do
    !     end if
    !     call mpi_barrier(mpi_comm_world,ier)
    ! end do
    ! stop
    call get_fslope(dm,source_ana)
!print*, "end"
    ! call get_favgm(dm)
    ! call get_favgm2(dm)
end subroutine

!subroutine prepare_barge_ir_tables()
!    use model_basic
!    implicit none
!    integer i,j, n,ier
!    real(8),external:: r_c, pd_dmless,jc_dmless
!    real(8) rmax,rc,x, jm,p_EJ_dmless,ra_dm,rp_dm
!    real(8) xsteps(dms%df_coe_bins)
!    !print*, "stop, finish prepare_barge_ir"
!    !stop
!    n=dms%df_coe_bins
!    call common_jc%init(log10emin_factor,log10emax_factor,n,sts_type_dstr)
!    !call set_gx_ranges_xb(common_jc%xb,xsteps,log10emin_factor,log10emax_factor,n, &
!    !    common_jc%bin_type,dms%jc,dms%pd,gx_func_max_step,gx_func_min_step)
!    !common_jc%xb=dms%mb(1)%star%barge_ir%xb
!    common_barp=dms%rc
!    call get_barp_xy(common_barp%xb,common_barp%fx,common_barp%nbin,&
!        dms%fphi_star,dms%fr_phi,mbh_dmless)
!    call set_nx_ranges_xb(common_jc%xb,common_jc%xsteps,common_jc%xmin,common_jc%xmax,common_jc%nbin,common_jc%bin_type,&
!        common_barp,gx_func_max_step,gx_func_min_step)
!    call common_pd%init(n,n,log10emin_factor,log10emax_factor,log10(jmin_value),log10(jmax_value),sts_type_dstr)
!    call common_pd%set_range()
!
!    common_pd%xcenter=common_jc%xb
!    common_rp=common_pd
!    common_ra=common_pd
!    common_rp%xcenter=common_jc%xb
!    common_ra%xcenter=common_jc%xb
!    do i=1, n
!        x=10**common_jc%xb(i)
!        call get_rmax_accurate(dms%fphi_star,  dms%fr_phi, common_jc%xb(i),rmax)
!        rc= r_c(mbh_dmless,dms%fphi_star,dms%fma_star,x,ier)
!        common_jc%fx(i)=jc_dmless(mbh_dmless,rc,dms%fma_star)
!        do j=1, n
!            jm=10**common_pd%ycenter(j)
!            call get_rpra_dmless(dms%fphi_star, x, jm, common_jc%fx(i), &
!				log10(rc), rmax, rp_dm,ra_dm)
!            common_pd%fxy(i,j)=p_EJ_dmless(dms%fphi_star,dms%all%all%fmden,dms%fma_star,&
!                x,jm,common_jc%fx(i), rp_dm,ra_dm)
!            common_rp%fxy(i,j)=rp_dm
!            common_ra%fxy(i,j)=ra_dm
!        end do
!    end do
!
!
!end subroutine


subroutine set_gx_nx_ranges_ir(dm)
    use md_dms
    use model_basic,only:ctl
    use MPI_comu,only:rid
    use md_coeff,only:ebin_type_lin,ebin_type_log
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
    use md_coeff,only:ebin_type_lin,ebin_type_log
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
    !if(rid.eq.0)then
    !    call cpu_time(t1)
    !end if
    !if(ctl%chattery.ge.1)then
   !     print*, "begin of update_re_tables"
   ! end if
    call get_diffuse_mspec_rtables(dms,spp_new)
    !   call check_rp_ra("at 4")
    call get_ini_ebounds()
    call get_diffuse_mspec_ebound()

    !call set_edstr_bound()
    call init_diffuse_mspec_etables(dms)
    !if(rid.eq.0)then
   !     print*, dms%mb(1)%dc%s2_dee%xcenter,dms%emin,dms%emax
    !end if
    !   call check_rp_ra("at 6")
    !call get_dms_dlnx()
    call get_frphi(spp_new%fphi_star,dms%fr_phi)
    call get_xrc(dms%frc_x,spp_new)
    !call set_barp(dms)
    !if(rid.eq.0.and.ctl%chattery.ge.1)then
    !    call dm%barp_ir%print("update_pre:barp_ir")
    !end if

    !call get_init_etables()
    if(ctl%chattery.ge.2.and.rid.eq.0)then
        print*, "end of update_re_tables"
    end if
    
    !if(rid.eq.0)then
    !    call cpu_time(t2)
    !    print*, "update_re_tables used time:", t2-t1, " s"
    !end if
end subroutine

subroutine set_nx_by_phi(phi,xb,xsteps,nbin,emax)
    !use com_main_gw
    use com_sts_type
    use md_star_pot
    implicit none
    integer n, i, nbin
    type(s1d_type)::phi
    !type(s1d_ird_type)::barp
    real(8) xb(nbin), xsteps(nbin)
    real(8) logr,phi_tmp,emax
    !call phi%print("phi")
    do i=1, nbin
        logr=(phi%xmax-phi%xmin)/real(nbin)*real(i-0.5)+phi%xmin
        call get_phi_star_full_range(spp_new,logr,phi_tmp)
        xb(nbin-i+1)=log10(10**phi_tmp+spp_new%mbh_dmless/10**logr)
        !print*, "logr,xb(i)=", logr, xb(i)
    end do
    xsteps(nbin)=(emax-xb(nbin))*2d0
    !print*, "emin=",emax
    
    do i=nbin-1, 1,-1
        xsteps(i)=(xb(i+1)-xb(i))*2d0-xsteps(i+1)
       ! print*, "xb(i),xsteps(i)=", xb(i),xsteps(i)
    end do
    !read(*,*)
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
        !print*, "barp%nbin=", barp%nbin
        !call barp%print("barp")
        !call set_nx_by_phi(dm%fphi_star,common_barp)
    case(ebin_type_lin)
        call get_barp_xy(log10(barp%xb),barp%fx,barp%nbin, &
        dms%fr_phi,spp%mbh_dmless,spp)
        !call set_nx_by_phi(dm%fphi_star,common_barp)
    end select

end subroutine
!subroutine get_dms_dlnx()
!    use com_main_gw
!    implicit none
!    integer i
!    real(8) mwithin,rtmp,x
!    !type(s1d_ird_type)::dlnphi_dlnr
!    call dms%dlnx%init(sample_logemin,sample_logemax,dms%dstr_bins,sts_type_dstr)
!    call dms%dlnx%set_range()
!    do i=1, dms%dlnx%nbin
!        call get_dms_dlnx_one(dms%dlnx%xb(i),dms%dlnx%fx(i))
!    end do
!    
!    !call dms%dlnx%print("dlnx")
!    call dms%dlnx_ir%init(sample_logemin,sample_logemax,dms%dstr_bins,sts_type_dstr)
!    associate(dlnx_ir=>dms%dlnx_ir)
!        call set_nx_ranges_xb_ir(dlnx_ir%xb,dlnx_ir%xsteps,dlnx_ir%xmin,dlnx_ir%xmax,dlnx_ir%nbin,sts_type_dstr,&
!            dms%dlnx%xb,dms%dlnx%fx,dms%dlnx%nbin,gx_func_max_step,gx_func_min_step)
!    end associate
!    
!    do i=1, dms%dlnx_ir%nbin
!        call get_dms_dlnx_one(dms%dlnx_ir%xb(i),dms%dlnx_ir%fx(i))
!    end do
!    if(rid.eq.0.and.ctl%chattery.ge.1)then
!        call dms%dlnx%print("dlnx")
!        call dms%dlnx_ir%print("dlnx_ir")
!    end if
!    !stop
!end subroutine
!subroutine get_dms_dlnx_one(logx,fx)
!    use com_main_gw
!    implicit none
!    integer i
!    real(8) mwithin,rtmp,x,logx,fx
!    
!    call get_rmax_accurate(dms%fphi_star,dms%fr_phi,logx,rtmp)
!    call get_beta_full_range(dms%fma_star,rtmp,mwithin)
!    x=10**logx
!    fx=(mwithin+mbh_dmless)/(x*10**rtmp)
!   
!end subroutine
subroutine set_barp_step_size(dm, barp_ir,spp)
    use com_main_gw
    implicit none
    integer i
    type(s1d_ird_type)::barp_ir
    type(s1d_type)::barp
    type(diffuse_mspec)::dm
    type(star_pot_para)::spp
    integer nt

    !call dms%barp_ir%init(sample_nxgx_logemin,sample_nxgx_logemax,&
    !dms%dstr_bins_e,sts_type_dstr)
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
    !barp=dms%jc
    !call set_barp_re(barp)
    !call barp%print("barp")
    !if(rid.eq.0.or.rid.eq.15)then
    !    print*, "set_barp_ir"
    !end if

    call set_barp_step_size(dm, barp_ir,spp)
    !barp_ir%xb=dms%dlnx_ir%xb
    !barp_ir%xsteps=dms%dlnx_ir%xsteps
    !call barp_ir%print("barp_ir")
    !call reset_barp_bin_size()
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
        !call set_nx_by_phi(dm%fphi_star,common_barp)
    case(ebin_type_lin)
        call get_barp_xy(log10(barp_ir%xb),barp_ir%fx,barp_ir%nbin, &
        dm%fr_phi,spp%mbh_dmless,spp)
        !call set_nx_by_phi(dm%fphi_star,common_barp)
    end select
end subroutine
subroutine reset_barp_bin_size()
    use com_main_gw
    implicit none
    integer nt
    !type(s1d_ird_type)::barp_ir

    call set_nx_1d_ranges_ir(dms)
    !do i=1, dm%n
        !call set_nx1d_ranges_mb(dm%mb(i),xb,xsteps,n, xmin,xmax)
    !end do
    !call set_nx1d_ranges_mb(dm%all,xb,xsteps,n, xmin,xmax)

    call get_nx_ir_simu(dms)

    call collection_and_avg_fx(dms%all%all%nx_ir%fxw,dms%barp_ir%nbin)

    call reallocate_steps_barp(dms%barp_ir, dms%all%all%nx_ir)

end subroutine
subroutine set_gx_nx_ranges_mb(mb)
    use com_main_gw
    implicit none
    integer i
    type(mass_bins)::mb
    !type(s1d_type)::barp
    
    !barp=dms%jc
    !call set_barp_re(barp)
    !call barp%print("barp")
    do i=1, n_tot_comp
        !call set_gx_ranges_obj(mb%dsp(i)%p%barge_ir)
        associate(b=>mb%dsp(i)%p%nx_ir, d=>mb%dsp(i)%p%barge_ir)
            !call set_nx_ranges_xb(b%xb,b%xsteps, b%xmin,b%xmax,b%nbin,b%bin_type,common_barp,&  
            ! gx_func_max_step,gx_func_min_step)
            !call set_nx_by_phi(dms%fphi_star,b%xb,b%xsteps,b%nbin,b%xmax)
            !call b%init(sample_logemin,sample_logemax,dms%dstr_bins,)

            call b%init(sample_logemin,sample_logemax,dms%dstr_bins_e,use_weight=.true.)

            !call set_nx_ranges_xb(b%xb,b%xsteps, b%xmin, &
            !b%xmax,b%nbin,b%bin_type,barp,&  
            !         gx_func_max_step,gx_func_min_step)

            !call b%print("b")
            !read(*,*)
            call d%init(sample_logemin,sample_logemax,dms%dstr_bins_e,sts_type_dstr)
             b%xb=dms%barp_ir%xb
             b%xsteps=dms%barp_ir%xsteps
             !d%xmin=b%xmin
             !d%xmax=b%xmax
             d%xb=b%xb
             d%xsteps=b%xsteps

        end associate
    end do
    associate(b=>mb%all%nx_ir,d=>mb%all%barge_ir)

        call b%init(sample_logemin,sample_logemax,dms%dstr_bins_e,use_weight=.true.)
        !call set_nx_ranges_xb(b%xb,b%xsteps, b%xmin, &
        !        b%xmax,b%nbin,b%bin_type,barp,&  
        !        gx_func_max_step,gx_func_min_step)
        call d%init(sample_logemin,sample_logemax,dms%dstr_bins_e,sts_type_dstr)
        b%xb=dms%barp_ir%xb
        b%xsteps=dms%barp_ir%xsteps
        !d%xmin=b%xmin
        !d%xmax=b%xmax
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
!    integer n
!    real(8) xb(n),xsteps(n),xmin,xmax
    
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
        !call set_nx_ranges_xb(b%xb,b%xsteps, b%xmin, &
        !b%xmax,b%nbin,b%bin_type,barp,&  
        !         gx_func_max_step,gx_func_min_step)

        !call b%print("b")
        !read(*,*)
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
   ! if(ctl%barge_grid_type.eq.barge_grid_type_regular)then
   !     print*, "ctl%barge_grid_type=regullar, this sub should not be called"
    !    stop
   ! end if
    !print*, dms%barp_ir%xb
    !print*, dms%fr_phi%xb
    !stop
    !dms%frmax%xb=dms%barp_ir%xb
    call dms%dlxb_ir%init(log10(dms%emin),log10(dms%emax),&
        dms%df_coe_bins,sts_type_dstr)
    n=dms%dlxb_ir%nbin
    call set_nx_ranges_xb_ir(dms%dlxb_ir%xb(1:n),dms%dlxb_ir%xsteps(1:n), dms%dlxb_ir%xmin, &
    dms%dlxb_ir%xmax,n,dms%dlxb_ir%bin_type,dms,spp_new,ctl%barge_grid_type)
    !call set_nx_ranges_xb_ir(dms%dlxb_ir%xb(1:n),dms%dlxb_ir%xsteps(1:n), dms%dlxb_ir%xmin, &
    !dms%dlxb_ir%xmax,n,dms%dlxb_ir%bin_type,dms,spp_new,barge_grid_type_iregular_phi)


    !if(rid.eq.0.or.rid.eq.15)then
    !    print*, "1",rid
    !end if
    dms%rc%xb=dms%dlxb_ir%xb
    !if(rid.eq.0)then
    !    print*, "xmin,xmax=",dms%dlxb_ir%xmin,dms%dlxb_ir%xmax
    !    print*, "rc%xb=",dms%rc%xb
    !end if
    !call mpi_barrier(mpi_comm_world,ierr)

    dms%jc%xb=dms%dlxb_ir%xb
    !dms%fr_phi%xb=dms%dlxb_ir%xb
    dms%pd%xcenter=dms%dlxb_ir%xb
    dms%rp%xcenter=dms%dlxb_ir%xb
    dms%ra%xcenter=dms%dlxb_ir%xb        
    !if(rid.eq.0)then
    !    print*, "dms%rp%xcenter=",dms%rp%xcenter
    !end if
    !print*, "2",rid
    !dms%frmax%xmax=dms%dlxb_ir%xmax
    dms%rc%xmax=dms%dlxb_ir%xmax
    dms%jc%xmax=dms%dlxb_ir%xmax
    !dms%fr_phi%xb=dms%dlxb_ir%xb
    dms%pd%xmax=dms%dlxb_ir%xmax
    dms%rp%xmax=dms%dlxb_ir%xmax
    dms%ra%xmax=dms%dlxb_ir%xmax  

    !dms%frmax%xmin=dms%dlxb_ir%xmin
    dms%rc%xmin=dms%dlxb_ir%xmin
    dms%jc%xmin=dms%dlxb_ir%xmin
    dms%pd%xmin=dms%dlxb_ir%xmin
    dms%rp%xmin=dms%dlxb_ir%xmin
    dms%ra%xmin=dms%dlxb_ir%xmin     
    !print*, "e",rid
    !if(rid.eq.0)then
    !    print*, dms%dlxb_ir%xmin,dms%dlxb_ir%xmax
    !    stop
    !endif
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


subroutine get_none_zero_s1d_two_side(barge, barge_nonzero)
    use com_sts_type
    implicit none
    type(s1d_type)::barge, barge_nonzero
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

    
    !if(rid.eq.0.or.rid.eq.15)then
    !    write(*,fmt="(A10,I5, 100F20.10)"), "rid,xb=",rid,xb(1:n)
    !    print*, "end of set nx ranges",rid
    !end if
end subroutine
! subroutine get_estimation_of_nx_from_data()
!     use com_main_gw
!     call set_nx_1d_ranges(dms)
!     call get_nx_simu(dms)
!     call collection_and_avg_fx(dms%all%all%nx%fxw,dms%all%all%nx%nbin)
!     !call dms%all%all%nx%print("nx")
!     !read(*,*)
! end subroutine
subroutine reallocate_steps_barp(barp,nxr)
    use com_main_gw
    implicit none
    type(s1d_ird_type)::barp
    type(s1d_hst_ird_type)::nxr
    integer nt
    real(8) xb_loop(barp%nbin), xstep_loop(barp%nbin), fxw_loop(barp%nbin)
    !print*, "norig,rid=",barp%nbin,rid
    !if(rid.eq.1)then
        
    !end if
    call reallocate_steps(barp%xb,barp%xsteps,barp%nbin,nxr%fxw,&
        xb_loop, xstep_loop,fxw_loop, nt)
    call barp%init(sample_nxgx_logemin,sample_nxgx_logemax,&
        nt,sts_type_dstr)
    barp%xb(1:nt)=xb_loop(1:nt)
    barp%xsteps(1:nt)=xstep_loop(1:nt)
    !print*, "nt,rid=",nt,rid
end subroutine
subroutine reallocate_steps(xb,xstep,n,fxw,xb_loop, xstep_loop,fxw_loop, nt)
    implicit none
    integer n,nt
    real(8) xb(n), xstep(n), fxw(n), xi, xe,x1,x2
    real(8) xb_tmp(n), xstep_tmp(n), fxw_tmp(n)
    real(8) xb_loop(n), xstep_loop(n), fxw_loop(n)
    integer i,j, ns
    integer current_ne
    logical all_nonzero
    nt=0
    !xi=xb(1)-xstep(1)/2d0
    !x1=xb(1)
    
    !if(.not.any(fxw.eq.0d0))then
    !!    print*, "fxw_loop(1:barp%nbin)=", fxw_loop(1:n)
    !    return
    !end if
    all_nonzero=.false.
    !do i=1, n
    !    print*, "xb,xstep,fx=",xb(i),xstep(i),fxw(i),xb(i)-xstep(i)/2d0,xb(i)+xstep(i)/2d0
    !end do
    nt=n
    xb_loop=xb
    xstep_loop=xstep
    fxw_loop=fxw
    do while(.not.all_nonzero)
        all_nonzero=.true.
        ns=0
        i=1
loop1:  do while(i<=nt)
            !i=i+1
            if(i.eq.nt)then
                if(fxw_loop(i).ne.0d0)then
                    ns=ns+1
                    xb_tmp(ns)=xb_loop(i)
                    xstep_tmp(ns)=xstep_loop(i)
                    fxw_tmp(ns)=fxw_loop(i)
                    !if(i+1.eq.nt)then
                    !x1=xb_loop(i)-xstep_loop(i)/2d0
                    !x2=xb_loop(i)+xstep_loop(i)/2d0
                   ! print*, "i,ns,x1,x2=",i,ns,x1,x2, xb_tmp(ns), xstep_tmp(ns), fxw_tmp(ns)
                else
                    ns=ns+1
                    all_nonzero=.false.
                    x1=xb_loop(i-1)-xstep_loop(i-1)/2d0
                    x2=xb_loop(i)+xstep_loop(i)/2d0
                    xb_tmp(ns)=(x1+x2)/2d0
                    xstep_tmp(ns)=(x2-x1)
                    fxw_tmp(ns)=fxw_loop(i)+fxw_loop(i-1)
                   ! print*, "i,ns,x1,x2=",i,ns,x1,x2, xb_tmp(ns), xstep_tmp(ns), fxw_tmp(ns)
                end if    
                exit loop1
            end if
            !print*, "i=",i
            if(fxw_loop(i).ne.0d0.and.fxw_loop(i+1).ne.0d0) then
                ns=ns+1
                xb_tmp(ns)=xb_loop(i)
                xstep_tmp(ns)=xstep_loop(i)
                fxw_tmp(ns)=fxw_loop(i)
                !if(i+1.eq.nt)then
                !x1=xb_loop(i)-xstep_loop(i)/2d0
                !x2=xb_loop(i)+xstep_loop(i)/2d0
                !print*, "i,ns,x1,x2=",i,ns,x1,x2, xb_tmp(ns), xstep_tmp(ns), fxw_tmp(ns)
                ns=ns+1
                xb_tmp(ns)=xb_loop(i+1)
                xstep_tmp(ns)=xstep_loop(i+1)
                fxw_tmp(ns)=fxw_loop(i+1)
                !x1=xb_loop(i+1)-xstep_loop(i+1)/2d0
                !x2=xb_loop(i+1)+xstep_loop(i+1)/2d0
                !print*, "i,ns,x1,x2=",i,ns,x1,x2, xb_tmp(ns), xstep_tmp(ns), fxw_tmp(ns)
                !    exit
                !end if
                !cycle
            else
                ns=ns+1
                all_nonzero=.false.
                x1=xb_loop(i)-xstep_loop(i)/2d0
                x2=xb_loop(i+1)+xstep_loop(i+1)/2d0
                xb_tmp(ns)=(x1+x2)/2d0
                xstep_tmp(ns)=(x2-x1)
                fxw_tmp(ns)=fxw_loop(i)+fxw_loop(i+1)
                !print*, "i,ns,x1,x2=",i,ns,x1,x2, xb_tmp(ns), xstep_tmp(ns), fxw_tmp(ns)
            end if
            i=i+2
        end do loop1
        xb_loop(1:ns)=xb_tmp(1:ns)
        xstep_loop(1:ns)=xstep_tmp(1:ns)
        fxw_loop(1:ns)=fxw_tmp(1:ns)
        !print*, "nt,ns=",nt,ns
        nt=ns        
        !do i=1, ns
        !    print*, "xb,xstep,fx=",xb_tmp(i),xstep_tmp(i),fxw_tmp(i),xb_tmp(i)-xstep_tmp(i)/2d0,&
        !        xb_tmp(i)+xstep_tmp(i)/2d0
        !end do
        !read(*,*)
    end do
    !do i=1, nt
    !    print*, "xb,xstep,fx=",xb_loop(i),xstep_loop(i),xb_loop(i)-xstep_loop(i)/2d0,xb_loop(i)+xstep_loop(i)/2d0
    !end do
    !stop
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

    !real(8) xstep_min,xstep_max
    if(n<1) return
    !print*, "set_nx_ranges_xb_ir: xmin,xmax=",xmin,xmax
    !print*, "dm%fr_phi%xmin,xmax=",dm%fr_phi%xmin,dm%fr_phi%xmax
    !read(*,*)
    !print*, "bg",rid
    !print*, "n=",n
    !call dm%frc_x%print("frc_x")
    call get_xrc(dm%frc_x,spp)
    call get_jc_minmax(xmin,xmax,rcmin,rcmax,jcmin,jcmax,spp)
    !print*, "jcmin,jcmax=",jcmin,jcmax
    !print*, "xmin,xmax, jcmin,jcmax=",xmin,xmax,rcmin,rcmax,jcmin,jcmax
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