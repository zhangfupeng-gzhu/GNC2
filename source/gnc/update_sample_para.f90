subroutine update_sample_para(sample,spp)
    use com_main_gw
    implicit none
    type(particle_sample_type)::sample
    type(star_pot_para)::spp
    real(8) ratio,  ecc_kpl, rp_dm, ra_dm,jc_dm,pd_dm
    real(8) phirp,  ac_dm,logx
    real(8) x, rc, rmax,jc_dmless,p_EJ_dmless, jphlc2
    real(8) rdx,rdy,phi_out
    integer ier, nbin
    if(spp%mbh_dmless.eq.0)then
        call set_ebound_for_samples(sample)
    end if

    sample%x=sample%en/ctl%energy0
    nbin=dms%df_coe_bins
    x=sample%x

    if(ctl%dc_grid_type.eq.dc_grid_irregular)then
        call get_ex_idx_ir(x, sample_table_idx,sample_table_rdx,sample_even)
    else
        call get_ex_idx(x, sample_table_idx,sample_table_rdx,sample_even)
    end if
    !print*, "spp%mbh_dmless=", spp%mbh_dmless
    if(ctl%all_exact.ge.1)then
        call get_sample_para_one(dms,sample,spp)
    else
        if(spp%mbh_dmless.ne.0)then
            logx=log10(x)
            ! call dms%alpha_r%get_value_l(log10(spp%mbh_dmless)-logx-0.3d0, alphax)

            ra_dm=sample%ra/r0_cl
            call get_phi_star_full_range(spp, log10(ra_dm), phi_out)
            sample_alpha=10**phi_out*ra_dm/spp%mbh_dmless

            
            if(sample_alpha<ctl%md_alpha_cri)then
                call get_sample_para_one_kpl(dms,sample)
                ctl%num_get_sample_para_kpl=ctl%num_get_sample_para_kpl+1
            else
                if(sample%rp<sample%r_lc*ctl%md_fine_rtd_range &
                    .and.ctl%include_loss_cone.ge.1)then
                    call get_sample_para_one_appd(dms,sample,spp)
                    ctl%num_get_sample_para_exact=ctl%num_get_sample_para_exact+1
                else
					call get_sample_para_one_grids(dms,sample,spp) 
                end if
            end if
        else 
			if(log10(sample%x)<dms%jc%xb(dms%jc%nbin))then
				call get_sample_para_one_grids(dms,sample,spp)
			else
				call get_sample_para_one(dms,sample,spp)
			end if 
        end if        
    end if
     
    if(ctl%include_loss_cone.ge.1 &
        .and.ctl%gw_radiation_otby.ge.1)then 
        jc_dm=sample%jc/(ctl%v0*r0_cl)
        rp_dm=sample%rp/r0_cl
        ecc_kpl=abs((sample%jm*jc_dm)**2/spp%mbh_dmless/rp_dm-1) 
        ac_dm=rp_dm/(1-ecc_kpl) 
        sample%byot%a_bin=ac_dm*r0_cl; sample%byot%e_bin=ecc_kpl 
    end if
    if(ctl%chattery.ge.4)then
        !if(sample%jm<1d-3)then
            print*
            print*, "start======================================"
    
    endif
    if(ctl%chattery.ge.3)then
            print*, "sample%x,jm,m=",sample%x,sample%jm,sample%m
    end if
    if(ctl%chattery.ge.4)then
        if(spp%mbh_dmless.ne.0)then
            print*, "alpha=",sample_alpha
        end if
        if(ctl%include_loss_cone.ge.1.and.ctl%gw_radiation_otby.ge.1)then
            print*, "keplerian a,e=",sample%byot%a_bin,sample%byot%e_bin
        end if
        !if(alphax<ctl%md_alpha_cri)then
            !   ra_dm=sample%ra/r0_cl
        !    jc_dm=sample%jc/(ctl%v0*r0_cl)
        !end if
        print*, "rpdm, radm, jcdm=",rp_dm, sample%ra/r0_cl, jc_dm
        print*, "sample%period, jph_dmless=",sample%period, sample%jph/(r0_cl*ctl%v0)
        print*, "sample%rp, ra, jc,rp_dmless=",sample%rp, sample%ra, sample%jc,&
            sample%rp/r0_cl
        print*, "sample_idx,idy,rdx,rdy=", sample_table_idx, sample_table_idy, sample_table_rdx, sample_table_rdy
    !    print*, "---- in update samples", rid
    !    read(*,*)
    !end if
    end if

end subroutine

subroutine get_sample_para_one_kpl(dm,sp)
    use com_main_gw
    implicit none
    type(particle_sample_type)::sp
    type(diffuse_mspec)::dm
    real(8) ex, logex, jc, jc_dm,jc_prev
    real(8) rmax, rc,jm,jc_xy,rp_dm,ra_dm!,X
    real(8) pd_xy,p_EJ_dmless
    integer ier
    if(ctl%chattery.ge.4)then
        print*, "==get_sample_para_one_kpl==================="
                !start=======================================
    end if
    ex=sp%x
    logex=log10(ex)
    sp%en=sp%x*ctl%energy0
    !x=sP%x
    jc_prev=sp%jc
    rc=spp_new%mbh_dmless/(2d0*ex)
    jc_dm=(rc*spp_new%mbh_dmless)**0.5
    sp%jc=jc_dm*ctl%v0*r0_cl

    !===========================
    if(ctl%chattery.ge.3)then
        print*, "state_emri_last=", sp%state_emri_last
    end if 
	sp%jm=sp%jph/sp%jc
    
    call set_jm_bound(sp%jm)
    !===========================
    jm=sp%jm
    call get_jm_idx(sp%jm, sample_table_idy,sample_table_rdy,sample_evjum)
    sp%jph=jm*sp%jc

    rp_dm=spp_new%mbh_dmless*(1-(1-jm**2)**0.5)/2d0/ex
    ra_dm=spp_new%mbh_dmless*(1+(1-jm**2)**0.5)/2d0/ex    
    sp%period=rc**1.5/spp_new%mbh_dmless**0.5*pi*2*r0_cl/ctl%v0
    sp%rp=rp_dm*r0_cl
    sp%ra=ra_dm*r0_cl

    if(ctl%chattery.ge.4)then
    	print*, "==end of sample_para_one_kpl==================="
     end if
end subroutine

subroutine get_sample_para_one(dm,sp,spp)
    use com_main_gw
    implicit none
    type(particle_sample_type)::sp
    type(diffuse_mspec)::dm
    type(star_pot_para)::spp
    real(8) ex, logex, jc, jc_dmless
    real(8) rmax, rc,jm,jc_xy,rp_xy,ra_xy
    real(8) pd_xy,p_EJ_dmless_fast,r_c_iter!, jph_dmless
    integer ier
    if(ctl%chattery.ge.4)then
        print*, "==get_sample_para_one==================="
                !start=======================================
    end if
    ex=sp%x
    logex=log10(ex)
    sp%en=sp%x*ctl%energy0 

    call get_rmax_accurate(spp ,  dm%fr_phi, logex,rmax)
 
    rc=r_c_iter(spp, ex,ier)

    jc_xy=jc_dmless(rc,spp)
    !=========================== 
    sp%jc=jc_xy*ctl%v0*r0_cl

    sp%jm=sp%jph/sp%jc 
    call set_jm_bound(sp%jm)
    !===========================
    jm=sp%jm 
    call get_jm_idx(sp%jm, sample_table_idy,sample_table_rdy,sample_evjum)
    sp%jph=jm*sp%jc
    if(jc_xy.eq.0)then
        sp%rp=10**dms%logrmin*r0_cl
        sp%ra=sp%rp
        sp%period=0d0
        return
    end if
    call get_rpra_dmless(spp, ex, jm, jc_xy, &
                log10(rc), rmax, rp_xy,ra_xy)
    sp%rp=rp_xy*r0_cl       
    sp%ra=ra_xy*r0_cl
    !print*, "ex,jm,rp_xy,ra_xy,jc_xy=",ex,jm,rp_xy,ra_xy,jc_xy,rc
    pd_xy=p_EJ_dmless_fast(spp, ex,jm,  jc_xy, rp_xy,ra_xy)
    sp%period=pd_xy*r0_cl/ctl%v0
    if(ctl%chattery.ge.3)then
        print*, "==end of sample_para_one_kpl==================="
        print*, "ex,period,pd_xy, jm=",ex,sp%period,pd_xy, sp%jm
    end if
end subroutine

subroutine get_sample_para_one_appd(dm,sp,spp)
    use com_main_gw
    implicit none
    type(particle_sample_type)::sp
    type(diffuse_mspec)::dm
    type(star_pot_para)::spp
    real(8) ex, logex, jc, jc_dmless
    real(8) rmax, rc,jm,jc_xy,rp_xy,ra_xy
    real(8) pd_xy,p_EJ_dmless,pd_dm,r_c_iter
    integer ier
    if(ctl%chattery.ge.3)then
        print*, "==get_sample_para_one_appd==================="
                !start=======================================
    end if
    ex=sp%x
    logex=log10(ex)
    sp%en=sp%x*ctl%energy0
    
    call get_rmax_accurate(spp ,  dm%fr_phi, logex,rmax) 
    rc=r_c_iter(spp,ex,ier)
    jc_xy=jc_dmless(rc,spp)
    sp%jc=jc_xy*ctl%v0*r0_cl
    !===========================
    sp%jm=sp%jph/sp%jc
    call set_jm_bound(sp%jm)
    !===========================
    
    call get_jm_idx(sp%jm, sample_table_idy,sample_table_rdy,sample_evjum)
    jm=sp%jm
    sp%jph=jm*sp%jc

    if(jc_xy.eq.0)then
        sp%rp=10**dms%logrmin*r0_cl
        sp%ra=sp%rp
        sp%period=0d0
        return
    end if
    call get_rpra_dmless_fast(spp, ex, jm, jc_xy, &
                log10(rc), rmax, rp_xy,ra_xy) 
    sp%rp=rp_xy*r0_cl       
    sp%ra=ra_xy*r0_cl 
    select case(ctl%method_interpolate)
    case(method_int_linear)
        call linear_int_2d_xy(sample_table_idx,sample_table_idy,sample_table_rdx,sample_table_rdy,&
        dm%pd%fxy,dms%df_coe_bins,dms%df_coe_bins,pd_dm)
    case(method_int_nearst)
        pd_dm=dm%pd%fxy(sample_table_idx,sample_table_idy)
    end select 
    sp%period=pd_dm*r0_cl/ctl%v0

end subroutine
    
