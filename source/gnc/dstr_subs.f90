subroutine dms_so_get_fj(so)
    use md_stellar_object
    implicit none
    type(dms_stellar_object)::so
    integer i
    !print*, "so%n=", so%n
    !print*, "size(so%nejw)=", size(so%nejw)
   ! print*, "so%nejw(1:2)%e=", so%nejw(1:2)%e
   ! print*, "so%nejw(1:2)%j=", so%nejw(1:2)%j
   ! print*, "so%nejw(1:2)%w=", so%nejw(1:2)%w
    
   
    !call d2to1(so%nejw(1:so%n)%e,so%nejw(1:so%n)%j**2,so%n,&
    !so%nxj%xmin,so%nxj%xmax,dj_n,so%nxj%ymin,so%nxj%ymax,so%nxj%da_y_proj)
    !do i=1, dj_n
    !    call get_sts_data_weight(so%nxj%da_y_proj(i)%y(:),&
    !    so%nejw(so%nxj%da_y_proj(i)%idx(:))%w,so%nxj%da_y_proj(i)%n,&
    !    so%nxj%ymin,so%nxj%ymax,dj_n2,fc_spacing_linear, so%nxj%fc_y_proj(i))
    !end do

    !call get_distr_of_j(so%nejw(1:so%n)%e,so%nejw(1:so%n)%j,&
    !     so%nejw(1:so%n)%w , so%n, so%barj)
end subroutine

subroutine get_nx_ir_simu(dm)
    use com_main_gw
    implicit none
    integer n,i,j
    type(diffuse_mspec)::dm
    has_made_trade_off=.false.
    dm%all%all%nx_ir%fx=0;dm%all%all%nx_ir%fxw=0
    dm%all%all%nx_ir%nb=0;dm%all%all%nx_ir%nbw=0
    do i=1, dm%n
        do j=1, n_tot_comp
            call dms_get_nx_ir(dm%mb(i)%dsp(j)%p)
            !call dms_get_nx_ir_skip(dm%mb(i)%dsp(j)%p)
        end do
        associate(mb=>dm%mb(i))
            mb%all%nx_ir%fx=0;mb%all%nx_ir%fxw=0
            mb%all%nx_ir%nb=0;mb%all%nx_ir%nbw=0
            do j=1, n_tot_comp
                mb%all%nx_ir%fx=mb%all%nx_ir%fx+mb%dsp(j)%p%nx_ir%fx
                mb%all%nx_ir%fxw=mb%all%nx_ir%fxw+mb%dsp(j)%p%nx_ir%fxw
                mb%all%nx_ir%nb=mb%all%nx_ir%nb+mb%dsp(j)%p%nx_ir%nb
                mb%all%nx_ir%nbw=mb%all%nx_ir%nbw+mb%dsp(j)%p%nx_ir%nbw
            end do
        dm%all%all%nx_ir%fx=dm%all%all%nx_ir%fx+mb%all%nx_ir%fx
        dm%all%all%nx_ir%fxw=dm%all%all%nx_ir%fxw+mb%all%nx_ir%fxw
        dm%all%all%nx_ir%nb=dm%all%all%nx_ir%nb+mb%all%nx_ir%nb
        dm%all%all%nx_ir%nbw=dm%all%all%nx_ir%nbw+mb%all%nx_ir%nbw
        end associate
    end do
    
    do j=1, n_tot_comp
        dm%all%dsp(j)%p%nx_ir%fx=0
        dm%all%dsp(j)%p%nx_ir%fxw=0
        dm%all%dsp(j)%p%nx_ir%nb=0
        dm%all%dsp(j)%p%nx_ir%nbw=0
        do i=1, dm%n
            dm%all%dsp(j)%p%nx_ir%fx=dm%all%dsp(j)%p%nx_ir%fx+dm%mb(i)%dsp(j)%p%nx_ir%fx
            dm%all%dsp(j)%p%nx_ir%fxw=dm%all%dsp(j)%p%nx_ir%fxw+dm%mb(i)%dsp(j)%p%nx_ir%fxw
            dm%all%dsp(j)%p%nx_ir%nb=dm%all%dsp(j)%p%nx_ir%nb+dm%mb(i)%dsp(j)%p%nx_ir%nb
            dm%all%dsp(j)%p%nx_ir%nbw=dm%all%dsp(j)%p%nx_ir%nbw+dm%mb(i)%dsp(j)%p%nx_ir%nbw
        end do
    end do
end subroutine

subroutine get_nx_simu(dm)
    use com_main_gw
    implicit none
    integer n,i,j
    type(diffuse_mspec)::dm
    has_made_trade_off=.false.
    dm%all%all%nx%fx=0;dm%all%all%nx%fxw=0
    dm%all%all%nx%nb=0;dm%all%all%nx%nbw=0
    do i=1, dm%n
        do j=1, n_tot_comp
            call dms_get_nx(dm%mb(i)%dsp(j)%p)
        end do
        associate(mb=>dm%mb(i))
            mb%all%nx%fx=0;mb%all%nx%fxw=0
            mb%all%nx%nb=0;mb%all%nx%nbw=0
            do j=1, n_tot_comp
                mb%all%nx%fx=mb%all%nx%fx+mb%dsp(j)%p%nx%fx
                mb%all%nx%fxw=mb%all%nx%fxw+mb%dsp(j)%p%nx%fxw
                mb%all%nx%nb=mb%all%nx%nb+mb%dsp(j)%p%nx%nb
                mb%all%nx%nbw=mb%all%nx%nbw+mb%dsp(j)%p%nx%nbw
            end do
        dm%all%all%nx%fx=dm%all%all%nx%fx+mb%all%nx%fx
        dm%all%all%nx%fxw=dm%all%all%nx%fxw+mb%all%nx%fxw
        dm%all%all%nx%nb=dm%all%all%nx%nb+mb%all%nx%nb
        dm%all%all%nx%nbw=dm%all%all%nx%nbw+mb%all%nx%nbw
        end associate
    end do
    !if(rid.eq.0)then
    !    call dm%all%all%nx%print("nx_all")
    !end if
end subroutine
subroutine get_fden_sample_particle(dm,so,fden)
	use com_main_gw
	implicit none
	type(s1d_type)::fden
    type(diffuse_mspec)::dm
    type(dms_stellar_object)::so
	integer i,n,j,ier
	real(8) ivr 
    real(8) ivrsum
    !type(s1d_type) aux
    real(8) w,logj,r,r1,r2,pd_xy, ra_xy,rp_xy,jc_xy
    real(8) rc,rmax,r_c, jc_dmless,ex,p_EJ_dmless, jm
    real(8) x1, x2, stmp,yout,t1,t2
    integer n_per_bin(fden%nbin),idid
    logical::first

    fden%fx=0
    n_per_bin=0

   
    do j=1, so%n
        if(mod(j,50000).eq.0) print*, "j=",j
        ex=10**so%nejw(j)%e
        jm=so%nejw(j)%j
        w=so%nejw(j)%w
        pd_xy=so%nejw(j)%pd
        rp_xy=so%nejw(j)%rp
        ra_xy=so%nejw(j)%ra
        jc_xy=so%nejw(j)%jc
        !stmp=0
        first=.true.
loop1:  do i=1, fden%nbin
            !print*, "i=",i
            r=10**fden%xb(i)
            r1=10**(fden%xb(i)-fden%xstep/2d0)
            r2=10**(fden%xb(i)+fden%xstep/2d0)
                 
            if(r1>=ra_xy.or.r2<=rp_xy) then
               ! print*, "r1>ra, r2<rp", r1>ra_xy, r2<rp_xy
                cycle
            end if    
            if(r1<=rp_xy.and.ra_xy<=r2)then
                fden%fx(i)=fden%fx(i)+w
                n_per_bin(i)=n_per_bin(i)+1
                exit loop1
            end if

            if(first)then
                first=.false.
                call get_aux_function_for_period_pi2(common_aux,spp_new,&
                    ex,jm,jc_xy,rp_xy,ra_xy)
            end if

            !call get_r_ana(dm%fphi_star,dm%fma_star,dm%all%all%fmden, &
            !    r,r1,r2, ex,jm,rp_xy,ra_xy,jc_xy,pd_xy, ivr)
            if(rp_xy<r1.and.r2<ra_xy)then
                x1=asin(((r1-rp_xy)/(ra_xy-rp_xy))**0.5); x2=asin(((r2-rp_xy)/(ra_xy-rp_xy))**0.5)
            elseif(r1<rp_xy.and.rp_xy<r2)then
                x1=0; x2=asin(((r2-rp_xy)/(ra_xy-rp_xy))**0.5)
            elseif(r1<ra_xy.and.ra_xy<r2)then
                x1=asin(((r1-rp_xy)/(ra_xy-rp_xy))**0.5); x2=pi/2d0
            end if
            yout=0
            call my_integral_acc(x1,x2,yout,1d-11,1d-6, fcn,idid)
            ivr=yout*2/pd_xy            
            n_per_bin(i)=n_per_bin(i)+1
            fden%fx(i)=fden%fx(i)+ivr*w
            !print*, "r1,rp, ra, r2", r1, rp_xy, ra_xy, r2, ivr
            
            !stmp=stmp+ivr
            !if(rp_xy<1e-2)then
            !    print*, "fx(i)=",i, ivr*w/pd_xy*2
            !end if
            if(pd_xy.eq.0)then
                print*,  "xb(i), fx(i), ivr, w, pd, rp,ra=", 10**fden%xb(i), fden%fx(i), ivr, w, pd_xy, rp_xy, ra_xy
                stop
            end if
        end do loop1
        !print*, "sum=",stmp
        !call fden%print("fden")
        !read(*,*)
        !print*, sum(fden%fx)-stmp, w
        !read(*,*)
        !ivrsum=sum(ivr/dr*w)
        !stmp=sum(fden%fx)
    end do 

    !if(rid.eq.0)then
    !    do i=1, fden%nbin
    !        write(*,fmt="(A20,2E15.6,I13, E15.6)")  "x,fx, n_per_bin=", fden%xb(i), &
    !        fden%fx(i)/(4*pi*10**(fden%xb(i)*2)+pi/3d0*(r2-r1)**2)/m0_cl, &
    !        n_per_bin(i), n_per_bin(i)*10**(fden%xb(i))
    !    end do
    !end if
    do i=1, fden%nbin
        !print*,r1,r2, 4*pi*10**(fden%xb(i)*2),pi/3d0*(r2-r1)**2
        
        if(n_per_bin(i)*10**(fden%xb(i))>ctl%den_bin_cri)then
            r1=10**(fden%xb(i)-fden%xstep/2d0)
            r2=10**(fden%xb(i)+fden%xstep/2d0)
            fden%fx(i)=fden%fx(i)/(4*pi/3d0*(r2**3-r1**3))/m0_cl
        else
            fden%fx(i)=0
        endif
        !read(*,*)
    end do
   
    !call fden%print("fden_simu")
    !read(*,*)
contains
	subroutine fcn(n, x, y, f, par, ipar)
		use, intrinsic :: ieee_arithmetic
		implicit none
		integer n, ipar(100)
		real(8) x, y(n), f(n), par(100)
		real(8) aux_tmp, sintheta,costheta, r,vr_tmp
		!real(8) 

		call common_aux%get_value_s(x, aux_tmp)
		f(1)=aux_tmp
		!print*, "f(1)=",x,f(1)
		if(ieee_is_nan(f(1)))then
			print*, "error! in P_EJ_dmless", x, aux_tmp
			call common_aux%print("aux")
			call spp_new%fphi_star%print("phi_star")
			call dm%all%all%fmden%print("frho")
			call spp_new%fma_star%print("fma")
			stop
		end if
	end subroutine
end subroutine
   
subroutine dms_get_nx(so)
    use md_stellar_object
    use model_basic,only:m0_cl,ctl,dms
    use mpi_comu,only:rid
    use md_star_pot
    implicit none
    type(dms_stellar_object)::so
    real(8) en(so%n), we(so%n), jm(so%n)
    real(8) emin,emax
    integer i
    if(so%n>0)then
        do i=1, so%n
            en(i)=so%nejw(i)%e
            !if(en(i)<ctl%energy_max.and.mbh_dmless.eq.0) en(i)=2*ctl%energy_max-en(i)
            select case(ctl%ebin_type)
            case(ebin_type_log)
                if(en(i)>log10emax_factor.and.spp_new%mbh_dmless.eq.0) then
                    print*, "dms_get_nx:setting en(i) to within:", en(i), " to", log10emax_factor*0.999999,rid
                    en(i)=log10emax_factor-0.00001
                    has_made_trade_off=.true.
                    !cycle
                end if
            case(ebin_type_lin)
                if(en(i)>emax_factor.and.spp_new%mbh_dmless.eq.0) then
                    en(i)=emax_factor*0.999999
                end if
            end select
        end do
        jm(1:so%n)=so%nejw(1:so%n)%j
        we(1:so%n)=so%nejw(1:so%n)%w
        !print*, "en=",en(1:10)
        call so%nx%get_hst(en, we, so%n)
        !call so%nx%print("nx")
        !read(*,*)
        so%nx%fxw=so%nx%fxw/(m0_cl)
        so%nx%fx=so%nx%fx/(m0_cl)
        so%nx%nbw=so%nx%nbw/(m0_cl)
        !so%n_real=sum(we)
    end if
end subroutine
 
subroutine dms_get_nx_ir(so)
    use md_stellar_object
    use model_basic,only:m0_cl,ctl,dms
    use mpi_comu,only:rid
    use md_star_pot
    implicit none
    type(dms_stellar_object)::so
    real(8) en(so%n), we(so%n), jm(so%n)
    real(8) emin,emax
    integer i,idx
    if(so%n>0)then
        do i=1, so%n
            en(i)=so%nejw(i)%e
            we(i)=so%nejw(i)%w
            !if(en(i)<ctl%energy_max.and.mbh_dmless.eq.0) en(i)=2*ctl%energy_max-en(i)
            !select case(ctl%ebin_type)
            !case(ebin_type_log)
            !    if(en(i)>log10emax_factor.and.mbh_dmless.eq.0) then
            !        print*, "dms_get_nx_ir:setting en(i) to within:", i,  en(i), " to", &
            !            log10emax_factor-1d-7*abs(log10emax_factor),current_weight_correct_factor, rid
            !        en(i)=log10emax_factor-1d-7*abs(log10emax_factor)
            !        !we(i)=we(i)*current_weight_correct_factor
            !        !current_weight_correct_factor=current_weight_correct_factor*weight_correct_factor
            !        !idx=i
            !        !has_made_trade_off=.true.
            !        !cycle
            !    end if
            !case(ebin_type_lin)
            !    if(en(i)>emax_factor.and.mbh_dmless.eq.0) then
            !        en(i)=emax_factor*0.999999
            !    end if
            !end select
        end do
        jm(1:so%n)=so%nejw(1:so%n)%j
        
        

        !if(has_made_trade_off)then
        !    !so%nx_ir%xmax=so%nx_ir%xmax+0.01
        !    print*, "en(i)=",en(idx)
        !    print*, "so%nx_ir%xmax, xmax=",so%nx_ir%xmax,&
        !        so%nx_ir%xb(so%nx_ir%nbin)+so%nx_ir%xsteps(so%nx_ir%nbin)/2, &
        !        log10emax_factor*0.999999
        !    call so%nx_ir%get_hst(en, we, so%n)
        !    call so%nx_ir%print("nx_ir bf")
       !
        !    en(1:so%n)=so%nejw(1:so%n)%e
        !    print*, "en(i)=",en(idx)
        !    call so%nx_ir%get_hst(en, we, so%n)
        !    !so%nx_ir%xmax=so%nx_ir%xmax+0.01
        !    call so%nx_ir%print("nx_ir af")
        !    stop
        !else
            call so%nx_ir%get_hst(en, we, so%n)
            !call so%nx_ir%print("nx_ir:1")
            !print*, "so%nx_ir%nbin=", so%nx_ir%nbin,size(so%nx_ir%xsteps),size(so%nx_ir%xb)
            !call get_s1d_ird_hst_weight_kernel(so%nx_ir,en,we,so%n,1)
            !call so%nx_ir%print("nx_ir:2")
            !stop
        !end if



        so%nx_ir%fxw=so%nx_ir%fxw/(m0_cl)
        so%nx_ir%fx=so%nx_ir%fx/(m0_cl)
        so%nx_ir%nbw=so%nx_ir%nbw/(m0_cl)
        !so%n_real=sum(we)
        !block
        !call so%nx_ir%print("nx_ir")
        !call ctl%ini_nx_tot%print("ini_nx")
        !read(*,*)
        !end block
    end if
end subroutine
 
subroutine get_dens_so(dm, spp, so, flag_fden_out)
    use md_dms
    use model_basic
    use md_star_pot
    implicit none
    type(diffuse_mspec)::dm
    type(dms_stellar_object)::so
    real(8)    emin
    integer i,flag_fden_out
    logical allzero
    type(star_pot_para)::spp
    
    !call so%fden%init(dm%logrmin,dm%logrmax,dm%dstr_bins,coeff_sts_type_dc)
    !call so%fden_simu%init(dm%logrmin,dm%logrmax, dm%dstr_bins, coeff_sts_type_dc)
    allzero=all(so%barge_ir%fx.eq.0)
    if(.not.allzero)then
        !select case(ctl%source_fden)
        !case(source_ana)
        !select case(source_ana)
        !case(source_ana)
        !select case (ctl%barge_grid_type)
        !case(barge_grid_type_iregular_phi,barge_grid_type_iregular_barp,barge_grid_type_regular)
            select case(ctl%fden_ana_est_method)
            case(fden_ana_est_method_1d_iso)
                !print*, "spt_rho_rmax=",so%spt_rho_rmax
                call get_fden_ir(dm,spp,so%barge_ir, so%fden, so%spt_rho_rmax, &
                 so%spt_rho_rmin,flag_fden_out)

            !call so%fden%print("fden")
            case(fden_ana_est_method_2d)
                !call get_fden_ir(dm,spp,so%barge_ir, so%fden, n0, &
                !v0, r0_cl, so%spt_rho_rmin)
               ! call so%fden%print("fden2")

                call get_fden_ir_2d(dm,spp,so%gxj_ir, so%fden,  so%spt_rho_rmin)
                !call so%fden%print("fden1")

                !read(*,*)
            end select

        !case(barge_grid_type_regular)
        !    call get_fden(dm,so%barge, so%fden, n0, &
        !        v0, r0_cl, weight_asym)
        !end select
        !case(source_simu)
        if(ctl%source_fden.eq.source_simu)then
            call get_fden_sample_particle(dm,so,so%fden_simu)
        end if
        !end select

        !select case(ctl%source_fden)
        !case(source_ana)
        !    call get_fna(so%fden, so%fNa)
        !case(source_simu)
        !    call get_fna(so%fden_simu, so%fNa)
        !end select

        !so%fden=so%fden_simu
        !print*, "emin=",emin
        !print*, "r0=", log10(r0_cl)
        !call so%fden%print("fden")
        !call so%fden_simu%print("fden_simu")
        !call so%barge%print("barge")
        !read(*,*)
    else
        so%fden%fx=0d0
        so%fden_simu%fx=0d0
    end if
end subroutine


subroutine get_dens_mb_simu(dm,spp, mb, n0, v0 )
    use md_dms
    use md_star_pot
    implicit none
    type(diffuse_mspec)::dm
    type(mass_bins)::mb
    type(star_pot_para)::spp
    real(8) n0, v0,   weight_asym
    integer i
    do i=1, n_tot_comp
        !call get_dens_so_simu(dm, spp, mb%dsp(i)%p, n0, v0, mb%dsp(i)%p%asymp)
        call get_fden_sample_particle(dm,mb%dsp(i)%p,mb%dsp(i)%p%fden_simu)
    end do
end subroutine	

  
subroutine set_gx_ranges_ir(dm)
    use md_dms
    implicit none
    type(diffuse_mspec)::dm
    integer i
    do i=1, dm%n
        call set_gx_ranges_mb(dm%mb(i))
    end do
    call set_gx_ranges_mb(dm%all)
end subroutine
subroutine set_gx_ranges_mb(mb)
    use com_main_gw
    implicit none
    integer i
    type(mass_bins)::mb
    do i=1, n_tot_comp
        !call set_gx_ranges_obj(mb%dsp(i)%p%barge_ir)
        associate(b=>mb%dsp(i)%p%barge_ir)
            call set_gx_ranges_xb(b%xb,b%xsteps, b%xmin,b%xmax,b%nbin,b%bin_type,dms%jc, &
            dms%pd, gx_func_max_step,gx_func_min_step)
        end associate
        !call mb%dsp(i)%p%barge%print("barge")
    end do
    !call set_gx_ranges_obj(mb%all%barge_ir)
    associate(b=>mb%all%barge_ir)
        call set_gx_ranges_xb(b%xb,b%xsteps, b%xmin,b%xmax,b%nbin,b%bin_type,dms%jc, &
            dms%pd, gx_func_max_step,gx_func_min_step)
    end associate
end subroutine
subroutine get_barge_stellar_direct(so,barge)
    use com_main_gw
    implicit none
    type(dms_stellar_object)::so
    type(s1d_type)::barge
    !real(8) sums
    integer i, j,ier
    real(8) x1,x,x2,int_out
    real(8) ex,jm,w,logex
    real(8) pd_xy,jc_xy
    integer n_per_bin(so%barge%nbin)

    if(so%n.eq.0)return
    !call so%barge%print("barge or")
    barge%fx=0
    n_per_bin=0
    do j=1, so%n
        !if(mod(j,50000).eq.0) print*, "j=",j
        logex=so%nejw(j)%e
        ex=10**logex
        jm=so%nejw(j)%j
        w=so%nejw(j)%w
        pd_xy=so%nejw(j)%pd
        jc_xy=so%nejw(j)%jc
        do i=1, barge%nbin
            !print*, "i=",i
            x=10**barge%xb(i)
            x1=10**(barge%xb(i)-barge%xstep/2d0)
            x2=10**(barge%xb(i)+barge%xstep/2d0)
            !print*, "i=",i
            !print*, "barge%xb,xstep=",barge%xb(i),barge%xstep
            !print*, "x, x1,x2=",x,x1,x2
            if(x1<=ex.and.ex<=x2) then
                n_per_bin(i)=n_per_bin(i)+1
                barge%fx(i)=barge%fx(i)+1d0/ex/jc_xy**2/pd_xy*w
            end if                
        end do
    end do
    do i=1, barge%nbin
        !call dms%jc%get_value_s(barge%xb(i),jc_xy)
        !call dms%pd%get_value_d(barge%xb(i),-0.5d0,pd_xy)
        
        !if(n_per_bin(i)*(10**barge%xb(i)*jc_xy**2*pd_xy)**0.5>0.1)then
            barge%fx(i)=barge%fx(i)/(barge%xstep*log(10d0))*pi**(-0.5d0)*2**(-0.5d0)/r0_cl**3/ctl%n0
        !else
        !    barge%fx(i)=0d0
        !end if
    end do
    
    !read(*,*)
    !call barge%print("barge af")
end subroutine
subroutine get_barge_stellar_direct_ird(so,barge)
    use com_main_gw
    implicit none
    type(dms_stellar_object)::so
    type(s1d_ird_type)::barge
    !real(8) sums
    integer i, j,ier,idx, nnz
    real(8) x1,x,x2,int_out
    real(8) ex,jm,w,logex
    real(8) pd_xy,jc_xy
    integer n_per_bin(barge%nbin)
    real(8) xb(barge%nbin), xsteps(barge%nbin), fx(barge%nbin)

    if(so%n.eq.0)return
    !call so%barge%print("barge or")
    barge%fx=0
    n_per_bin=0
    do j=1, so%n
        if(mod(j,50000).eq.0) print*, "j=",j
        logex=so%nejw(j)%e
        ex=10**logex
        jm=so%nejw(j)%j
        w=so%nejw(j)%w
        pd_xy=so%nejw(j)%pd
        jc_xy=so%nejw(j)%jc
        !print*, "ex,jm,w, pd_xy,jc_xy=",ex,jm,w, pd_xy,jc_xy
        do i=1, barge%nbin
            !print*, "i=",i
            x=10**barge%xb(i)
            x1=10**(barge%xb(i)-barge%xsteps(i)/2d0)
            x2=10**(barge%xb(i)+barge%xsteps(i)/2d0)
                 
            if(x1<=ex.and.ex<=x2) then
                n_per_bin(i)=n_per_bin(i)+1
                barge%fx(i)=barge%fx(i)+1d0/jc_xy**2/pd_xy*w
                !call return_idx(log10(jm), log10(jmin_value),log10(jmax_value), dms%df_coe_bins, idx, &
                !    coeff_sts_type_dc)
                !if(idx.ge.common_pd%ny)idx=common_pd%ny
                !barge%fx(i)=barge%fx(i)+1d0/common_jc%fx(i)**2/common_pd%fxy(i,idx)*w
            end if                
            !if(rid.eq.0)then
            !    write(*,fmt="(A20,2F12.6,I5,2F12.6)")  "xb,fx,nbin=",barge%xb(i),barge%fx(i),n_per_bin(i),&
            !        barge%xsteps(i)
            !end if
        end do
    end do
    do i=1, barge%nbin
        !call dms%jc%get_value_s(barge%xb(i),jc_xy)
        !call dms%pd%get_value_d(barge%xb(i),-0.5d0,pd_xy)
        x1=10**(barge%xb(i)-barge%xsteps(i)/2d0)
        x2=10**(barge%xb(i)+barge%xsteps(i)/2d0)
        !if(n_per_bin(i)*(10**barge%xb(i)*jc_xy**2*pd_xy)**0.5>0.1)then
            barge%fx(i)=barge%fx(i)/(x2-x1)*pi**(-0.5d0)*2**(-0.5d0)/r0_cl**3/ctl%n0
        !else
        !    barge%fx(i)=0d0
        !end if
    end do

    !do j=1, ctl%ntasks
    !    if(rid.eq.j-1)then
    !        print*, "j=",j
    !        do i=1, barge%nbin
    !            write(*,fmt="(A20,2F12.6,I5,2F12.6)")  "xb,fx,nbin=",barge%xb(i),barge%fx(i),n_per_bin(i)
    !        end do
    !    end if
    !    call mpi_barrier(mpi_comm_world,ier)
    !end do
    !read(*,*)
    !call barge%print("barge af")
end subroutine

subroutine get_barge0(dm)
    use model_basic
    implicit none
    type(diffuse_mspec)::dm
    integer i,j 
    select case(ctl%barge_evl_method)
    case(barge_evl_method_grid_2d) 
            dm%all%all%barge_ir%fx=0
            do i=1, dm%n
                associate(mb=>dm%mb(i))
                    mb%all%barge_ir%fx=0
                    do j=1, n_tot_comp
                        call get_barge_stellar_2d_ir(mb%dsp(j)%p,dm%jbin_type)
                        mb%all%barge_ir%fx=mb%all%barge_ir%fx+mb%dsp(j)%p%barge_ir%fx
                        !read(*,*)
                    end do
                end associate
                associate(all=>dm%all)
                    all%all%barge_ir%fx=all%all%barge_ir%fx+dm%mb(i)%all%barge_ir%fx*dms%mb(i)%mc
                end associate
            end do
      

    case(barge_evl_method_direct)
        dm%all%all%barge%fx=0
        do i=1, dm%n
            associate(mb=>dm%mb(i))
                mb%all%barge%fx=0
                do j=1, n_tot_comp
                    call get_barge_stellar_direct(mb%dsp(j)%p,mb%dsp(j)%p%barge)
                    mb%all%barge%fx=mb%all%barge%fx+mb%dsp(j)%p%barge%fx
                    !read(*,*)
                end do
            end associate
            associate(all=>dm%all)
                all%all%barge%fx=all%all%barge%fx+dm%mb(i)%all%barge%fx*dms%mb(i)%mc
            end associate
        end do
        !call prepare_barge_ir_tables()
        dm%all%all%barge_ir%fx=0
        do i=1, dm%n
            associate(mb=>dm%mb(i))
                mb%all%barge_ir%fx=0
                do j=1, n_tot_comp
                    call get_barge_stellar_direct_ird(mb%dsp(j)%p,mb%dsp(j)%p%barge_ir)
                    mb%all%barge_ir%fx=mb%all%barge_ir%fx+mb%dsp(j)%p%barge_ir%fx
                    !read(*,*)
                end do
            end associate
            associate(all=>dm%all)
                all%all%barge_ir%fx=all%all%barge_ir%fx+dm%mb(i)%all%barge_ir%fx*dms%mb(i)%mc
            end associate
        end do


    case default
        print*, "error! defile barge method"
        stop
    end select
    
end subroutine
  
subroutine get_barge_stellar_ir(so)
    use md_dms
    use model_basic,only:dms,ctl,ctl
    use MPI_comu, only:rid,mpi_comm_world
    use md_coeff
    implicit none
    type(dms_stellar_object)::so
    !real(8) sums
    integer i, j, idid,ierr
    real(8) x,int_out,sums,x1,x2,x3
    real(8) f1,f2,f3,barp_avg
    !type(s1d_type)::barp

    if(so%n.eq.0)return
    !barp=so%barge
    !barp%xb=so%barge_ir%xb
        
    !call get_barp_xy(barp%xb, barp%fx, barp%nbin,dms%fphi_star,dms%fr_phi,mbh_dmless)
    !do i=1, ctl%ntasks
    !    if(rid.eq.i-1)then
    !        print*, "nx_ir%xsteps=",so%nx_ir%xsteps(1:10)
    !    end if
    !    call mpi_barrier(mpi_comm_world,ierr)
    !end do
    !call dms%barp_ir%print("barp_ir")
    do i=1, so%barge_ir%nbin
        int_out=0
        !print*,"xmax=",mb%barge%xmax
        if(so%barge_ir%nbin.ne.so%nx_ir%nbin.or.so%barge_ir%nbin.ne.dms%barp_ir%nbin)then
            print*, "error! barge%nbin should = so%nx_ir"
            stop
        end if
        !sums=0
        !x1=so%nx_ir%xb(i)-so%nx_ir%xsteps(i)/2d0
        !call barp%get_value_l(x1,f1)
        !call get_barp_one(x1,f1,dms%fphi_star,dms%fr_phi,mbh_dmless)
        !x2=so%nx_ir%xb(i)
        !f2=barp%fx(i)
        !call get_barp_one(x2,f2,dms%fphi_star,dms%fr_phi,mbh_dmless)
        !x3=so%nx_ir%xb(i)+so%nx_ir%xsteps(i)/2d0
        !call get_barp_one(x3,f3,dms%fphi_star,dms%fr_phi,mbh_dmless)
        !call barp%get_value_l(x3,f3)
        
        !sums=sums+(f1*0.5+f2+f3*0.5)/2d0
        !print*, "sums=",barp%xb(i),sums,barp%fx(i)
        !read(*,*)
        !so%barge_ir%fx(i)=2**(-1.5d0)*pi**(-0.5)*so%nx_ir%fxw(i)/sums/(10**so%nx_ir%xb(i)*log(10d0))
        select case(ctl%ebin_type)
        case(ebin_type_log)
            call dms%barp_ir%get_value_s(so%nx_ir%xb(i), barp_avg)
            !if(rid.eq.0)then
                !print*, "so%nx_ir%xb(i),barp_avg=", so%nx_ir%xb(i),barp_avg
                !call dms%barp_ir%print('barp_ir')
            !end if
            if(barp_avg.ne.0)then
                !call dms%barp_ir%get_value_s(so%nx_ir%xb(i)-so%nx_ir%xsteps(i)/2,f1)
                !f2=dms%barp_ir%fx(i)
                !call dms%barp_ir%get_value_s(so%nx_ir%xb(i)+so%nx_ir%xsteps(i)/2,f3)
                !barp_avg=(f1*0.5+f2+f3*0.5)/2d0
                !barp_avg=dms%barp_ir%fx(i)
                so%barge_ir%fx(i)=2**(-1.5d0)*pi**(-0.5)*so%nx_ir%fxw(i)&
                    /barp_avg/(10**so%nx_ir%xb(i)*log(10d0))   
               ! print*, so%nx_ir%xb(i), barp_avg, f2, so%nx_ir%fxw(i)
            else
                so%barge_ir%fx(i)=0d0
            end if
        case(ebin_type_lin)
            print*, "stop    xx"
            stop
        end select
        !print*, "fx=",so%barge_ir%fx(i)
        !so%barge_ir%fx(i)=2**(-1.5d0)*pi**(-0.5)*so%nx_ir%fxw(i)/barp%fx(i)/&
        !    ((10**(so%nx_ir%xb(i)+so%nx_ir%xsteps(i)/2d0)-10**(so%nx_ir%xb(i)-so%nx_ir%xsteps(i)/2d0))/&
        !    (so%nx_ir%xsteps(i)))
        !print*, "fx=",so%barge_ir%fx(i)            
    end do
    !call so%barge_ir%print("barge_ir")
    !read(*,*)
    call so%barge_ir%prepare_spline()
 
    !call so%nx%print("nx")
    !if(rid.eq.0)then
    !    call dms%fphi_star%print("phi")
    !    call dms%fr_phi%print("fr_phi")
    !    call barp%print("barp")
    !end if
    !read(*,*)
end subroutine
  

subroutine get_fden_ir(dm,spp,gx, fden, rho_rmax,rho_rmin,flag_fden_out)
	use constant
	use com_sts_type
	use my_intgl
    use md_dms
    use model_basic,only:ctl,sample_logrmin
    use md_star_pot
    use md_coeff
    use MPI_comu,only:rid,mpi_comm_world
	use, intrinsic :: ieee_arithmetic
    implicit none
    type(s1d_ird_type)::common_gx_ir
    type(diffuse_mspec)::dm
	type(s1d_type)::fden
    type(s1d_ird_type)::gx
	real(8) w_asymp,n0, rho_rmin, xb(gx%nbin),gx_out,rho_rmax
	integer i,ierr,flag_out,flag_fden_out
    integer::debug=0
    type(star_pot_para)::spp
    interface 
        subroutine get_fden_ird_one_log(logr, fx,common_gx_ir,dm,spp)
            use com_main_gw
            implicit none
            real(8) logr, fx
            type(s1d_ird_type),target::common_gx_ir
            type(diffuse_mspec)::dm
            type(star_pot_para)::spp
        end subroutine
        subroutine get_fden_ird_one(logr, fx,common_gx_ir,dm,spp,rho_rmax,flag_out)
            use com_main_gw
            implicit none
            real(8) logr, fx,rho_rmax
            type(s1d_ird_type),target::common_gx_ir
            type(diffuse_mspec)::dm
            type(star_pot_para)::spp
            integer flag_out
        end subroutine 
    end interface
    !print*, "using fen ir"
    common_gx_ir=gx
    flag_fden_out=0
	do i=1, fden%nbin
		!call get_fden_ird_one_log(fden%xb(i),fden%fx(i),common_gx_ir,dm,spp)
        call get_fden_ird_one(fden%xb(i),fden%fx(i),common_gx_ir,dm,spp,rho_rmax,flag_out)
        if(flag_out.eq.1)then
            if(rid.eq.0)then
                print*, "in bin i,xb=",i, fden%xb(i)
            end if
            flag_fden_out=1
        end if
        
	end do
    !print*, "fden%xb(2), fx(2),rid=",fden%xb(2),fden%fx(2),rid
    !do i=1, ctl%ntasks
    !    if(rid.eq.i-1)then
          !  call fden%print("fden")
    !    end if
    !    call mpi_barrier(mpi_comm_world,ierr)
    !end do
    !read(*,*)
   ! call get_fden_ird_one_log(sample_logrmin,rho_rmin,common_gx_ir,dm,spp)
    call get_fden_ird_one(sample_logrmin,rho_rmin,common_gx_ir,dm,spp,rho_rmax,flag_out)
    if(flag_out.eq.1)then
        if(rid.eq.0)then
            print*, "in sample_logrmin=",sample_logrmin
        end if
        flag_fden_out=1
    end if
    !call get_fden_ird_one_fast(sample_logrmin,rho_rmin,common_gx_ir,dm,spp)
    !print*, "rid,rho_rmin=",rid,rho_rmin
    !print*, "rho_rmax=",rho_rmax
   !
    !read(*,*)
end subroutine
subroutine get_fden_ird_one(logr, fx,common_gx_ir,dm,spp,rho_rmax,flag_out)
    use com_main_gw
    implicit none
    real(8) logr, fx, int_out,phi_tot,exmax,exmin
    real(8) xmin,xmax,rho_rmax
    type(s1d_ird_type),target::common_gx_ir
   ! type(s1d_ird_type)::gx_ir_tmp
    integer idid,ierr,debug
    type(diffuse_mspec)::dm
    type(star_pot_para)::spp
    integer i,flag_out

    int_out=0
    debug=0
    flag_out=0
    !fc_ir_share=>common_gx_ir
    !xmin=10**fc_ir_share%xmin
    !xmax=10**fc_ir_share%xmax
!    if(rid.eq.0)then
!        call common_gx_ir%print("common_gx_ir")
        !stop
!    end if
    !gx_ir_tmp=common_gx_ir
    !!
    !do i=1, gx_ir_tmp%nbin
    !    if(common_gx_ir%fx(i)>0)then
    !        gx_ir_tmp%fx(i)=log10(common_gx_ir%fx(i))
    !    else
    !        gx_ir_tmp%fx(i)=-100
    !    end if
    !end do
    !call gx_ir_tmp%prepare_spline(0d0,0d0)
    !!call gx_ir_tmp%print("gx_ir_tmp")
    !call gx_ir_tmp%print("gx_ir_tmp")

    select case(ctl%ebin_type)
    case(ebin_type_log)
        exmin=10**(common_gx_ir%xmin)
    case(ebin_type_lin)
        exmin=common_gx_ir%xmin
    case default
        print*, "get_fden_ird_one:stop"
        stop
    end select

    call get_phi_star_full_range(spp, logr, exmax)
    !call get_dms_starpt_one(fden%xb(i),xmax,dms%all%all%fmden)
    select case(ctl%ebin_type)
    case(ebin_type_log)
        phi_tot=min(10**exmax+spp%mbh_dmless/10**logr,10**common_gx_ir%xmax)
        !if(rid.eq.0)then
        !    print*, "phi_tot=", phi_tot, 10**exmax+mbh_dmless/10**logr,10**common_gx_ir%xmax
        !end if
        !phi_tot= 10**exmax+mbh_dmless/10**fden%xb(i) 
    case(ebin_type_lin)
        phi_tot=min(exmax+spp%mbh_dmless/10**logr, common_gx_ir%xmax)
        !phi_tot=exmax+mbh_dmless/10**fden%xb(i) 
    end select
    !call common_gx_ir%print("gx")
    !debug=1
    !print*, "i=",i
    !print*, "exmin+1d-12<phi_tot=",exmin+1d-12<phi_tot, exmin,phi_tot
    if(exmin+1d-12<phi_tot)then
        !debug=1
        !==test====
        !phi_tot=1
        !common_gx_ir%xmin=-2d0
        !common_gx_ir%xmax=0d0
        !!call common_gx_ir%set_range()
        !call set_range(common_gx_ir%xb,common_gx_ir%nbin,common_gx_ir%xmin,common_gx_ir%xmax,sts_type_dstr)
        !do i=1, common_gx_ir%nbin
        !    common_gx_ir%fx(i)=10**(common_gx_ir%xb(i)*2)
        !end do
        !call my_integral_acc(0.01d0,phi_tot, int_out,&
        !    fr_funcs_int_acc_a,fr_funcs_int_acc_r, fcn,idid)
        !print*, "int_out=",int_out
        !block
        !    real(8) t,fth
        !    t=(1-0.01)**0.5
        !    fth=2d0/3d0*t**3-4d0/5d0*t**5+2d0/7d0*t**7
        !    print*, "fth, d=", fth, (fth-int_out)/fth
        !end block
        !
        !stop
        !-----------
        call my_integral_acc(exmin,phi_tot, int_out,&
            fr_funcs_int_acc_a,fr_funcs_int_acc_r, fcn,idid)
            !read(*,*)
        if(idid<0)then
            if(rid.eq.0)then
                print*, "get_fden_ird_one:xmin,phi_tot,int_out=",exmin,phi_tot,int_out
                print*, "logr,xmax=",logr,exmax,spp%mbh_dmless/10**logr
                print*, "idid=",idid
                print*, "sample_rmin=",sample_logrmin
                call common_gx_ir%print("common_gx_ir")
            end if
            ! debug=1
            ! call my_integral_acc(exmin,phi_tot, int_out,&
            !     fr_funcs_int_acc_a,fr_funcs_int_acc_r, fcn,idid)
            flag_out=1
            ! return
        end if
    else
        int_out=0d0
    end if
    !call my_integral_acc(xmin,phi_tot,int_out, fcn,idid)
    
    fx=max(2/pi**0.5*int_out,0d0)+rho_rmax !*n0 ! in unit of r0_cl^{-3}
    !print*, "r,fx=",logr,fx
    !print*, "xmin,xmax,r=",10**xmin,phi_tot,10**fden%xb(i), fden%fx(i)
    !if(fden%fx(i).eq.0)then
     !   read(*,*)
    !end if
    !if(rid.eq.0)then
    !    print*, fden%fx(i), int_out,xmax,xmin,phi_tot
    !end if
    !read(*,*)

    !do j=1, ctl%ntasks
    !    if(rid.eq.j-1.and.ctl%chattery.ge.1)then
    !        print*, "i=",j
    !        print*, "xmin,phi_tot=", xmin,phi_tot,int_out
    !    end if
    !    call mpi_barrier(mpi_comm_world,ierr)
    !end do
contains
	subroutine fcn(n, x, y, f, par, ipar)
		implicit none
		integer n, ipar(100)
		real(8) x, y(n), f(n), par(100), yout,xin

		select case(ctl%ebin_type)
        case(ebin_type_log)
            xin=max(log10(x),common_gx_ir%xmin )
        case(ebin_type_lin)
            xin=max(x,exmin )
        end select
        !print*, "1"
        !call get_gx_full_range_ir_log(gx_ir_tmp,xin,yout)
        call get_gx_full_range_ir(common_gx_ir,xin,yout)
        !=====================
        !xin=x        
        !call fgx_mb_ir(xin,yout,xmin,xmax)
        !=====================
       ! print*, "x,xmin,yout2=",x,xin,yout
        !read(*,*)
        !call get_phi_star_full_range(dms%fphi_star,)
        yout=max(yout,0d0)
		f(1)=yout*sqrt(max(phi_tot-x,0d0))
        if(debug.eq.1)then
            write(*,fmt="(A20,4E30.20)")  "x, yout, f=", x,xin, yout, f(1)
        endif
        !write(*,fmt="(A20,5E32.20,I5,L5)"), "x,f,yout=",x,xin,f(1),yout,common_gx_ir%xmin, &
        !    rid,xin>common_gx_ir%xmin
        !if(f(1).eq.0)then
        !    print*, "xmin=",xmin
        !    call common_gx_ir%print("common_gx_ir")
        !    read(*,*)
        !end if
		if(phi_tot-x<-1e-1)then
			print*, "error, r0_cl/r-x<0"
            print*, "10**exmax,1d0/r, r0_cl, logr, x=", &
                10**exmax, 1d0/10**logr, r0_cl, logr, x
			stop
		end if
		if(ieee_is_nan(f(1)).or..not.(ieee_is_finite(f(1))))then
			print*, "in fden"
			print*, "x, f=", x, f(1), yout, r0_cl/10**logr-x, r0_cl, logr
			call common_gx_ir%print("common_gx")
			stop
		end if
		!call gx%print("gx")
		
		!read(*,*)
	end subroutine
end subroutine
 

subroutine get_fden_ird_one_log(logr, fx,common_gx_ir,dm,spp)
    use com_main_gw
    implicit none
    real(8) logr, fx, int_out,phi_tot,exmax,exmin
    real(8) xmin,xmax
    type(s1d_ird_type),target::common_gx_ir
   ! type(s1d_ird_type)::gx_ir_tmp
    integer idid,ierr,debug
    type(diffuse_mspec)::dm
    type(star_pot_para)::spp
    integer i

    int_out=0
    debug=0
    !fc_ir_share=>common_gx_ir
    !xmin=10**fc_ir_share%xmin
    !xmax=10**fc_ir_share%xmax
!    if(rid.eq.0)then
!        call common_gx_ir%print("common_gx_ir")
        !stop
!    end if
    !gx_ir_tmp=common_gx_ir
    !!
    !do i=1, gx_ir_tmp%nbin
    !    if(common_gx_ir%fx(i)>0)then
    !        gx_ir_tmp%fx(i)=log10(common_gx_ir%fx(i))
    !    else
    !        gx_ir_tmp%fx(i)=-100
    !    end if
    !end do
    !call gx_ir_tmp%prepare_spline(0d0,0d0)
    !!call gx_ir_tmp%print("gx_ir_tmp")
    !call gx_ir_tmp%print("gx_ir_tmp")

    select case(ctl%ebin_type)
    case(ebin_type_log)
        exmin=10**(common_gx_ir%xmin)
    case(ebin_type_lin)
        exmin=common_gx_ir%xmin
    case default
        print*, "get_fden_ird_one:stop"
        stop
    end select

    call get_phi_star_full_range(spp, logr, exmax)
    !call get_dms_starpt_one(fden%xb(i),xmax,dms%all%all%fmden)
    select case(ctl%ebin_type)
    case(ebin_type_log)
        phi_tot=min(10**exmax+spp%mbh_dmless/10**logr,10**common_gx_ir%xmax)
        !if(rid.eq.0)then
        !    print*, "phi_tot=", phi_tot, 10**exmax+mbh_dmless/10**logr,10**common_gx_ir%xmax
        !end if
        !phi_tot= 10**exmax+mbh_dmless/10**fden%xb(i) 
    case(ebin_type_lin)
        phi_tot=min(exmax+spp%mbh_dmless/10**logr, common_gx_ir%xmax)
        !phi_tot=exmax+mbh_dmless/10**fden%xb(i) 
    end select
    !call common_gx_ir%print("gx")
    !debug=1
    !print*, "i=",i
    !print*, "xmin,xmax,r=",10**xmin,phi_tot,10**fden%xb(i), fden%xb(i)
    if(exmin+1d-12<phi_tot)then
        !debug=1
         !==test====
        !phi_tot=1
        !common_gx_ir%xmin=-2d0
        !common_gx_ir%xmax=0d0
        !!call common_gx_ir%set_range()
        !call set_range(common_gx_ir%xb,common_gx_ir%nbin,common_gx_ir%xmin,common_gx_ir%xmax,sts_type_dstr)
        !do i=1, common_gx_ir%nbin
        !    common_gx_ir%fx(i)=common_gx_ir%xb(i)*2
        !end do
        !call my_integral_acc(0.01d0,1d0, int_out,&
        !    fr_funcs_int_acc_a,fr_funcs_int_acc_r, fcn,idid)
        !print*, "int_out:log=",int_out
        !block
        !    real(8) t,fth
        !    t=(1-0.01)**0.5
        !    fth=2d0/3d0*t**3-4d0/5d0*t**5+2d0/7d0*t**7
        !    print*, "fth, d=", fth, (int_out-fth)/fth
        !end block
        !
        !stop

        call my_integral_acc(exmin,phi_tot, int_out,&
            fr_funcs_int_acc_a,fr_funcs_int_acc_r, fcn,idid)
            !read(*,*)
        if(idid<0)then
            print*, "get_fden_ird_one_log:xmin,phi_tot,int_out=",exmin,phi_tot,int_out
            print*, "logr,xmax=",logr,exmax,spp%mbh_dmless/10**logr
            debug=1
            call my_integral_acc(exmin,phi_tot, int_out,&
                fr_funcs_int_acc_a,fr_funcs_int_acc_r, fcn,idid)
            stop
        end if
    else
        int_out=0d0
    end if
    !call my_integral_acc(xmin,phi_tot,int_out, fcn,idid)
    
    fx=max(2/pi**0.5*int_out,0d0) !*n0 ! in unit of r0_cl^{-3}
    !print*, "r,fx=",logr,fx
    !print*, "xmin,xmax,r=",10**xmin,phi_tot,10**fden%xb(i), fden%fx(i)
    !if(fden%fx(i).eq.0)then
     !   read(*,*)
    !end if
    !if(rid.eq.0)then
    !    print*, fden%fx(i), int_out,xmax,xmin,phi_tot
    !end if
    !read(*,*)

    !do j=1, ctl%ntasks
    !    if(rid.eq.j-1.and.ctl%chattery.ge.1)then
    !        print*, "i=",j
    !        print*, "xmin,phi_tot=", xmin,phi_tot,int_out
    !    end if
    !    call mpi_barrier(mpi_comm_world,ierr)
    !end do
contains
	subroutine fcn(n, x, y, f, par, ipar)
		implicit none
		integer n, ipar(100)
		real(8) x, y(n), f(n), par(100), yout,xin

		select case(ctl%ebin_type)
        case(ebin_type_log)
            xin=max(log10(x),common_gx_ir%xmin )
        case(ebin_type_lin)
            xin=max(x,exmin )
        end select
        !print*, "1"
        call get_gx_full_range_ir_log(common_gx_ir,xin,yout)
        !call get_gx_full_range_ir(common_gx_ir,xin,yout)
        !=====================
        !xin=x        
        !call fgx_mb_ir(xin,yout,xmin,xmax)
        !=====================
       ! print*, "x,xmin,yout2=",x,xin,yout
        !read(*,*)
        !call get_phi_star_full_range(dms%fphi_star,)
        yout=max(yout,0d0)
		f(1)=yout*sqrt(max(phi_tot-x,0d0))
        if(debug.eq.1)then
            write(*,fmt="(A20,20E30.20)")  "x, yout, f=", x,xin, yout, f(1),y
        endif
        !write(*,fmt="(A20,5E32.20,I5,L5)"), "x,f,yout=",x,xin,f(1),yout,common_gx_ir%xmin, &
        !    rid,xin>common_gx_ir%xmin
        !if(f(1).eq.0)then
        !    print*, "xmin=",xmin
        !    call common_gx_ir%print("common_gx_ir")
        !    read(*,*)
        !end if
		if(phi_tot-x<-1e-1)then
			print*, "error, r0_cl/r-x<0"
            print*, "10**exmax,1d0/r, r0_cl, logr, x=", &
                10**exmax, 1d0/10**logr, r0_cl, logr, x
			stop
		end if
		if(ieee_is_nan(f(1)).or..not.(ieee_is_finite(f(1))))then
			print*, "in fden"
			print*, "x, f=", x, f(1), yout, r0_cl/10**logr-x, r0_cl, logr
			call common_gx_ir%print("common_gx")
			stop
		end if
		!call gx%print("gx")
		
		!read(*,*)
	end subroutine
end subroutine


subroutine 	get_fmden(fden, fmden, mass)
    use com_sts_type
    implicit NONE
    type(s1d_type)::fden, fmden
    real(8) mass
    integer i
    
    do i=1, fden%nbin
        !fmden%fx(i)=fden%fx(i)+log10(mass)
        fmden%fx(i)=fden%fx(i)*mass
    end do	
end subroutine   
subroutine get_fmden_mb(mb,source)
    use model_basic
    use MPI_comu,only:rid
    use md_mass_bins
    implicit none
    type(mass_bins)::mb
    integer i,source
    logical allzero
    type(s1d_type)::fden
    
    do i=1, n_tot_comp
        !call so%fmden%init(ctl%log10rmin_factor, ctl%log10rmax_factor, dms%dstr_bins, coeff_sts_type_dc)
        !if(mb%dsp(i)%p%n>0)then
        allzero=all(mb%dsp(i)%p%barge_ir%fx.eq.0)
        !if(rid.eq.0.or.rid.eq.15)then
        !    print*, "rid, i, allzero=",rid, i,  allzero
        !end if
        if(.not.allzero)then
            select case (source)
            case(source_ana)
                fden=mb%dsp(i)%p%fden
            case(source_simu)
                fden=mb%dsp(i)%p%fden_simu
            case default
                print*, "error define source in fmden_mb"
                stop
            end select
            !print*, "mb%dsp(i)%p%n, mb%mc=", mb%dsp(i)%p%n, mb%mc
            call get_fMden(fden, mb%dsp(i)%p%fMden, mb%mc)
            !call  mb%dsp(i)%p%fden%print("mb%dsp(i)%p%fden")
        else
            mb%dsp(i)%p%fMden%fx=0
        end if
    end do		
end subroutine	
  
subroutine get_fma0(dm)
    use com_main_gw
    implicit none
    type(diffuse_mspec)::dm
    integer i,j
    !call dms%mb(1)%star%fden%print("fden before")
    !read(*,*)
    do i=1, dm%n

        call get_fma_mb(dm%mb(i))
        
        associate(mb=>dm%mb(i))
            mb%all%fma%fx=0
            mb%all%fma_simu%fx=0
            do j=1, n_tot_comp
                mb%all%fma%fx=mb%all%fma%fx+mb%dsp(j)%p%fma%fx
                mb%all%fma_simu%fx=mb%all%fma_simu%fx+mb%dsp(j)%p%fma_simu%fx
            end do
        end associate
    end do

    associate(all=>dm%all)
        do j=1, n_tot_comp
            all%dsp(j)%p%fma%fx=0
            all%dsp(j)%p%fma_simu%fx=0
        end do    
        all%all%fma%fx=0
        all%all%fma_simu%fx=0
        do i=1, dm%n
            do j=1, n_tot_comp
                all%dsp(j)%p%fma%fx=all%dsp(j)%p%fma%fx+dm%mb(i)%dsp(j)%p%fma%fx
                all%dsp(j)%p%fma_simu%fx=all%dsp(j)%p%fma_simu%fx+dm%mb(i)%dsp(j)%p%fma_simu%fx
            end do
        end do
        do j=1, n_tot_comp
            all%all%fma%fx=all%all%fma%fx+all%dsp(j)%p%fma%fx
            all%all%fma_simu%fx=all%all%fma_simu%fx+all%dsp(j)%p%fma_simu%fx
        end do
    end associate
end subroutine

subroutine get_dens0(dm,spp)
    use com_main_gw
    implicit none
    type(diffuse_mspec)::dm
    integer i,j,k,ier,flag_fden_out
    real(8) t1,t2
    type(star_pot_para)::spp
    !call dms%mb(1)%star%fden%print("fden before")
    !read(*,*)
    !if(rid.eq.ctl%ntasks-1)then
    !    call cpu_time(t1)
   ! end if

    do i=1, dm%n
        !call get_dens_mb(dm,spp, dm%mb(i), dm%n0, dm%v0,flag_fden_out )

        do j=1, n_tot_comp_sg
            ! if(ctl%bin_mass_simulation_particle_number_comp_tot(i,j)>ctl%min_sample_in_mass_bin)then
                call get_dens_so(dm, spp, dm%mb(i)%dsp(j)%p, flag_fden_out)            
                if(flag_fden_out.eq.1)then
                    if(rid.eq.0)then
                        print*, "get_dens_mbh:component number=", i
                    end if
                    print*, "n,rid=", ctl%bin_mass_simulation_particle_number_comp_tot(i,j),rid
                end if
            ! end if
        end do

        if(flag_fden_out.eq.1)then
            call mpi_BARRIER(mpi_comm_world,ier)
            if(rid.eq.0)then
                print*, "get_dens0, i, bintot_number=", i, ctl%bin_mass_simulation_particle_number_tot(i)
                stop
            end if
        end if
        call get_fmden_mb(dm%mb(i),ctl%source_fden)


        ! do j=1, ctl%ntasks
        !     if(rid.eq.j-1)then
        !         print*, "rid=",rid
        !         if(i.eq.7)then
        !            print*, "mass i=",i
        !            ! call dm%mb(j)%all%fden%print("all fden")
        !            ! print*, "all%n=",dm%mb(j)%all%n
        !            print*, "sample_logrmin,sample_logrmax=",sample_logrmin,sample_logrmax
        !             do k=1, n_tot_comp_sg
        !                 print*, 'type=',k
        !                 print*, "so%n=",dm%mb(i)%dsp(k)%p%n
        !                 print*, "rho_rmax,rmin=", dm%mb(i)%dsp(k)%p%spt_rho_rmax,dm%mb(i)%dsp(k)%p%spt_rho_rmin
        !                 if(dm%mb(i)%dsp(k)%p%n>0)then
        !                     call dm%mb(i)%dsp(k)%p%fden%print("fden")
        !                 end if
        !             end do
        !         endif
        !     end if
        !     call mpi_barrier(mpi_comm_world,ier)
        ! end do
        

        associate(mb=>dm%mb(i))
            mb%all%fden%fx=0; 
            mb%all%fden_simu%fx=0
            mb%all%fmden%fx=0; 
            mb%all%spt_rho_rmin=0
          !  mb%all%fma%fx=0
            do j=1, n_tot_comp
                mb%all%fden%fx= mb%all%fden%fx+ mb%dsp(j)%p%fden%fx 
                mb%all%fmden%fx= mb%all%fmden%fx+ mb%dsp(j)%p%fmden%fx 
                mb%all%fden_simu%fx= mb%all%fden_simu%fx+ mb%dsp(j)%p%fden_simu%fx 
                mb%all%spt_rho_rmin=mb%all%spt_rho_rmin+mb%dsp(j)%p%spt_rho_rmin
            end do

        end associate
    end do
    !stop

    associate(all=>dm%all)
        do j=1, n_tot_comp
            all%dsp(j)%p%fden%fx= 0
            all%dsp(j)%p%fden_simu%fx= 0
            all%dsp(j)%p%fmden%fx= 0
        end do    
        all%all%fden%fx= 0
        all%all%fden_simu%fx= 0
        all%all%fmden%fx=0
        do i=1, dm%n
            do j=1, n_tot_comp
                all%dsp(j)%p%fden%fx= all%dsp(j)%p%fden%fx+ dm%mb(i)%dsp(j)%p%fden%fx 
                all%dsp(j)%p%fden_simu%fx= all%dsp(j)%p%fden_simu%fx+ dm%mb(i)%dsp(j)%p%fden_simu%fx 
                all%dsp(j)%p%fmden%fx= all%dsp(j)%p%fmden%fx+ dm%mb(i)%dsp(j)%p%fmden%fx 
            end do
        end do
        do j=1, n_tot_comp
            all%all%fden%fx= all%all%fden%fx+ all%dsp(j)%p%fden%fx
            all%all%fden_simu%fx= all%all%fden_simu%fx+ all%dsp(j)%p%fden_simu%fx
            all%all%fmden%fx= all%all%fmden%fx+ all%dsp(j)%p%fmden%fx
        end do
    end associate


end subroutine
subroutine get_fna0(dm)
    use com_main_gw
    implicit none
    type(diffuse_mspec)::dm
    integer i,j,output_flag
    
    !call dms%mb(1)%star%fden%print("fden before")
    !read(*,*)
    do i=1, dm%n
        call get_fna_mb(dm%mb(i),output_flag)
        if(output_flag.ne.0)then
            print*, "error in mass bin idx=", i
            stop
        end if
        associate(mb=>dm%mb(i))
            mb%all%fna_simu%fx=0
            mb%all%fna%fx=0
            do j=1, n_tot_comp
                mb%all%fna%fx=mb%all%fna%fx+mb%dsp(j)%p%fna%fx
                mb%all%fna_simu%fx=mb%all%fna_simu%fx+mb%dsp(j)%p%fna_simu%fx
            end do
        end associate
    end do

    associate(all=>dm%all)
        do j=1, n_tot_comp
            all%dsp(j)%p%fna%fx=0
            all%dsp(j)%p%fna_simu%fx=0            
        end do    
        all%all%fna%fx=0
        all%all%fna_simu%fx=0

        do i=1, dm%n
            do j=1, n_tot_comp
                all%dsp(j)%p%fna%fx=all%dsp(j)%p%fna%fx+dm%mb(i)%dsp(j)%p%fna%fx
                all%dsp(j)%p%fna_simu%fx=all%dsp(j)%p%fna_simu%fx+dm%mb(i)%dsp(j)%p%fna_simu%fx
            end do
        end do
        do j=1, n_tot_comp
            all%all%fna%fx=all%all%fna%fx+all%dsp(j)%p%fna%fx
            all%all%fna_simu%fx=all%all%fna_simu%fx+all%dsp(j)%p%fna_simu%fx
        end do
    end associate
end subroutine

subroutine dms_so_get_nxj_from_nejw(so, jbtype)
    use md_stellar_object
    implicit none
    type(dms_stellar_object)::so
    integer ntasks, jbtype
    real(8) en(so%n), jm(so%n),we(so%n)
    if(so%n>0)then
        en(1:so%n)=so%nejw(1:so%n)%e
        jm(1:so%n)=so%nejw(1:so%n)%j
        we(1:so%n)=so%nejw(1:so%n)%w
        select case(jbtype)
        case(Jbin_type_lin)
            call so%nxj%get_stats_weight(en, jm, we, so%n)
        case(Jbin_type_log)
            !print*, so%nejw(1:10)%e
            !print*, so%nejw(1:10)%j
            !print*, so%nejw(1:10)%w
            !print*, so%n
            call so%nxj%get_stats_weight(en, log10(jm), we, so%n)
            !call so%nxj%print("nxj")
            !stop
        case(Jbin_type_sqr)
            call so%nxj%get_stats_weight(en, jm**2d0, we, so%n)
        case default
            print*, "dms_nxj_newj:error! define jbtype", jbtype
            stop
        end select
        
        !so%n_real=sum(so%nejw(1:so%n)%w)
    else
        !so%n_real=0
        return
    end if
end subroutine

subroutine dms_so_get_nxj_from_nejw_ir(so, jbtype)
    use md_stellar_object
    implicit none
    type(dms_stellar_object)::so
    integer ntasks, jbtype
    real(8) en(so%n), jm(so%n),we(so%n)
    integer i
    if(so%n>0)then
        en(1:so%n)=so%nejw(1:so%n)%e
        jm(1:so%n)=so%nejw(1:so%n)%j
        we(1:so%n)=so%nejw(1:so%n)%w
        select case(jbtype)
        case(Jbin_type_lin)
            !call so%nxj%get_stats_weight(en, jm, we, so%n)
            print*, "dms_so_get_nxj_from_nejw_ir, not finished"
            stop
        case(Jbin_type_log)
            

            !====================================
            call so%nxj_ir%get_hst(en, log10(jm), we, so%n)
            !====================================
            !so%nxj_ir%nbw
            
            
            !call get_s2d_ird_hst_weight_kernel(so%nxj_ir,en,log10(jm),we, so%n,1)

        case(Jbin_type_sqr)
            call so%nxj_ir%get_hst(en, jm**2d0, we, so%n)
        case default
            print*, "dms_nxj_newj:error! define jbtype", jbtype
            stop
        end select
        
        !so%n_real=sum(so%nejw(1:so%n)%w)
    else
        !so%n_real=0
        return
    end if
end subroutine

subroutine get_fMa_mb(mb)
    use md_mass_bins
    implicit none
    type(mass_bins)::mb
    integer i
    do i=1, n_tot_comp
        if(mb%dsp(i)%p%n>0)then
            call get_fMa(mb%dsp(i)%p%fNa, mb%dsp(i)%p%fMa, mb%mc)
            call get_fMa(mb%dsp(i)%p%fNa_simu, mb%dsp(i)%p%fMa_simu, mb%mc)
        else
            mb%dsp(i)%p%fMa%fx=0
            mb%dsp(i)%p%fMa_simu%fx=0
        end if
    end do		
end subroutine	

subroutine 	get_fma(fna, fma, mass)
    use md_mass_bins
    implicit NONE
    type(s1d_type)::fna, fma
    real(8) mass
    integer i
    
    do i=1, fna%nbin
        if(i>1)then
            if(fna%fx(i)<fna%fx(i-1)) fna%fx(i)=fna%fx(i-1)
        end if
        fma%fx(i)=fna%fx(i)*mass        
    end do	
end subroutine   

subroutine get_fna_mb(mb,output_flag)
    use md_dms
    use model_basic
    use MPI_comu,only:rid
    type(mass_bins)::mb
    integer i,output_flag
    do i=1, n_tot_comp
        !select case(ctl%source_fden)
        !case(source_ana)
            call get_fna(mb%dsp(i)%p%fden,mb%dsp(i)%p%fna,output_flag)
            if(output_flag.ne.0)then
                if(rid.eq.0)then
                    print*, "error in component idx: ", i
                    return
                end if
            end if
        !case(source_simu)
        !    call get_fna(mb%dsp(i)%p%fden_simu,mb%dsp(i)%p%fna_simu)
        !case default
        !    print*, "error! get_fna_mb", ctl%source_fden
        !    stop
        !end select
    end do
end subroutine


subroutine 	get_fna(fden, fna,output_flag)
	use constant
	use com_sts_type
	use my_intgl
    use MPI_comu,only:rid
	use, intrinsic :: ieee_arithmetic
    use model_basic
	implicit NONE
	type(s1d_type)::fden, fna
	integer i, idid,output_flag
	real(8) int_out
    output_flag=0
	do i=1, fna%nbin
        int_out=0
        call my_integral_acc(fden%xmin,fna%xb(i), int_out, fr_funcs_int_acc_a,fr_funcs_int_acc_r, fcn,idid)
        fna%fx(i)=4*pi*int_out*log(10d0)
        if(fna%fx(i)<0)then
            if(rid.eq.0)then
                print*, "warnning, fna%fx(i)<0, set fna%fx(i)=0", i, fna%xb(i), fna%fx(i),int_out
                call fden%print("fden")
            end if
            fna%fx(i)=0
            ! stop
            !output_flag=2
            !return
        end if

        if(idid<0)then
            print*, "idid,int_out=",idid,int_out
            call fden%print("fden")
            stop
        end if
	end do
	!fna%nsam=fden%nsam
    !call fden%print('fden')
    !call fna%print("fna")
    !read(*,*)
contains
	subroutine fcn(n, x, y, f, par, ipar)
		implicit none
		integer n, ipar(100)
		real(8) x, y(n), f(n), par(100), y_out
        !call get_rho_full_range()
        call fden%get_value_s(x,y_out)
        if(y_out>0)then
            !call get_value_at_x_fc(fden, x, y_out, 1)
            !if(y_out<0)y_out=0
            !f(1)=10**y_out*(10**x)**3
            f(1)=y_out*(10**x)**3
        else
            f(1)=0
        end if
        ! print*,"f=",f(1),x,y
		if(ieee_is_nan(f(1)).or..not.(ieee_is_finite(f(1))))then
			print*, "in fna"
			print*, "x, f=", x, f(1), y_out
            call fden%print("fden")
			stop
		end if

	end subroutine	
end subroutine

subroutine get_nxj(dm)
    use md_dms
    implicit none
    type(diffuse_mspec)::dm
    integer i, j, k,l
    !integer flag
    do i=1, dm%n
        associate(mb=>dm%mb(i))
            do j=1, n_tot_comp
                call dms_so_get_nxj_from_nejw(mb%dsp(j)%p,dm%jbin_type)
            end do
           
            !if(mb%star%n.ne.0.or.mb%sbh%n.ne.0.or.mb%bbh%n.ne.0.or.mb%bstar%n.ne.0)then
                !print*, "mc=", mb%mc
                !if(mb%star%n.ne.0) print*, "num_star=", mb%star%n, mb%star%n_real
                !if(mb%sbh%n.ne.0) print*, "num_sbh=", mb%sbh%n, mb%sbh%n_real
                !if(mb%bbh%n.ne.0) print*, "num_bbh=", mb%bbh%n, mb%bbh%n_real
                !if(mb%bstar%n.ne.0) print*, "num_bstar=", mb%bstar%n, mb%bstar%n_real
            !end if
            mb%all%nxj%nxyw=0
            do j=1, mb%all%nxj%nx
                do k=1, mb%all%nxj%ny
                    do l=1, n_tot_comp
                        mb%all%nxj%nxyw(j,k)=mb%all%nxj%nxyw(j,k)+mb%dsp(l)%p%nxj%nxyw(j,k)
                    end do
                end do 
            end do
        end associate
    end do
   ! print*,"get_nxj:"
   ! print*, dm%mb(1)%star%nxj%nxyw(:,3)
    
    !associate(star0=>dm%all%star%nxj, bstar0=>dm%all%bstar%nxj, &
    !    sbh0=>dm%all%sbh%nxj, bbh0=>dm%all%bbh%nxj, all0=>dm%all%all%nxj, &
    !	wd0=>dm%all%wd%nxj, ns0=>dm%all%ns%nxj, bd0=>dm%all%bd%nxj)
    do j=1, n_tot_comp
        dm%all%dsp(j)%p%nxj%nxyw=0
        !dm%all%dsp(j)%p%n=0
        !dm%all%dsp(j)%p%n_real=0
    end do
    dm%all%all%nxj%nxyw=0
    
    do i=1, dm%n
        !dm%mb(i)%all%n_real=0
        !dm%mb(i)%all%n=0
        do j=1, n_tot_comp
            dm%all%dsp(j)%p%nxj%nxyw=dm%all%dsp(j)%p%nxj%nxyw+dm%mb(i)%dsp(j)%p%nxj%nxyw
            !dm%mb(i)%all%n_real=dm%mb(i)%all%n_real+dm%mb(i)%dsp(j)%p%n_real
            !dm%mb(i)%all%n=dm%mb(i)%all%n+dm%mb(i)%dsp(j)%p%n
            !dm%all%dsp(j)%p%n=dm%all%dsp(j)%p%n+dm%mb(i)%dsp(j)%p%n
            !dm%all%dsp(j)%p%n_real=dm%all%dsp(j)%p%n_real+dm%mb(i)%dsp(j)%p%n_real
        end do
        dm%all%all%nxj%nxyw=dm%all%all%nxj%nxyw+dm%mb(i)%all%nxj%nxyw
        !dm%all%all%n=dm%all%all%n+dm%mb(i)%all%n
    end do    
   
end subroutine

subroutine get_nxj_ir(dm)
    use md_dms
    implicit none
    type(diffuse_mspec)::dm
    integer i, j, k,l
    do i=1, dm%n
        associate(mb=>dm%mb(i))

            do j=1, n_tot_comp
                call dms_so_get_nxj_from_nejw_ir(mb%dsp(j)%p,dm%jbin_type)
            end do

            mb%all%nxj_ir%nxyw=0
            do j=1, mb%all%nxj_ir%nx
                do k=1, mb%all%nxj_ir%ny
                    do l=1, n_tot_comp
                        mb%all%nxj_ir%nxyw(j,k)=mb%all%nxj_ir%nxyw(j,k)+mb%dsp(l)%p%nxj_ir%nxyw(j,k)
                    end do
                end do 
            end do
        end associate
    end do
   ! print*,"get_nxj_ir:"
   ! print*, dm%mb(1)%star%nxj_ir%nxyw(:,3)
    
    !associate(star0=>dm%all%star%nxj_ir, bstar0=>dm%all%bstar%nxj_ir, &
    !    sbh0=>dm%all%sbh%nxj_ir, bbh0=>dm%all%bbh%nxj_ir, all0=>dm%all%all%nxj_ir, &
    !	wd0=>dm%all%wd%nxj_ir, ns0=>dm%all%ns%nxj_ir, bd0=>dm%all%bd%nxj_ir)
    do j=1, n_tot_comp
        dm%all%dsp(j)%p%nxj_ir%nxyw=0
        !dm%all%dsp(j)%p%n=0
        !dm%all%dsp(j)%p%n_real=0
    end do
    dm%all%all%nxj_ir%nxyw=0
    
    do i=1, dm%n
        !dm%mb(i)%all%n_real=0
        !dm%mb(i)%all%n=0
        do j=1, n_tot_comp
            dm%all%dsp(j)%p%nxj_ir%nxyw=dm%all%dsp(j)%p%nxj_ir%nxyw+dm%mb(i)%dsp(j)%p%nxj_ir%nxyw
            !dm%mb(i)%all%n_real=dm%mb(i)%all%n_real+dm%mb(i)%dsp(j)%p%n_real
            !dm%mb(i)%all%n=dm%mb(i)%all%n+dm%mb(i)%dsp(j)%p%n
            !dm%all%dsp(j)%p%n=dm%all%dsp(j)%p%n+dm%mb(i)%dsp(j)%p%n
            !dm%all%dsp(j)%p%n_real=dm%all%dsp(j)%p%n_real+dm%mb(i)%dsp(j)%p%n_real
        end do
        dm%all%all%nxj_ir%nxyw=dm%all%all%nxj_ir%nxyw+dm%mb(i)%all%nxj_ir%nxyw
        !dm%all%all%n=dm%all%all%n+dm%mb(i)%all%n
    end do    
   
end subroutine


