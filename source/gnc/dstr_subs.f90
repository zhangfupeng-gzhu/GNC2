 

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
    real(8) rc,rmax, jc_dmless,ex,p_EJ_dmless, jm
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
    end do  
    do i=1, fden%nbin 
        if(n_per_bin(i)*10**(fden%xb(i))>ctl%den_bin_cri)then
            r1=10**(fden%xb(i)-fden%xstep/2d0)
            r2=10**(fden%xb(i)+fden%xstep/2d0)
            fden%fx(i)=fden%fx(i)/(4*pi/3d0*(r2**3-r1**3))/m0_cl
        else
            fden%fx(i)=0
        endif
        !read(*,*)
    end do 
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
             
        end do
        jm(1:so%n)=so%nejw(1:so%n)%j
        
         
		call so%nx_ir%get_hst(en, we, so%n)
         
        so%nx_ir%fxw=so%nx_ir%fxw/(m0_cl)
        so%nx_ir%fx=so%nx_ir%fx/(m0_cl)
        so%nx_ir%nbw=so%nx_ir%nbw/(m0_cl)
         
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
    allzero=all(so%barge_ir%fx.eq.0)
    if(.not.allzero)then 
            select case(ctl%fden_ana_est_method)
            case(fden_ana_est_method_1d_iso)
                !print*, "spt_rho_rmax=",so%spt_rho_rmax
                call get_fden_ir(dm,spp,so%barge_ir, so%fden, so%spt_rho_rmax, &
                 so%spt_rho_rmin,flag_fden_out)
 
            case(fden_ana_est_method_2d) 

                call get_fden_ir_2d(dm,spp,so%gxj_ir, so%fden,  so%spt_rho_rmin) 
            end select
 
        if(ctl%source_fden.eq.source_simu)then
            call get_fden_sample_particle(dm,so,so%fden_simu)
        end if 
    else
        so%fden%fx=0d0
        so%fden_simu%fx=0d0
    end if
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
            if(x1<=ex.and.ex<=x2) then
                n_per_bin(i)=n_per_bin(i)+1
                barge%fx(i)=barge%fx(i)+1d0/ex/jc_xy**2/pd_xy*w
            end if                
        end do
    end do
    do i=1, barge%nbin 
		barge%fx(i)=barge%fx(i)/(barge%xstep*log(10d0))*pi**(-0.5d0)*2**(-0.5d0)/r0_cl**3/ctl%n0
    end do 
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
            end if               
        end do
    end do
    do i=1, barge%nbin
        
        x1=10**(barge%xb(i)-barge%xsteps(i)/2d0)
        x2=10**(barge%xb(i)+barge%xsteps(i)/2d0) 
            barge%fx(i)=barge%fx(i)/(x2-x1)*pi**(-0.5d0)*2**(-0.5d0)/r0_cl**3/ctl%n0
        !else
        !    barge%fx(i)=0d0
        !end if
    end do
 
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
        call get_fden_ird_one(fden%xb(i),fden%fx(i),common_gx_ir,dm,spp,rho_rmax,flag_out)
        if(flag_out.eq.1)then
            if(rid.eq.0)then
                print*, "in bin i,xb=",i, fden%xb(i)
            end if
            flag_fden_out=1
        end if 
	end do 
    call get_fden_ird_one(sample_logrmin,rho_rmin,common_gx_ir,dm,spp,rho_rmax,flag_out)
    if(flag_out.eq.1)then
        if(rid.eq.0)then
            print*, "in sample_logrmin=",sample_logrmin
        end if
        flag_fden_out=1
    end if 
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
     
    select case(ctl%ebin_type)
    case(ebin_type_log)
        exmin=10**(common_gx_ir%xmin) 
    case default
        print*, "get_fden_ird_one:stop"
        stop
    end select

    call get_phi_star_full_range(spp, logr, exmax)
    !call get_dms_starpt_one(fden%xb(i),xmax,dms%all%all%fmden)
    select case(ctl%ebin_type)
    case(ebin_type_log)
        phi_tot=min(10**exmax+spp%mbh_dmless/10**logr,10**common_gx_ir%xmax) 
    end select
     
    if(exmin+1d-12<phi_tot)then
       
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
            flag_out=1 
        end if
    else
        int_out=0d0
    end if
    !call my_integral_acc(xmin,phi_tot,int_out, fcn,idid)
    
    fx=max(2/pi**0.5*int_out,0d0)+rho_rmax !*n0 ! in unit of r0_cl^{-3}
     
contains
	subroutine fcn(n, x, y, f, par, ipar)
		implicit none
		integer n, ipar(100)
		real(8) x, y(n), f(n), par(100), yout,xin

		select case(ctl%ebin_type)
        case(ebin_type_log)
            xin=max(log10(x),common_gx_ir%xmin ) 
        end select
        !print*, "1"
        !call get_gx_full_range_ir_log(gx_ir_tmp,xin,yout)
        call get_gx_full_range_ir(common_gx_ir,xin,yout)
        !=====================
         
        yout=max(yout,0d0)
		f(1)=yout*sqrt(max(phi_tot-x,0d0))
        if(debug.eq.1)then
            write(*,fmt="(A20,4E30.20)")  "x, yout, f=", x,xin, yout, f(1)
        endif 
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
        case(Jbin_type_log) 
            call so%nxj%get_stats_weight(en, log10(jm), we, so%n)
        case default
            print*, "dms_nxj_newj:error! define jbtype", jbtype
            stop
        end select 
    else 
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
        case(Jbin_type_log)
            !====================================
            call so%nxj_ir%get_hst(en, log10(jm), we, so%n)
            !==================================== 
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
contains
	subroutine fcn(n, x, y, f, par, ipar)
		implicit none
		integer n, ipar(100)
		real(8) x, y(n), f(n), par(100), y_out
        !call get_rho_full_range()
        call fden%get_value_s(x,y_out)
        if(y_out>0)then 
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
    do j=1, n_tot_comp
        dm%all%dsp(j)%p%nxj_ir%nxyw=0 
    end do
    dm%all%all%nxj_ir%nxyw=0
    
    do i=1, dm%n 
        do j=1, n_tot_comp
            dm%all%dsp(j)%p%nxj_ir%nxyw=dm%all%dsp(j)%p%nxj_ir%nxyw+dm%mb(i)%dsp(j)%p%nxj_ir%nxyw 
        end do
        dm%all%all%nxj_ir%nxyw=dm%all%all%nxj_ir%nxyw+dm%mb(i)%all%nxj_ir%nxyw 
    end do    
   
end subroutine


