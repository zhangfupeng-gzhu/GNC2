
subroutine get_fxj0(dm)
	use com_main_gw
	implicit none
	type(diffuse_mspec)::dm
	integer i

	do i=1, dm%n
		call get_fxj(dm,dm%mb(i),dm%jbin_type)
	end do
   ! print*, "get_fxj0:"
   ! print*, dm%mb(1)%star%nxj%nxyw(:,3)
	call get_fxj(dm, dm%all,dm%jbin_type)
end subroutine
 

subroutine get_fxj0_ir(dm)
	use com_main_gw
	implicit none
	type(diffuse_mspec)::dm
	integer i,j,k,l

	do i=1, dm%n
		call get_fxj_ir(dm,dm%mb(i),dm%jbin_type)
        associate(mb=>dm%mb(i))
            mb%all%gxj_ir%fxy=0
            do j=1, mb%all%gxj_ir%nx
                do k=1, mb%all%gxj_ir%ny
                    do l=1, n_tot_comp
                        mb%all%gxj_ir%fxy(j,k)=mb%all%gxj_ir%fxy(j,k)+mb%dsp(l)%p%gxj_ir%fxy(j,k)
                    end do
                end do 
            end do
        end associate
	end do
    do j=1, n_tot_comp
        dm%all%dsp(j)%p%gxj_ir%fxy=0
    end do
    dm%all%all%gxj_ir%fxy=0
    
    do i=1, dm%n
        
        do j=1, n_tot_comp
            dm%all%dsp(j)%p%gxj_ir%fxy=dm%all%dsp(j)%p%gxj_ir%fxy+dm%mb(i)%dsp(j)%p%gxj_ir%fxy
        end do
        dm%all%all%gxj_ir%fxy=dm%all%all%gxj_ir%fxy+dm%mb(i)%all%gxj_ir%fxy
    end do    
   ! print*, "get_fxj0:"
   ! print*, dm%mb(1)%star%nxj%nxyw(:,3)
	!call get_fxj(dm, dm%all,dm%jbin_type)
end subroutine
subroutine get_fxj_ir(dm,mb, jbtype)
	use md_dms
	implicit none
	type(mass_bins)::mb
    type(diffuse_mspec)::dm
	integer i, j, jbtype
	real(8) jm ,x

    do i=1, n_tot_comp
        call dms_so_get_fxj_spt_ir(mb%dsp(i)%p, mb%n0, dm%pd,dm%jc,  mb%r0_cl,jbtype)
    end do
    !mb%all%gxj_ir%
    !call dms_so_get_fxj_spt_ir(mb%all, mb%n0, dm%pd,dm%jc, mb%r0_cl,jbtype)
end subroutine
subroutine get_fxj(dm,mb, jbtype)
	use md_dms
	implicit none
	type(mass_bins)::mb
    type(diffuse_mspec)::dm
	integer i, j, jbtype
	real(8) jm ,x


    call dms_so_get_fxj_spt(mb%all, mb%n0, dm%pd,dm%jc, mb%r0_cl,jbtype)
    do i=1, n_tot_comp
        call dms_so_get_fxj_spt(mb%dsp(i)%p, mb%n0, dm%pd,dm%jc,  mb%r0_cl,jbtype)
    end do

end subroutine



subroutine get_ge(dm, sms_arr_single)
    use com_main_gw
    implicit none
    type(diffuse_mspec)::dm
    type(particle_samples_arr_type)::sms_arr_single
    type(particle_samples_arr_type)::sms_arr_stellar_sg(n_tot_comp_sg)     
    logical::norm_in
    if(rid.eq.0)then
        print*, "get ge"
    end if
    call separate_to_species(sms_arr_single, sms_arr_stellar_sg) 
    call init_diffuse_mspec_gxtables(dm)
     
    call conv_dms_nejw(dm, sms_arr_stellar_sg)        
    call get_dms_particle_numbers(dm)
    call update_asymptot(dm)
      
    if(ctl%include_stellar_evolution.ge.1)then
        
        call get_mass_distributions(dm,sms_arr_single)
        call get_rad_distributions(dm, sms_arr_single)
        call collect_dms_fmr(dm)
    end if
    call ignore_bins_with_few_samples(dm)

    if(rid.eq.0)then
        print*, "finished get_ge"
    end if
end subroutine
subroutine get_dms_particle_numbers(dm)
    use com_main_gw
    implicit none
    type(diffuse_mspec)::dm
    integer i,j 
    
    do i=1, dm%n
        dm%mb(i)%all%n_real=0
        dm%mb(i)%all%n=0
        do j=1, n_tot_comp
            dm%mb(i)%all%n_real=dm%mb(i)%all%n_real+dm%mb(i)%dsp(j)%p%n_real
            dm%mb(i)%all%n=dm%mb(i)%all%n+dm%mb(i)%dsp(j)%p%n 
        end do
        dm%all%all%n=dm%all%all%n+dm%mb(i)%all%n
    end do     

end subroutine
subroutine update_asymptot(dm)
    use com_main_gw
    implicit none
    type(diffuse_mspec)::dm
    real(8) np(dm%n, n_tot_comp_sg), nptot_mbin(dm%n),nptot
    real(8) mp(dm%n, n_tot_comp_sg), mptot_mbin(dm%n),mptot
    integer nptot_simu(dm%n), ns(dm%n,n_tot_comp_sg)
    integer i,j
    np=0
    do i=1, dm%n
        do j=1, n_tot_comp_sg
            np(i,j)=dm%mb(i)%dsp(j)%p%n_real !call dm%mb(i)%dsp(j)%p%get_n()
        end do
        nptot_mbin(i)=sum(np(i,1:n_tot_comp_sg))
    end do
    nptot=sum(nptot_mbin(1:dm%n))
    call collection_and_avg_s2d(np,dm%n,n_tot_comp_sg)
    call collection_and_avg_fx(nptot_mbin,dm%n)
    call collection_and_avg_real(nptot)
    ctl%num_particle_tot=nptot
    
    ctl%asymptot_now(1,1:dms%n)=nptot_mbin(1:dm%n)/nptot
    do i=1, n_tot_comp_sg
        do j=1, dm%n
            if(nptot_mbin(j).ne.0)then
                ctl%asymptot_now(i+1, j)=np(j,i)/nptot_mbin(j)
            else
                ctl%asymptot_now(i+1, j)=0
            end if
        end do
    end do
    
    do i=1, n_tot_comp_sg
        ctl%type_num_particle_tot(i)=sum(np(1:dm%n,i))
    end do

    call update_particle_existence()
    if(rid.eq.0)then
        write(*,fmt="(A10,20I4)") "exists(sg):", ctl%exist_stellar_type(1:n_tot_comp_sg)
    end if
    mp=0
    do i=1, dm%n
        do j=1, n_tot_comp_sg
            call dm%mb(i)%dsp(j)%p%get_mtot(mp(i,j))
        end do
        mptot_mbin(i)=sum(mp(i,1:n_tot_comp_sg))
    end do
    mptot=sum(mptot_mbin(1:dm%n))
    
    call collection_and_avg_s2d(mp,dm%n,n_tot_comp_sg)
    call collection_and_avg_fx(mptot_mbin,dm%n)
    call collection_and_avg_real(mptot)
    ctl%mass_particle_tot=mptot
    
    ctl%bin_fracmass_now(1:dms%n)=mptot_mbin(1:dm%n)/mptot
    do i=1, n_tot_comp_sg
        ctl%type_mass_particle_tot(i)=sum(mp(1:dm%n,i))
    end do    

    do i=1, dm%n
        do j=1, n_tot_comp_sg
            ctl%bin_mass_simulation_particle_number_comp_tot(i,j)=dm%mb(i)%dsp(j)%p%n
        end do
    end do
    ns=ctl%bin_mass_simulation_particle_number_comp_tot(1:dm%n,1:n_tot_comp_sg)
    call collection_fxy_int(ns(1:dm%n,1:n_tot_comp_sg),dm%n, n_tot_comp_sg)
    ctl%bin_mass_simulation_particle_number_comp_tot(1:dm%n,1:n_tot_comp_sg)=ns(1:dm%n,1:n_tot_comp_sg)
    do i=1, dm%n
        ctl%bin_mass_simulation_particle_number_tot(i)=sum(ctl%bin_mass_simulation_particle_number_comp_tot(i,1:n_tot_comp_sg))
    end do
    
end subroutine

subroutine get_mass_distributions(dm,sams)
    use com_main_gw
    implicit none
    type(diffuse_mspec)::dm
    type(particle_samples_arr_type)::sams
    integer i
    real(8) m1, m2 
    m1=minval(sams%sp(:)%m)
    m2=maxval(sams%sp(:)%m)
    !print*, "m1,m2=",m1,m2,rid
    call collection_and_get_min_real(m1)
    call collection_and_get_max_real(m2)
    
    !print*, "m1,m2=",m1,m2
    do i=1, n_tot_comp
        !print*, "i,so%n=", i, dm%all%dsp(i)%p%n
        call get_mass_distr_one(dm%all%dsp(i)%p, log10(m1), log10(m2))
    end do
    dm%all%all%fm_dstr=dm%all%star%fm_dstr
    dm%all%all%fm_dstr%fx=0
    dm%all%all%fm_dstr%fxw=0
    do i=1, n_tot_comp
        dm%all%all%fm_dstr%fxw=dm%all%all%fm_dstr%fxw+dm%all%dsp(i)%p%fm_dstr%fxw
    end do
    !print*, "---"

end subroutine
subroutine get_mass_distr_one(so, logm1, logm2)
    use md_stellar_object
    implicit none
    type(dms_stellar_object)::so
    real(8) logm1, logm2
    real(8) m(so%n),w(so%n)
    integer n

    call so%fm_dstr%init(logm1, logm2, 50, use_weight=.true.)
    call so%fm_dstr%set_range()
    if(so%n>0)then
        
        if(size(so%nejw(:)).ne.so%n)then
            print*, "warnning! size(so%nejw(:)).ne.so%n=", size(so%nejw(:)), so%n, logm1, logm2
            print*, "m=",so%nejw(1:min(10,size(so%nejw(:))))%m
        end if
        n=min(so%n,size(so%nejw))
        m(1:n)=log10(so%nejw(1:n)%m)
        w(1:n)=so%nejw(1:n)%w
        call so%fm_dstr%get_hst(m(1:n),w(1:n),n)
    end if
end subroutine


subroutine get_rad_distributions(dm,sams)
    use com_main_gw
    implicit none
    type(diffuse_mspec)::dm
    type(particle_samples_arr_type)::sams
    integer i
    real(8) r1, r2 
    r1=minval(sams%sp(:)%byot%ms%radius)
    r2=maxval(sams%sp(:)%byot%ms%radius)
    !print*, "m1,m2=",m1,m2,rid
    call collection_and_get_min_real(r1)
    call collection_and_get_max_real(r2)
    
    !print*, "m1,m2=",m1,m2
    do i=1, n_tot_comp
        !print*, "i,so%n=", i, dm%all%dsp(i)%p%n
        call get_rad_distr_one(dm%all%dsp(i)%p, log10(r1), log10(r2))
    end do
    dm%all%all%fr_dstr=dm%all%star%fr_dstr
    dm%all%all%fr_dstr%fx=0
    dm%all%all%fr_dstr%fxw=0
    do i=1, n_tot_comp
        dm%all%all%fr_dstr%fxw=dm%all%all%fr_dstr%fxw+dm%all%dsp(i)%p%fr_dstr%fxw
    end do
    !print*, "---"

end subroutine
subroutine get_rad_distr_one(so, logr1, logr2)
    use md_stellar_object
    implicit none
    type(dms_stellar_object)::so
    real(8) logr1, logr2
    real(8) r(so%n),w(so%n)

    call so%fr_dstr%init(logr1, logr2, 50, use_weight=.true.)
    call so%fr_dstr%set_range()
    r(1:so%n)=log10(so%nejw(1:so%n)%rad)
    w(1:so%n)=so%nejw(1:so%n)%w
    call so%fr_dstr%get_hst(r(1:so%n),w(1:so%n),so%n)
    
end subroutine

subroutine ignore_bins_with_few_samples(dm)
    use com_main_gw
     
    implicit none
    type(diffuse_mspec)::dm
    integer i,j
    real(8) bin_n_simu
    do i=1, dm%n
        bin_n_simu=ctl%bin_mass_simulation_particle_number_tot(i)
        if(bin_n_simu.lt.ctl%min_sample_in_mass_bin.and.bin_n_simu.gt.0)then
            if(rid.eq.0)then
                print*, "bin=", i, " number=", bin_n_simu, " too small, set to zero"
            end if
            do j=1, n_tot_comp
                dm%mb(i)%dsp(j)%p%n=0
                dm%mb(i)%dsp(j)%p%n_real=0
            end do
            dm%mb(i)%all%n=0
            dm%mb(i)%all%n_real=0
        end if
    end do
    call update_dms_numbers(dm)
end subroutine
subroutine set_jm_bound(jm)
    use com_main_gw,only:jmin_value,jmax_value,ctl
    implicit none
    integer niter
    real(8) jm
    niter=0
100   if(jm>jmax_value)then
        !print*, "jm=",jm, "->", jmax_value
        if(ctl%chattery.ge.3)then
            print*, "set_jm_bound: jm=", jm, " -> ", jmax_value
        end if
        jm=jmax_value!*2-jm 
        
    elseif(jm<jmin_value)then
        if(ctl%chattery.ge.3)then
            print*, "set_jm_bound: jm=", jm, " -> ", jmin_value*2-jm
        end if 
        jm=jmin_value*2-jm 
    end if
    if(jm<jmin_value.or.jm>jmax_value) then
        niter=niter+1
        if(niter>1000)then
            print*, "error! jm=",jm, jmin_value, jmax_value
            stop
        end if
        goto 100
    end if
end subroutine
subroutine get_fxr_ir(dm,mb)
    use com_main_gw
    implicit none
    type(diffuse_mspec)::dm
    type(mass_bins)::mb
    integer i
    do i=1, n_tot_comp
        call get_fxr_ir_one(mb%dsp(i)%p)
    end do
    
end subroutine
subroutine get_fxr_ir_one(so)
    use com_main_gw
    implicit none
    type(dms_stellar_object)::so
    call so%gxr_ir%init(ctl%diff_coeff_bins,ctl%dstr_bins_r,so%gxj_ir%xmin,so%gxj_ir%xmax,&
    dms%logrmin,dms%logrmax,sts_type_dstr)
    so%gxr_ir%xcenter=so%gxj_ir%xcenter
    so%gxr_ir%ycenter=so%fden%xb

end subroutine
subroutine get_cluster_density(dm)
    use com_main_gw
    implicit none
    type(s1d_type)::fphi
    real(8) cri
    type(diffuse_mspec)::dm
    real(8) logrmin,logrmax
    real(8) t1, t2
   ! call cpu_time(t1)

    call get_sample_rrange(logrmin,logrmax)
    if(rid.eq.0.and.ctl%chattery.ge.2)then
        print*, "logrmin,logrmax=",logrmin,logrmax
    end if
    call get_rho_rmax()
    call init_stellar_obj_rtables(dm) 
    if(ctl%chattery.ge.2.and.rid.eq.0)then
        print*, "start get_dens0"
    end if
    call get_dens0(dm,spp_new)
    if(ctl%chattery.ge.2.and.rid.eq.0)then
        call dm%all%all%fmden%print("fmden")
        print*, "start get_fna0"        
    end if
    call get_fna0(dm)
    if(ctl%chattery.ge.2.and.rid.eq.0)then
        print*, "start get_fma0"
    end if
    call get_fma0(dm)
    !print*, "2"
    call get_slope0(dm)
   ! print*, "10"

end subroutine

subroutine mb_get_gxjcr_mpi(mb)
    use com_main_gw
    implicit none
	integer n,i,j, ierr
	type(coeff_type)::cej
	real(8) kappa, sigma32, n0
	type(mass_bins)::mb
	integer nbg, ned, bks

	n=mb%all%gxjcr%nx
	bks=int(ctl%diff_coeff_bins/ctl%ntasks)
	nbg=bks*rid+1
	ned=bks*(rid+1)
	!print*, "n=",n
	!print*,"nbg, ned=", nbg, ned
	do i=nbg, ned
		do j=1, n
			!print*,"x,y=", mb%all%gxjcr%xcenter(i),mb%all%gxjcr%ycenter(j)
			call get_gxjcr(mb%all%gxj, mb%all%gxjcr%fxy(i,j), mb%all%gxjcr%xcenter(i), &
				mb%all%gxjcr%ycenter(j))
			!print*, "gxjcr=",mb%all%gxjcr%fxy(i,j)
		end do
		!print*, "gxjcr=",mb%all%gxjcr%fxy(i,:)
	end do

	!kappa=(4*pi*mb%mc)**2*log(mb%mtot/mb%mc)!*factor
    call mpi_BARRIER(mpi_comm_world, ierr)
    call collect_data_mpi_y(mb%all%gxjcr%fxy, n, nbg, ned, bks, ctl%ntasks)
	!print*, "finished"
end subroutine

subroutine dm_get_gxjcr_mpi(dm)
    use com_main_gw
    implicit none
    type(diffuse_mspec):: dm
    integer i 
    do i=1, dm%n 
        call mb_get_gxjcr_mpi(dm%mb(i))
    end do 

end subroutine 
subroutine set_dm_init(dm)
    use com_main_gw
    implicit none
    type(diffuse_mspec)::dm
    logical e_iregular

    call dm%init(ctl%m_bins) 
	e_iregular=.true. 
    call dm%set_diffuse_mspec(ctl%diff_coeff_bins,ctl%dstr_bins_r,ctl%dstr_bins_e, &
        ctl%log10rmin_factor, ctl%log10rmax_factor, &
        jmin_value ,jmax_value, m0_cl,  ctl%v0, ctl%n0, r0_cl,ctl%x_boundary, &
        ctl%idx_ref, ctl%ebin_type,ctl%jbin_type,e_iregular)

    !print*, "set_dm_init:", dm%dstr_bins_e,ctl%dstr_bins_e
end subroutine

   
subroutine get_n_from_particlex(mstar,xstar, n_star,&
    m1,  m2, idx, nsam)
    use model_basic
    implicit none
    integer n_star
    real(8):: mstar(n_star), xstar(n_star)
    integer i, nsam, idx(n_star)
    real(8) jmax, m1,m2
    nsam=0
    if(n_star.eq.0) return
    !print*, "emin_factor, emax_factor=",emin_factor,emax_factor
    do i=1, n_star
        !print*, ps_arr%sp(i)%m, m1, m2
        if(mstar(i).ge.m1.and.mstar(i).le.m2)then
            nsam=nsam+1
            idx(nsam)=i
        end if
    end do
    !print*, "nsam, n_star=",nsam, n_star
end subroutine
subroutine get_ejw_from_particle(estar, jstar, wstar,rpstar,rastar,pdstar,jcstar, &
    mstar, radstar,n_star,  m1,  m2, mbh_in, v0, xb, nejw, nsam, nwsam)
    use md_stellar_object
    use model_basic,only:r0_cl,ctl
    use md_coeff,only:ebin_type_log
    implicit none
    integer n_star
    real(8):: estar(n_star),jstar(n_star),wstar(n_star),mstar(n_star), xstar(n_star)
    real(8) rpstar(n_star),rastar(n_star),pdstar(n_star), jcstar(n_star),radstar(n_star)
    integer i, nsam, idx(n_star)
    type(nejw_type),allocatable::nejw(:)
    real(8) jmax, m1,m2, v0, mbh_in,xb, nwsam
    !integer jb_type

    xstar=abs(estar(1:n_star))/v0**2
    !print*, "v0=", v0
    !print*, "xstar=", xstar(1:10)
    call get_n_from_particlex(mstar,xstar, n_star, m1, m2, idx, nsam)
    !call get_n_from_particle(mstar, n_star, m1, m2, idx, nsam)
    if(allocated(nejw)) deallocate(nejw)
    allocate(nejw(nsam))
    nwsam=0
    do i=1, nsam
        select case(ctl%ebin_type)
        case(ebin_type_log)
            nejw(i)%e=log10(xstar(idx(i))) 
        end select
        if(ieee_is_nan(nejw(i)%e))then
            print*, "nejw=NaN:", estar(idx(i)), v0
            read(*,*)
        end if
 
        nejw(i)%j=jstar(idx(i))!/jmax 
        nejw(i)%rp=rpstar(idx(i))/r0_cl
        nejw(i)%ra=rastar(idx(i))/r0_cl
        nejw(i)%jc=jcstar(idx(i))/(v0*r0_cl)
        nejw(i)%pd=pdstar(idx(i))/(r0_cl/v0)
        nejw(i)%m=mstar(idx(i))
        nejw(i)%w=wstar(idx(i))
        nejw(i)%rad=radstar(idx(i))
        nejw(i)%idx=idx(i)
        if(isnan(wstar(idx(i))).or.wstar(idx(i))>1d90)then
            print*, "isnan wstar=",wstar(idx(i)),idx(i)
            stop
        end if
        nwsam=nwsam+wstar(idx(i))
    end do
     
end subroutine    
subroutine conv_dms_nejw_obj(dm,en,jm,m,rp,ra,pd,jc,rad,w_real,nobj, obj_type)
    use com_main_gw
    implicit none
    type(diffuse_mspec)::dm
    integer nobj
    real(8) en(nobj),jm(nobj), w_real(nobj), wobj(nobj), m(nobj)
    real(8) rp(nobj), ra(nobj), pd(nobj), jc(nobj), rad(nobj)
    integer i, obj_type,dsp_idx
    real(8) m1, m2
    !integer,parameter::obj_type_star=1, obj_type_sbh=2, obj_type_bstar=3
	!integer,parameter::obj_type_bbh=4, obj_type_wd=5, obj_type_ns=6, obj_type_bd=7
    interface
        subroutine get_ejw_from_particle(estar, jstar, wstar,rpstar,rastar,pdstar,jcstar, &
            mstar, radstar, n_star,  m1,  m2, mbh, v0, xb, nejw, nsam,nwsam)
            use md_stellar_object
            use model_basic,only:r0_cl
            implicit none
            integer n_star
            real(8):: estar(n_star),jstar(n_star),wstar(n_star),mstar(n_star), xstar(n_star)
            real(8) rpstar(n_star),rastar(n_star),pdstar(n_star), jcstar(n_star),radstar(n_star)
            integer nsam, idx(n_star)
            type(nejw_type),allocatable::nejw(:)
            real(8) m1,m2, v0, mbh,xb, nwsam
    end subroutine
    end interface
    call get_dsp_idx_by_type(obj_type,dsp_idx)
    if(nobj>0)then
        wobj(1:nobj)=w_real(1:nobj)!/dble(ctl%ntasks)
	else
        do i=1, dm%n
            dm%mb(i)%dsp(dsp_idx)%p%n=0;  dm%mb(i)%dsp(dsp_idx)%p%n_real=0
        end do
		return
    end if
    call get_ejw_from_particle(en(1:nobj),jm(1:nobj),wobj(1:nobj),rp(1:nobj),&
        ra(1:nobj),pd(1:nobj), jc(1:nobj), m(1:nobj), rad(1:nobj), &
        nobj, minval(ctl%bin_mass_m1(1:ctl%m_bins)), maxval(ctl%bin_mass_m2(1:ctl%m_bins)), dm%mtot, dm%v0, dm%x_boundary, &
        dm%all%dsp(dsp_idx)%p%nejw,dm%all%dsp(dsp_idx)%p%n,dm%all%dsp(dsp_idx)%p%n_real) 
    do i=1, dm%n
        associate(mb=>dm%mb(i))
            m1=mb%m1; m2=mb%m2           
            
            call get_ejw_from_particle(en(1:nobj),&
            jm(1:nobj),wobj(1:nobj),rp(1:nobj),&
            ra(1:nobj),pd(1:nobj), jc(1:nobj), m(1:nobj), rad(1:nobj), &
            nobj, m1, m2, dm%mtot, dm%v0, dm%x_boundary, mb%dsp(dsp_idx)%p%nejw,mb%dsp(dsp_idx)%p%n, &
            mb%dsp(dsp_idx)%p%n_real) 
        end associate      
    end do
    
end subroutine


subroutine conv_dms_newj_obj_one(dm, bk, i_obj)
    use com_main_gw
    implicit none
    type(particle_samples_arr_type)::bk
    type(diffuse_mspec)::dm
    integer i_obj,i
    real(8) en(bk%n), jm(bk%n), mass(bk%n), weight_real(bk%n)
    real(8) rp(bk%n), ra(bk%n), pd(bk%n), jc(bk%n), rad(bk%n)
    
    en=bk%sp(1:bk%n)%en
    jm=bk%sp(1:bk%n)%jm
    mass=bk%sp(1:bk%n)%m
    rp=bk%sp(1:bk%n)%rp
    ra=bk%sp(1:bk%n)%ra
    pd=bk%sp(1:bk%n)%period
    jc=bk%sp(1:bk%n)%jc
    rad=bk%sp(1:bk%n)%byot%ms%radius
    weight_real=bk%sp(1:bk%n)%weight_real 
    call conv_dms_nejw_obj(dm,en, jm , mass, rp,ra,pd,jc, rad, weight_real, bk%n, i_obj)
end subroutine

subroutine conv_dms_nejw(dm, bkstellar_sg)
    use com_main_gw
    implicit none
    type(diffuse_mspec)::dm
    type(particle_samples_arr_type)::bkstellar_sg(n_tot_comp_sg)
    integer i,j,obj_type
 
    do i=1, n_tot_comp_sg
        call get_type_from_ctl_obidx_sg(i,obj_type) 
        call conv_dms_newj_obj_one(dm, bkstellar_sg(i),obj_type) 
    end do
    call update_dms_numbers(dm)
end subroutine
subroutine update_dms_numbers(dm)
    use com_main_gw
    implicit none
    type(diffuse_mspec)::dm
    integer i,j
    do i=1, dm%n
        associate(mb=>dm%mb(i))
            mb%all%n=0
            mb%all%n_real=0
            do j=1, n_tot_comp
                mb%all%n=mb%all%n+mb%dsp(j)%p%n
                mb%all%n_real=mb%all%n_real+mb%dsp(j)%p%n_real
                
                if(mb%dsp(j)%p%n>0)then
                    if(allocated(mb%dsp(j)%p%nejw))then
                        if(size(mb%dsp(j)%p%nejw).ne.mb%dsp(j)%p%n)then
                            print*, "i,j,rid=",i,j,rid
                            print*, "size(newj)=",size(mb%dsp(j)%p%nejw)
                            print*, "error in update dms numbers:"
                            print*, "n=", mb%dsp(j)%p%n
                            print*, "n != size(nejw)"
                            stop
                        endif
                    else
                        print*, "i,j,rid=",i,j,rid
                            !print*, "size(newj)=",size(mb%dsp(j)%p%nejw)
                            print*, "n=", mb%dsp(j)%p%n
                            print*, "error in update dms numbers:"
                        print*, "array not allocated"
                        stop
                    end if
                  
                end if
            end do 
        end associate            
    end do
    associate(all=>dm%all)
         

        do i=1, n_tot_comp
            all%all%n=0
            all%all%n_real=0
            do j=1,dm%n
                all%all%n=all%all%n+dm%mb(j)%all%n
                all%all%n_real=all%all%n_real+dm%mb(j)%all%n_real
            end do
        end do
    end associate
end subroutine
subroutine separate_to_species(bks, bkstellar_obj)
    use com_main_gw
    implicit none
    type(particle_samples_arr_type)::bks, bkstellar_obj(n_tot_comp_sg)
    integer stellar_type,i
   
   do i=1, n_tot_comp_sg
        call get_type_from_ctl_obidx_sg(i,stellar_type)
        call sams_arr_select_type_single(bks, bkstellar_obj(i), stellar_type) 
   end do
end subroutine
 

subroutine gen_nxj_fxj_gx_iregular(dm,spp)
    use com_main_gw
    implicit none
    type(diffuse_mspec)::dm
    type(star_pot_para)::spp
    type(particle_samples_arr_type)::bkstar, bksbh, bkwd, bkns
	type(particle_samples_arr_type)::bkbd
    integer ierr,i
    logical::norm_in
    real(8) t1, t2

    if(rid.eq.0)then
        call cpu_time(t1)
    end if


    !print*, "sample_nxgx_logemin,sample_nxgx_logemax=",sample_nxgx_logemin,sample_nxgx_logemax
    call set_barp(dm,spp)
    if(rid.eq.0.and.ctl%chattery.ge.2)then
        call dm%barp_ir%print("nxgx:barp_ir")
    end if
    if(ctl%barge_evl_method.eq.barge_evl_method_grid_2d)then
        call set_gx_nx2d_ranges_ir(dm)
        !print*, "1"
        call set_nx_1d_ranges_ir(dm)
        !print*, "2"
        call get_nx_ir_simu(dm)
        !print*, "3"
        call get_nxj_ir(dm)
        !print*, "4"
        call get_fxj0_ir(dm)
         
        call get_barge0(dm)
    else
        call set_gx_nx_ranges_ir(dm)
        call get_nx_ir_simu(dm) 
        call get_barge0(dm)
    end if
    !if(rid.eq.0.and.ctl%chattery.ge.2)then
    !    call dm%all%all%barge_ir%print("all:barge")
    !end if
	
    if(ctl%barge_evl_method.eq.barge_evl_method_direct)then
        call set_gx_ranges_ir(dm)
    end if
    if(rid.eq.0)then
        call cpu_time(t2)
        print*, "gen_nxj_fxj_gx_iregular, used time=", t2-t1, " s"
    end if
end subroutine

subroutine print_num_all(dm)
    use com_main_gw
    implicit none
    type(diffuse_mspec)::dm
    integer i
    real(8) nb(ctl%m_bins,10),nbw(ctl%m_bins,10)
    
    nb=0
    write(*,fmt=*)  "print_num_all=========================="
    write(*,fmt="(A3, 10A13)") "i","mc", "star", "sbh", "ns", "wd","bd", "rg", "nHe", "dm"!, "bbh"
    do i=1, dm%n
        call get_num_all(nb(i,1),  nbw(i,1),dm%mb(i)%star)
        call get_num_all(nb(i,2),  nbw(i,2),  dm%mb(i)%sbh)
        ! call get_num_all(nb(i,3),  nbw(i,3), dm%mb(i)%bbh)
		call get_num_all(nb(i,3),  nbw(i,3), dm%mb(i)%ns)
		call get_num_all(nb(i,4),  nbw(i,4), dm%mb(i)%wd)
		call get_num_all(nb(i,5),  nbw(i,5), dm%mb(i)%bd)
        call get_num_all(nb(i,6),  nbw(i,6), dm%mb(i)%rg)
        call get_num_all(nb(i,7),  nbw(i,7), dm%mb(i)%nakedHe)
        call get_num_all(nb(i,8),  nbw(i,8), dm%mb(i)%dark_matter)
	end do	
	do i=1, dm%n
        if(sum(nb(i,1:8)).ne.0)then
    	    write(*,fmt="(I3,F13.2, 10I13)") i, dm%mb(i)%mc, ctl%bin_mass_simulation_particle_number_comp_tot(i,1:n_tot_comp_sg)
        end if
	end do
	write(*,fmt="(A3, 8A13)") "i","mc", "starw", "sbhw", "nsw", "wdw","bdw", "rgw", "nHe", "dmw"!, "bbhw"
	do i=1, dm%n
        if(sum(nbw(i,1:8)).ne.0)then
            write(*,fmt="(I3,10E13.2)") i, dm%mb(i)%mc, nbw(i,1:8)
        end if
    end do
	!write(*,fmt="(A3, 8A13)") "i","mc","j>0.5","","","j<0.5"
    !do i=1, dm%n
    !    write(*,fmt="(I3,7F13.2)") i, dm%mb(i)%mc, nbj1(i,1:3), nbj2(i,1:3)
    !end do
    write(*,fmt=*) "end of print_num_all==================="
end subroutine
subroutine print_norm_dms(dm,out_unit)
    use md_dms
    use model_basic,only:ctl
    implicit none 
    type(diffuse_mspec)::dm
    integer i,out_unit
    real(8) w
    write(out_unit, fmt="(2A10,12A14)") "N", "mc", "rNstar", "rNsbh","rNns", "rNwd", "rNbd", "rRg" 
         !,"rNbstar", "rNbbh"
    do i=1, dm%n
        write(out_unit, fmt="(I10, F10.2, 12F14.2)") i, dm%mb(i)%mc,&
            dm%mb(i)%star%n_real, dm%mb(i)%sbh%n_real, &
            dm%mb(i)%ns%n_real, dm%mb(i)%wd%n_real, dm%mb(i)%bd%n_real, dm%mb(i)%rg%n_real 
    end do
    w=0
    !print*, "1/0d0>1d99", 1d0/w>1d90
    write(out_unit, fmt="(2A10,12A14)") "N", "mc", "sNstar", "sNsbh", &
        "sNns", "sNwd", "sNbd", "sNrg" 
    do i=1, dm%n
        write(out_unit, fmt="(I10, F10.2, 12I14)") i, dm%mb(i)%mc,&
            ctl%bin_mass_simulation_particle_number_comp_tot(i,1:6)
    end do
    write(out_unit, fmt="(A8,12I10)") "N:", ctl%bin_mass_simulation_particle_number_tot(1:dm%n)
end subroutine
 
subroutine print_num_boundary(dm)
    use com_main_gw
    implicit none
    type(diffuse_mspec)::dm
    integer i
    real(8) nb(ctl%m_bins,10),nbw(ctl%m_bins,10),nbj1(ctl%m_bins,10),nbj2(ctl%m_bins,10)
    
    nb=0
    write(*,fmt=*)  "print_num_boundary=========================="
    write(*,fmt="(A3, 10A13)") "i","mc", "star", "sbh", "ns", "wd","bd", "rg"
    do i=1, dm%n
        call get_num_boundary(nb(i,1), nbj1(i,1), nbj2(i,1), nbw(i,1),dm%mb(i)%star)
        call get_num_boundary(nb(i,2), nbj1(i,2), nbj2(i,2), nbw(i,2),  dm%mb(i)%sbh) 
		call get_num_boundary(nb(i,3), nbj1(i,3), nbj2(i,3), nbw(i,3), dm%mb(i)%ns)
		call get_num_boundary(nb(i,4), nbj1(i,4), nbj2(i,4), nbw(i,4), dm%mb(i)%wd)
		call get_num_boundary(nb(i,5), nbj1(i,5), nbj2(i,5), nbw(i,5), dm%mb(i)%bd)
        call get_num_boundary(nb(i,6), nbj1(i,6), nbj2(i,6), nbw(i,6), dm%mb(i)%rg)
        call get_num_boundary(nb(i,7), nbj1(i,7), nbj2(i,7), nbw(i,7), dm%mb(i)%nakedHe)
	end do	
	do i=1, dm%n
    	write(*,fmt="(I3,7F13.2)") i, dm%mb(i)%mc, nb(i,1:8)
	end do
	write(*,fmt="(A3, 8A13)") "i","mc", "starw", "sbhw",  "nsw", "wdw","bdw", "rgw", "nHew"
	do i=1, dm%n
        write(*,fmt="(I3,7F13.2)") i, dm%mb(i)%mc, nbw(i,1:8)
    end do 
    write(*,fmt=*) "end of print_num_boundary==================="
end subroutine
subroutine get_num_all(nb,nbw, so)
    use com_main_gw
    implicit none
    real(8) nb, nbw
    type(dms_stellar_object)::so
    integer i
    nb=0;nbw=0; 
    do i=1, so%n
        if(so%nejw(i)%e>log10(ctl%x_boundary).and.so%nejw(i)%e<log10(emax_factor))then
            nbw=nbw+so%nejw(i)%w
            nb=nb+1
					
        end if
    end do
end subroutine
subroutine get_num_boundary(nb,nbj1,nbj2,nbw, so)
    use com_main_gw
    implicit none
    real(8) nb, nbw, nbj1,nbj2
    type(dms_stellar_object)::so
    integer i
    nb=0;nbw=0;nbj1=0; nbj2=0
    do i=1, so%n
        if(so%nejw(i)%e<log10(ctl%x_boundary))then
            nbw=nbw+so%nejw(i)%w
            nb=nb+1
			if(so%nejw(i)%j>0.5d0)then
				nbj1=nbj1+1
			else
				nbj2=nbj2+1
			end if	
        end if
		
        if(so%nejw(i)%e<log10emin_factor)then
            print*, "error!??:", so%nejw(i)%e, log10emin_factor
        end if
    end do
end subroutine
 
  
 

subroutine update_weights()
    use com_main_gw
    implicit none

    call set_clone_weight(bksams)
    !call set_asym_weight(sms_single)
    call set_real_weight(bksams)
    
end subroutine

subroutine update_arrays_single(clean_after)
    use com_main_gw
	implicit none
    logical clean_after
    !print*, "1"
	call all_chain_to_arr_single(bksams,bksams_arr)
    print*, "bksams_arr_all%n=", bksams_arr%n,rid
	call set_sample_arr_indexs_rid_particle(bksams_arr,rid)
    !print*, "3"
	call convert_sams_pointer_arr(bksams, bksams_pointer_arr,type=1)
    !print*, "4"
	call bksams_arr%select(bksams_arr_norm, exit_normal, -1, -1d0, -1d0)
    print*, "bksams_arr_norm%n=", bksams_arr_norm%n,rid 
    if(clean_after)then
        if(allocated(bksams_arr%sp))deallocate(bksams_arr%sp)
    end if
    !call check_bksams()
end subroutine	
 
subroutine set_spt_y2(spp)
    !use model_basic
    use constant
    use md_star_pot
    implicit none
    !type(diffuse_mspec)::dm
    type(star_pot_para)::spp
    integer n, i
    real(8) dphidr,dphidr2,beta_tmp
    real(8) x,rho_tmp,r,phi, dlnphi2dlnr2, dlnphidlnr
    
    
    n=spp%fphi_star%nbin
    call spp%fphi_star%prepare_spline()
    !print*, dm%fphi_star%y2(1:n)
    do i=1, n
        x=spp%fphi_star%xb(i)
        r=10**x
        phi=10**spp%fphi_star%fx(i)
        call get_beta_full_range(spp,x, beta_tmp)
        call get_rho_full_range_spp(spp,x,rho_tmp)
        dphidr=-beta_tmp/r**2
        dphidr2=2*beta_tmp/r**3-rho_tmp*4*pi
     !   if(i>1.and.i<n-1)then
    !        print*,r, log(10d0)*r**2/phi*(dphidr/r-dphidr**2/phi+dphidr2),dm%fphi_star%y2(i), &
    !        (dm%fphi_star%fx(i+1)+dm%fphi_star%fx(i-1)-dm%fphi_star%fx(i)*2)/dm%fphi_star%xstep**2
    !    end if
        spp%fphi_star%y2(i)=log(10d0)*r**2/phi*(dphidr/r-dphidr**2/phi+dphidr2)
    end do
    !stop
end subroutine

subroutine set_fma_y2(spp)
    !use model_basic
    use constant
    use md_star_pot
    implicit none
    !type(diffuse_mspec)::dm
    type(star_pot_para)::spp
    integer n, i
    real(8) dphidr,dphidr2,beta_tmp
    real(8) x,rho_tmp,r,fma, dlnphi2dlnr2, dlnphidlnr
    
    
    n=spp%fma_star%nbin
    call spp%fma_star%prepare_spline()
    !print*, dm%fphi_star%y2(1:n)
    do i=1, n
        x=spp%fma_star%xb(i)
        r=10**x
       ! fma=10**spp%fma_star%fx(i)
        call get_rho_full_range_spp(spp,x,rho_tmp)
        spp%fma_star%y2(i)=log(10d0)*r**3*4*pi*rho_tmp
    end do
    !stop
end subroutine

subroutine get_diffuse_mspec_rtables(dm,spp)
    use md_dms
    use model_basic,only:ctl
    use md_star_pot
    use MPI_comu,only:rid
    implicit none
    type(diffuse_mspec)::dm
    type(star_pot_para)::spp
 
    if(spp%mbh.ne.0)then
        
        call get_dms_alpha(dm,spp)
        !print*, "3"
        if(ctl%chattery.ge.2.and.rid.eq.0)then
            call dm%alpha_r%print("alpha")
        end if
    end if
    call get_radius_energy_transition(dm,spp)
	!call get_spp_beta(dm)
    
end subroutine
subroutine get_radius_energy_transition(dm,spp)
    use md_dms
    use md_star_pot
    use model_basic
    implicit none
    type(diffuse_mspec)::dm
    type(star_pot_para)::spp
    type(s1d_type)::inv_alpha
    integer niter,ier
    real(8) emax,rtbis_yacc,par(50), xl, xh

    if(spp%mbh_dmless.eq.0)then
        call get_phi_star_full_range(spp,ctl%log10rmin_factor,emax)
        energy_transition=10**emax
    else
        !inv_alpha=dm%alpha_r
        !inv_alpha%fx=dm%alpha_r%xb
        !inv_alpha%xb=dm%alpha_r%fx
        !call inv_alpha%get_value_s(1d0,radius_transition)
        xl=ctl%log10rmin_factor; xh=ctl%log10rmax_factor
        radius_transition= rtbis_yacc(func,xl,xh,1d-10,par,niter,1000,ier,.true.)
        energy_transition=spp%mbh_dmless/10**radius_transition
    end if
contains 
real(8) function  func(x,par)
    implicit none
    real(8) x, par(50), phi_tmp
    call get_phi_star_full_range(spp,x,phi_tmp)
    func=spp%mbh_dmless/10**x-10**phi_tmp
end function
end subroutine
 
subroutine collect_dms_density_simu(dm)
    use com_main_gw
	implicit none
	integer ierr
	integer i,j, nx
    type(diffuse_mspec)::dm
    !real(8) t1, t2

    !if(rid.eq.0)then
    !    call cpu_time(t1)
    !end if

    nx=dm%mb(1)%all%fden_simu%nbin
    
    do i=1, dm%n
        call collection_and_avg_fx(dm%mb(i)%all%fden_simu%fx,nx)
        call collection_and_avg_fx(dm%mb(i)%all%fna_simu%fx,nx)
        call collection_and_avg_fx(dm%mb(i)%all%fma_simu%fx,nx)
        do j=1, n_tot_comp
            call collection_and_avg_fx(dm%mb(i)%dsp(j)%p%fden_simu%fx,nx)
            call collection_and_avg_fx(dm%mb(i)%dsp(j)%p%fna_simu%fx,nx)
            call collection_and_avg_fx(dm%mb(i)%dsp(j)%p%fma_simu%fx,nx)
        end do
    end do
    call collection_and_avg_fx(dm%all%all%fden_simu%fx,nx)
    call collection_and_avg_fx(dm%all%all%fna_simu%fx,nx)
    call collection_and_avg_fx(dm%all%all%fma_simu%fx,nx)

    !if(rid.eq.0)then
    !    call cpu_time(t2)
    !    print*, "density colletion used time: ", t2-t1, " s"
    !end if
end subroutine	

subroutine collect_dms_fmr(dm)
    use com_main_gw
	implicit none
	integer ierr
	integer i,j, nx
    type(diffuse_mspec)::dm
    !real(8) t1, t2

    !if(rid.eq.0)then
    !    call cpu_time(t1)
    !end if

    nx=dm%all%all%fm_dstr%nbin
    !do i=1, ctl%ntasks
    !    if(rid.eq.i-1)then
    !        print*, "i=",i
    !        call dm%all%all%fden%print("fden")
    !    end if
    !    call mpi_barrier(mpi_comm_world,ierr)
    !end do
    !stop

    !do i=1, ctl%ntasks
    !    if(rid.eq.i)then
    !        do j=1, dm%n
    !            call dm%mb(j)%all%fmden%print("fmden")
    !        end do
    !    end if
    !end do


    do j=1, n_tot_comp
        call collection_and_avg_fx(dm%all%dsp(j)%p%fm_dstr%fxw,nx)
        call collection_and_avg_fx(dm%all%dsp(j)%p%fr_dstr%fxw,nx)
    end do
    call collection_and_avg_fx(dm%all%all%fm_dstr%fxw,nx)
    call collection_and_avg_fx(dm%all%all%fr_dstr%fxw,nx)
end subroutine	

subroutine collect_dms_density(dm)
    use com_main_gw
	implicit none
	integer ierr
	integer i,j, nx
    type(diffuse_mspec)::dm
    !real(8) t1, t2

    !if(rid.eq.0)then
    !    call cpu_time(t1)
    !end if

    nx=dm%mb(1)%all%fmden%nbin
    !do i=1, ctl%ntasks
    !    if(rid.eq.i-1)then
    !        print*, "i=",i
    !        call dm%all%all%fden%print("fden")
    !    end if
    !    call mpi_barrier(mpi_comm_world,ierr)
    !end do
    !stop

    !do i=1, ctl%ntasks
    !    if(rid.eq.i)then
    !        do j=1, dm%n
    !            call dm%mb(j)%all%fmden%print("fmden")
    !        end do
    !    end if
    !end do

    do i=1, dm%n
        call collection_and_avg_fx(dm%mb(i)%all%fmden%fx, nx)
        call collection_and_avg_fx(dm%mb(i)%all%fden%fx,nx)
        call collection_and_avg_fx(dm%mb(i)%all%fna%fx,nx)
        call collection_and_avg_fx(dm%mb(i)%all%fma%fx,nx)
        do j=1, n_tot_comp
            call collection_and_avg_fx(dm%mb(i)%dsp(j)%p%fmden%fx,nx)
            call collection_and_avg_fx(dm%mb(i)%dsp(j)%p%fden%fx,nx)
            call collection_and_avg_fx(dm%mb(i)%dsp(j)%p%fna%fx,nx)
            call collection_and_avg_fx(dm%mb(i)%dsp(j)%p%fma%fx,nx)
        end do
    end do
    call collection_and_avg_fx(dm%all%all%fmden%fx,nx)
    call collection_and_avg_fx(dm%all%all%fden%fx,nx)
    

    call collection_and_avg_fx(dm%all%all%fna%fx,nx)
    call collection_and_avg_fx(dm%all%all%fma%fx,nx)
    if(ctl%barge_evl_method.eq.barge_evl_method_direct)then
        call collect_dms_density_simu(dm)
    end if

    !if(rid.eq.0)then
    !    call cpu_time(t2)
    !    print*, "density colletion used time: ", t2-t1, " s"
    !end if
end subroutine	

subroutine collect_dms_gx(dm)
    use com_main_gw
	implicit none
	integer ierr
	integer i,j,k,nx, ny
    type(diffuse_mspec)::dm
    real(8),allocatable::fx(:)
    integer  nb(dms%barp_ir%nbin)
    !real(8) t1, t2
    !if(rid.eq.0)then
    !    call cpu_time(t1)
    !end if
    !print*, "start of collection, rid=",rid
    select case(ctl%barge_evl_method) 
    case(barge_evl_method_grid_2d)
        nx=dm%mb(1)%all%nx_ir%nbin 
        do i=1, dm%n
            call collection_and_avg_fx(dm%mb(i)%all%nx_ir%fxw,nx)
            !print*, "1"
            do j=1, n_tot_comp
                call collection_and_avg_fx(dm%mb(i)%dsp(j)%p%nx_ir%fxw,nx)
            end do
        end do
        call collection_and_avg_fx(dm%all%all%nx_ir%fxw,nx)
         
        nx=dm%mb(1)%all%gxj_ir%nx
        ny=dm%mb(1)%all%gxj_ir%ny
        !print*, "nx,ny=",nx,ny
        do i=1, dm%n
            call collection_and_avg_s2d(dm%mb(i)%all%gxj_ir%fxy,nx, ny)

            do j=1, n_tot_comp
                call collection_and_avg_s2d(dm%mb(i)%dsp(j)%p%gxj_ir%fxy,nx, ny)
            end do
        end do
        call collection_and_avg_s2d(dm%all%all%gxj_ir%fxy,nx, ny)

        do i=1, dm%n
            call collection_and_avg_s2d(dm%mb(i)%all%nxj_ir%nxyw,nx, ny)
            do j=1, n_tot_comp
                call collection_and_avg_s2d(dm%mb(i)%dsp(j)%p%nxj_ir%nxyw,nx, ny)
            end do
        end do
        call collection_and_avg_s2d(dm%all%all%nxj_ir%nxyw,nx, ny)
  
        do i=1, dm%n 
            call collection_fx_int(dm%mb(i)%all%nx_ir%nb,nx)
            !call collection_int(dms%mb(i)%all%n,sample_threads_number)
            call collection_and_avg_fx(dm%mb(i)%all%barge_ir%fx,nx) 
            do j=1, n_tot_comp 
                call collection_and_avg_fx(dm%mb(i)%dsp(j)%p%barge_ir%fx,nx) 
            end do
        end do 
        nb=dm%all%all%nx_ir%nb
        call collection_fx_int(nb,nx)
         
        call collection_and_avg_fx(dm%all%all%barge_ir%fx,nx)
        
		if(rid.eq.0.and.ctl%chattery.ge.2)then
			if(ctl%chattery.ge.3)then
				call dm%mb(1)%all%nxj_ir%print("nxj ir m1")
			end if
			do i=1, dm%n
				print*, "i=",i
				call dm%mb(i)%all%nx_ir%print("nx m1")
				call dm%mb(i)%all%barge_ir%print("barge ir m1")
			end do
			
			!
		end if 
    end select
end subroutine	

