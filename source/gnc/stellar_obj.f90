module md_stellar_object
    use com_sts_type
    ! use md_bk_species
    use md_coeff
	use md_chain_pointer
    use md_sts_d21
    use, intrinsic :: ieee_arithmetic
    type nejw_type
        real(8) e, j, w, rp, ra, pd, jc, m, rad   
        integer idx
    end type
    type dms_stellar_object
        integer:: n=0
        real(8)::n_real=0
        type(nejw_type),allocatable::nejw(:)
        type(s1d_type)::fden_simu, fden  
		type(s1d_type)::fmden  !fmden represents the mass density
        type(s1d_type)::fNa, fNa_simu, fMa, fMa_simu, faniso, fslope
		type(s1d_hst_type)::fm_dstr   ! mass distribution function, only for models with mass spectrums
		type(s1d_hst_type)::fr_dstr   ! stellar radius distribution function, only for models with mass spectrums
		type(s1d_ird_type)::barge_ir!, barge_norm
		type(s1d_type)::barge
		type(s2d_type):: gxj, gxjcr, gxj_ir, gxr_ir
		type(sts_fc_type)::barj
		type(s2d_hst_type)::nxj
		type(s1d_hst_type)::nx
		type(s1d_hst_ird_type)::nx_ir
		type(s2d_hst_ird_type)::nxj_ir
        real(8) asymp, spt_rho_rmin, spt_rho_rmax
        contains
			procedure::init_r=>init_dms_stellar_object_r
			procedure::init_e=>init_dms_stellar_object_e
		!procedure::init=>init_dms_stellar_object
            procedure::deallocation=>deallocate_dms_stellar_object            
			procedure::get_n=>dms_so_get_n
			procedure::get_mtot=>dms_so_get_mtot
    end type
	type dso_pointer
		type(dms_stellar_object),pointer::p
	end type
    private::deallocate_dms_stellar_object,dms_so_get_n,dms_so_get_mtot
	private::init_dms_stellar_object_r,init_dms_stellar_object_e
    integer,parameter::dj_n=8, dj_n2=8
contains
	subroutine init_dms_stellar_object_r(so, df_coe_bins, dstr_bins,&
		rmin,rmax, jmin,jmax, r0_cl,jb_type)
		implicit none
		class(dms_stellar_object)::so
		integer df_coe_bins, dstr_bins, i
		real(8) emin, emax , r0_cl
		real(8) rmin, rmax,jmin,jmax,tmin,tmax
		integer jb_type
		select case(jb_type) 
		case(jbin_type_log)
			tmin=log10(jmin)
			tmax=log10(jmax) 
		end select

		call so%fden%init(rmin, rmax, dstr_bins, coeff_sts_type_dc)
		call so%fmden%init(rmin, rmax, dstr_bins, coeff_sts_type_dc)
		call so%faniso%init(rmin,rmax,dstr_bins,coeff_sts_type_dc)
		call so%fslope%init(rmin,rmax,dstr_bins,coeff_sts_type_dc)
		call so%fden_simu%init(rmin, rmax, dstr_bins, coeff_sts_type_dc)

		call so%fden%set_range()
		call so%fmden%set_range()
		call so%faniso%set_range()
		call so%fslope%set_range()
		call so%fden_simu%set_range()

		call so%fNa%init(rmin, rmax, dstr_bins, coeff_sts_type_dc)
		call so%fNa%set_range()
		call so%fNa_simu%init(rmin, rmax, dstr_bins, coeff_sts_type_dc)
		call so%fNa_simu%set_range()
		!so%fden%nsam=1;
		call so%fMa%init(rmin, rmax, dstr_bins, coeff_sts_type_dc)
		call so%fMa%set_range()
		call so%fMa_simu%init(rmin, rmax, dstr_bins, coeff_sts_type_dc)
		call so%fMa_simu%set_range()
		!so%fMa%nsam=1
		!print*, "2"

		call so%barj%init(0d0, 1d0, dstr_bins, fc_spacing_linear, use_weight=.true.)
		call so%barj%set_range()

		!so%n=0
		!so%n_real=0
		!print*, "xxx"
	end subroutine
    
    subroutine init_dms_stellar_object_e(so, df_coe_bins, dstr_bins,&
        emin_in, emax_in,   jmin, jmax, r0_cl,eb_type, jb_type)
        implicit none
        class(dms_stellar_object)::so
        integer df_coe_bins, dstr_bins, i
        real(8) emin_in, emax_in, jmin, jmax, r0_cl
        real(8) rmin, rmax, smin,smax
		integer jb_type,eb_type
		real(8) tmin, tmax

		select case(jb_type) 
        case(jbin_type_log)
            tmin=log10(jmin)
            tmax=log10(jmax) 
        end select

		if(eb_type.eq.ebin_type_log)then
			smin=log10(emin_in)
			smax=log10(emax_in) 
		end if
        !rmin=log10(0.5d0*r0_cl/(10**smax));rmax=log10(0.5d0*r0_cl/(10**emin))
		!print*, "emin,smax,rmin,rmax,tmin,tmax=",emin,smax,rmin,rmax,tmin,tmax
        call so%barge%init( smin, smax, dstr_bins,sts_type_dstr)
        call so%barge%set_range()
		
		!call so%barge_ir%init( smin, smax, dstr_bins,sts_type_dstr)

		call so%nx%init(smin,smax,dstr_bins,use_weight=.true.)
		call so%nx%set_range()

		!call so%nx_ir%init(smin,smax,dstr_bins,use_weight=.true.)

        call so%gxj%init(df_coe_bins, df_coe_bins, smin, smax, tmin,tmax, sts_type_dstr)
        call so%gxj%set_range()
		!print*, "2.1"
		!call so%gxjcr%init(df_coe_bins*2, df_coe_bins*2, smin, smax, tmin,tmax, sts_type_dstr)
		!call so%gxjcr%set_range()
		!print*, "2.2"
		!call so%gxr_ir%init(df_coef_bins, dstr_bins, smin,smax, ,sts_type_dstr )
		!call so%gxr_ir%set_range()

		call so%nxj%init(df_coe_bins, df_coe_bins, smin, smax, tmin,tmax, use_weight=.true.)
        call so%nxj%set_range()
		!print*, "3"
        deallocate(so%nxj%da_y_proj, so%nxj%fc_y_proj)
        allocate(so%nxj%da_y_proj(dj_n),so%nxj%fc_y_proj(dj_n))
        do i=1, dj_n
            call so%nxj%fc_y_proj(i)%init(smin,smax,so%gxj%ny,fc_spacing_linear)
        end do
        call so%barj%init(0d0, 1d0, df_coe_bins, fc_spacing_linear, use_weight=.true.)
        call so%barj%set_range()

		!print*, "xxx"
    end subroutine
    
	! subroutine dms_so_get_faniso(so)
	! 	implicit none
	! 	class(dms_stellar_object)::so
	! 	!real(8) e0
	! 	integer i
    !     real(8) en(so%n), jm(so%n), wm(so%n)

	! 	if(so%n>0)then
    !         en(1:so%n)=10**so%nejw(1:so%n)%e!*e0
    !         jm(1:so%n)=so%nejw(1:so%n)%j
    !         wm(1:so%n)=so%nejw(:)%w
	! 		call get_anisotropy_mc(en, jm ,wm,so%n,&
	! 		so%faniso%xmin,so%faniso%xmax, so%faniso%nbin, &
	! 		so%faniso%xb(:),so%faniso%fx(:))
	! 	end if
	! end subroutine
	
    subroutine deallocate_dms_stellar_object(so)
        implicit none
        class(dms_stellar_object)::so
        call so%fden%deallocate()
        call so%fna%deallocate()
        call so%fma%deallocate()
        call so%barj%deallocate()
        call so%barge%deallocate()
		call so%fslope%deallocate()
		call so%faniso%deallocate()
    end subroutine
    

    subroutine get_asymp_norm_factor_one(dso,  x_boundary, norm)
        implicit none
        class(dms_stellar_object)::dso
        real(8) cnorm, norm,x_boundary
!
        !call dm%all%star%barge%print("asymp_norm_allstar")
        !print*, log10(dm%x_boundary),dm%x_boundary
        if(dso%n_real>0)then
            !print*, log10(x_boundary),x_boundary
			!call dso%barge%print("norm dso%barge")
            !call get_value_at_x_fc(dso%barge, log10(x_boundary), cnorm, 1)
			call dso%barge%get_value_l(log10(x_boundary), cnorm)
            if( ieee_is_nan(cnorm).or. cnorm.eq.0)then
                print*, "star:cnorm is nan or 0", cnorm
                call dso%barge%print("dso%barge")
                stop
            !    call dm%barge0%print()
            !    print*, "va=", log10(0.25d0), log10(0.25d0)>dm%barge0%xmin, &
            !    log10(0.25d0)<dm%barge0%xmax
            !    !stop
            !    maxidx=maxloc(dm%barge0%fx)
            !    norm=maxval(dm%barge0%fx)*(dm%barge0%xb(maxidx(1)))**(-1d0/4d0)
            !    !print*, "norm=", norm
            else
                norm=dso%asymp/cnorm
				!print*, "norm=", norm, dso%asymp
				!stop
                !call dso%barge%print()
            end if
        else
            print*, "error! dso%n_real=", dso%n_real
            stop
        end if
    end subroutine
    subroutine normalize_barge_one(dso, norm)
        implicit none
        class(dms_stellar_object)::dso
        real(8) norm
        integer i
        
        !dso%barge%nsam=dso%n
        dso%barge%fx=dso%barge%fx*norm
        
    end subroutine
	subroutine normalize_gxj_one(dso, norm)
        implicit none
        class(dms_stellar_object)::dso
        real(8) norm
        integer i
        
       ! dso%gxj%nsam=dso%n
        dso%gxj%fxy=dso%gxj%fxy*norm
        
    end subroutine
	subroutine dms_so_get_fxj_spt(so, n0, pd,jc, r0,jbtype)
        implicit none
        class(dms_stellar_object)::so
        integer i, j,jbtype
        real(8) jm ,x, n0,  r0
		type(s2d_type)::pd
        type(s1d_type) jc
		real(8) jc_xy, pd_xy

		if(so%n.eq.0) return
		select case(jbtype) 
		case(Jbin_type_log)
			do i=1, so%nxj%nx
				x=10**so%nxj%xcenter(i)
				do j=1, so%nxj%ny
					jm=10**so%nxj%ycenter(j)
					call jc%get_value_s(so%nxj%xcenter(i),jc_xy)
					call pd%get_value_l(so%nxj%xcenter(i),so%nxj%ycenter(j),pd_xy)
					!so%gxj%fxy(i,j)=so%nxj%nxyw(i,j)/(x*log(10d0))&
					!/so%nxj%xstep/so%nxj%ystep &
					!*pi**(-0.5d0)*2**(-1.5d0)/r0**3/(jm**2*log(10d0))/n0/jc%fx(i)**2/pd%fxy(i,j)
					so%gxj%fxy(i,j)=so%nxj%nxyw(i,j)/(x*log(10d0))&
					/so%nxj%xstep/so%nxj%ystep &
					*pi**(-0.5d0)*2**(-1.5d0)/r0**3/(jm**2*log(10d0))/n0/jc_xy**2/pd_xy

					!print*, "xystep=",so%nxj%nxyw(i,j), so%nxj%xstep, so%nxj%ystep
					!read(*,*)
				end do 
			end do 
		case default
			print*, "fxj error!"
			stop
		end select
    end subroutine
	subroutine dms_so_get_n(so,n)
		implicit none
		class(dms_stellar_object)::so
		real(8) n
		integer i
		n=0
		do i=1, so%n
			n=n+so%nejw(i)%w
		end do
	end subroutine
	subroutine dms_so_get_mtot(so,mtot)
		implicit none
		class(dms_stellar_object)::so
		real(8) mtot
		integer i
		mtot=0
		do i=1, so%n
			mtot=mtot+so%nejw(i)%w*so%nejw(i)%m
		end do
	end subroutine
    subroutine dms_so_get_fxj(so, n0, mbhin, v0,jbtype)
        implicit none
        class(dms_stellar_object)::so
        integer i, j,jbtype
        real(8) jm ,x, n0, mbhin, v0
		if(so%n.eq.0) return
		select case(jbtype) 
		case(Jbin_type_log)
			do i=1, so%nxj%nx
				x=10**so%nxj%xcenter(i)
				do j=1, so%nxj%ny
					jm=10**so%nxj%ycenter(j)
					so%gxj%fxy(i,j)=so%nxj%nxyw(i,j)/(x*log(10d0))&
					/so%nxj%xstep/so%nxj%ystep &
					*pi**(-1.5d0)*v0**6*x**2.5d0/(jm**2*log(10d0))/n0/mbhin**3
					!print*, "xystep=",so%nxj%xstep, so%nxj%ystep
					!read(*,*)
				end do 
			end do 
		case default
			print*, "fxj error!"
			stop
		end select
    end subroutine
    
end module
module md_mass_bins 
	use md_stellar_object

	integer,parameter:: n_tot_comp_sg=8, n_tot_comp_by=0
	integer,parameter:: n_tot_comp=n_tot_comp_sg+n_tot_comp_by
	type mass_bins
		type(dms_stellar_object):: star, sbh, ns, all, wd, bd, rg, dark_matter, nakedHe
		type(dso_pointer)::dsp(n_tot_comp)

		real(8) mc, m1, m2 !m1<mass<m2
		integer df_coe_bins, dstr_bins
		real(8) frac   ! the fraction of the mass over all mass bins
		real(8) emin, emax, jmin, jmax, mtot, v0, n0, barmin, r0_cl
		!type(sts_fc_type):: n_collid_11 ! expected collision rates
		!type(sts_fc_type):: n_collid_gw ! expected collision rates
		type(diffuse_coeffient_type)::dc
		!type(s2d_hst_type)::dt
		contains
		procedure::init_r=>init_mass_bins_r
		procedure::init_e=>init_mass_bins_e
		!procedure::init=>init_mass_bins
		procedure::write_mb=>write_info_mass_bin
		procedure::read_mb=>read_info_mass_bin
		!procedure::get_nEJ0_fEJ0
		!procedure::get_barge
	end type
contains
	subroutine init_mass_bins_r(mb, df_coe_bins, dstr_bins,  rmin,rmax, &
		jmin,jmax,mbhin, v0, n0, r0_cl,jb_type)
		implicit none
		class(mass_bins),target::mb
		integer df_coe_bins, dstr_bins
		real(8) emin, emax, jmin, jmax, mbhin, n0, v0, rmin, rmax, r0_cl
		integer  i,jb_type
		integer icout

		mb%df_coe_bins=df_coe_bins; mb%dstr_bins=dstr_bins
		mb%jmin=jmin; mb%jmax=jmax
		mb%mtot=mbhin; mb%v0=v0; mb%n0=n0; mb%r0_cl=r0_cl

		!rmin=log10(0.5d0*r0_cl/(10**emax));rmax=log10(0.5d0*r0_cl/(10**emin))
 
		call mb%all%init_r(df_coe_bins, dstr_bins,  rmin,rmax, jmin,jmax,r0_cl,jb_type) 

		icout=1
		mb%dsp(icout)%p=>mb%star
		icout=icout+1
		mb%dsp(icout)%p=>mb%sbh
		icout=icout+1
		mb%dsp(icout)%p=>mb%ns
		icout=icout+1
		mb%dsp(icout)%p=>mb%wd
		icout=icout+1
		mb%dsp(icout)%p=>mb%bd  
		icout=icout+1
		mb%dsp(icout)%p=>mb%rg
		icout=icout+1
		mb%dsp(icout)%p=>mb%dark_matter
		icout=icout+1
		mb%dsp(icout)%p=>mb%nakedHe
		! icout=icout+1
		! mb%dsp(icout)%p=>mb%bstar
		! icout=icout+1
		! mb%dsp(icout)%p=>mb%bbh		
		

		
		!print*, "mass bin 2"
		do i=1, n_tot_comp
			call mb%dsp(i)%p%init_r(df_coe_bins, dstr_bins, rmin,rmax,  jmin,jmax, r0_cl,jb_type)
		end do
		!call mb%n_collid_11%init(rmin, rmax, mb%dstr_bins, fc_spacing_linear)
		!call mb%n_collid_gw%init(rmin, rmax, mb%dstr_bins, fc_spacing_linear)
		!call mb%n_collid_11%set_range()
		!call mb%n_collid_gw%set_range()

	end subroutine
	subroutine init_mass_bins_e(mb, df_coe_bins, dstr_bins, emin, emax,  jmin, jmax,&
	mbhin, v0, n0, r0_cl,eb_type,jb_type)
		implicit none
		class(mass_bins),target::mb
		integer df_coe_bins, dstr_bins
		real(8) emin, emax, jmin, jmax, mbhin, n0, v0, rmin, rmax, r0_cl
		integer jb_type, i, eb_type,icout

		mb%df_coe_bins=df_coe_bins; mb%dstr_bins=dstr_bins
		mb%emin=emin; mb%emax=emax; mb%jmin=jmin; mb%jmax=jmax
		mb%mtot=mbhin; mb%v0=v0; mb%n0=n0; mb%r0_cl=r0_cl
 
		call mb%all%init_e(df_coe_bins, dstr_bins,emin, emax,  jmin, jmax, r0_cl, eb_type, jb_type) 

		icout=1
		mb%dsp(icout)%p=>mb%star
		icout=icout+1
		mb%dsp(icout)%p=>mb%sbh
		icout=icout+1
		mb%dsp(icout)%p=>mb%ns
		icout=icout+1
		mb%dsp(icout)%p=>mb%wd
		icout=icout+1
		mb%dsp(icout)%p=>mb%bd  
		icout=icout+1
		mb%dsp(icout)%p=>mb%rg
		icout=icout+1
		mb%dsp(icout)%p=>mb%dark_matter
		icout=icout+1
		mb%dsp(icout)%p=>mb%nakedHe
		! icout=icout+1
		! mb%dsp(icout)%p=>mb%bstar
		! icout=icout+1
		! mb%dsp(icout)%p=>mb%bbh

		!print*, "mass bin 2"
		do i=1, n_tot_comp
			call mb%dsp(i)%p%init_e(df_coe_bins, dstr_bins,emin, emax, jmin, jmax, r0_cl, eb_type, jb_type)
		end do 

	end subroutine

	
    subroutine write_info_mass_bin(mb, funit)
        implicit none
        class(mass_bins)::mb
        integer funit, i

        call mb%dc%write_grid(funit)        
        write(funit) mb%mc, mb%m1, mb%m2, mb%df_coe_bins, mb%dstr_bins
        write(funit) mb%emin, mb%emax, mb%jmin, mb%jmax, &
             mb%mtot, mb%v0, mb%n0, mb%r0_cl
		!write(funit) mb%star%barj
		write(funit) mb%all%n,mb%star%n, mb%sbh%n, mb%wd%n, mb%ns%n, mb%bd%n, mb%rg%n, mb%nakedHe%n
		! write(funit) mb%bstar%n, mb%bbh%n	 
        write(funit) mb%all%barge, mb%all%fden,mb%all%fden_simu, mb%all%asymp
		write(funit) mb%all%fmden, mb%all%fma, mb%all%fna, mb%all%fma_simu, mb%all%fna_simu
		write(funit) mb%all%barge_ir, mb%all%spt_rho_rmin
		!call mb%all%barge_ir%print("barge_ir")
		do i=1, n_tot_comp
			write(funit) mb%dsp(i)%p%n, mb%dsp(i)%p%spt_rho_rmin
			!print*, "mb%dsp(i)%p%n=",mb%dsp(i)%p%n
			if(mb%dsp(i)%p%n>0)then
				write(funit) mb%dsp(i)%p%fden, mb%dsp(i)%p%fden_simu, mb%dsp(i)%p%barge,mb%dsp(i)%p%asymp
				write(funit) mb%dsp(i)%p%barge_ir
				!call mb%dsp(i)%p%barge_ir%print("i barge_ir")
				write(funit) mb%dsp(i)%p%fmden, mb%dsp(i)%p%fma, mb%dsp(i)%p%fna, mb%dsp(i)%p%fma_simu,&
					 mb%dsp(i)%p%fna_simu
			end if
		end do		

    end subroutine
    subroutine read_info_mass_bin(mb, funit)
        implicit none
        class(mass_bins)::mb
        integer funit,i
		
        call mb%dc%read_grid(funit)        
		!print*, "grid"
        read(funit) mb%mc, mb%m1, mb%m2, mb%df_coe_bins, mb%dstr_bins
        read(funit) mb%emin, mb%emax, mb%jmin, mb%jmax, &
            mb%mtot, mb%v0, mb%n0, mb%r0_cl
		!read(funit) mb%star%barj
		read(funit) mb%all%n,mb%star%n, mb%sbh%n, mb%wd%n, mb%ns%n, mb%bd%n,mb%rg%n,  mb%nakedHe%n
		! read(funit) mb%bstar%n, mb%bbh%n	 
		read(funit) mb%all%barge, mb%all%fden,mb%all%fden_simu, mb%all%asymp
		read(funit) mb%all%fmden, mb%all%fma, mb%all%fna, mb%all%fma_simu, mb%all%fna_simu
		read(funit) mb%all%barge_ir, mb%all%spt_rho_rmin
		!call mb%all%barge_ir%print("barge_ir")
		mb%all%spt_rho_rmax=0
		do i=1, n_tot_comp
			read(funit) mb%dsp(i)%p%n, mb%dsp(i)%p%spt_rho_rmin
			mb%dsp(i)%p%spt_rho_rmax=0
			!print*, "mb%dsp(i)%p%n=",mb%dsp(i)%p%n
			if(mb%dsp(i)%p%n>0)then
				read(funit) mb%dsp(i)%p%fden, mb%dsp(i)%p%fden_simu, mb%dsp(i)%p%barge,mb%dsp(i)%p%asymp
				read(funit) mb%dsp(i)%p%barge_ir
				!call mb%dsp(i)%p%barge_ir%print("i barge_ir")
				read(funit) mb%dsp(i)%p%fmden, mb%dsp(i)%p%fma, mb%dsp(i)%p%fna, mb%dsp(i)%p%fma_simu,&
						mb%dsp(i)%p%fna_simu
			end if
		end do		
        !read(funit) mb%all%barge, mb%all%fden, mb%star%fden, mb%sbh%fden, &
        !    mb%bstar%fden, mb%bbh%fden, mb%ns%fden, mb%wd%fden
        !read(funit) mb%all%fden_simu, mb%star%fden_simu, mb%sbh%fden_simu, &
        !    mb%bstar%fden_simu, mb%bbh%fden_simu, mb%ns%fden_simu, mb%wd%fden_simu
        !read(funit) mb%star%barge, mb%sbh%barge, mb%bbh%barge, mb%ns%barge, &
		!	mb%wd%barge
        !read(funit) mb%all%asymp,mb%star%asymp,mb%bstar%asymp,mb%sbh%asymp, &
        !    mb%bbh%asymp, mb%ns%asymp, mb%wd%asymp

    end subroutine
	subroutine get_dsp_idx_by_type(star_type,dsp_idx)
		! use md_sample,only:bytype_msb, bytype_bhb
		implicit none
		integer star_type,dsp_idx
		select case(star_type)
		case(star_type_ms)
			dsp_idx=1
		case(star_type_bh)
			dsp_idx=2
		case(star_type_ns)
			dsp_idx=3
		case(star_type_wd)
			dsp_idx=4
		case(star_type_bd)
			dsp_idx=5
		case(star_type_rg)
			dsp_idx=6
		case(star_type_dark_matter)
			dsp_idx=7
		case(star_type_nakedHe)
			dsp_idx=8
		! case(bytype_msb)
		! 	dsp_idx=9
		! case(bytype_bhb)
		! 	dsp_idx=10
		case default
			print*, "error! define obj_type", star_type
			stop
		end select
	end subroutine
	subroutine mb_get_fslope(mb,source)
		implicit none
		class(mass_bins)::mb
		integer i,source

		do i=1, n_tot_comp
			!if(mb%dsp(i)%p%n>0)then
			!	print*, "p=",i
			!	call mb%dsp(i)%p%fden%print("fden")
			!end if

			call dms_so_get_fslope(mb%dsp(i)%p,source)
		end do
		!print*, "2"
		!if(mb%all%n>0)then
		!	!print*, "all"
		!	call mb%all%fden%print("all fden")
		!end if

		! if(mb%all%n.eq.1536)then
		! 	print*, "type=",i
		! end if
		call dms_so_get_fslope(mb%all,source)
		! print*, "finished"
	end subroutine
	

	
	
end module