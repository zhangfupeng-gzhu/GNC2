
module md_dms
	use md_mass_bins    
    use com_sts_type
   !integer,parameter::n_massbin=2
   type diffuse_mspec
        integer n
       type(mass_bins),allocatable:: mb(:)
       type(mass_bins):: all
       !real(8) weight_asym
       integer idx_ref
       !type(s1d_type)::barge0, barge_bh, barge_star
       type(diffuse_coeffient_type)::dc0
       !type(diffuse_coeffient_type)::dc0_bk   !dc in negative energy
       !type(s1d_type)::n_collid_11 ! expected collision rates
    !    type(s1d_type)::n_collid_gw ! expected collision rates
    !    type(s1d_type)::n_collid_gw_simu
    !    type(s1d_type)::favg_mass, favg_mass2   ! average mass at given r for all components
       type(s2d_type)::fvm    ! v_M(a,j) the orbital precession due to background masses
	   type(s2d_type)::pd, rp,ra
	   !type(s1d_type)::fphi_star ! in unit of sigma_h^2     x=log(r/r0_cl)
       type(s1d_type)::fphi_tot ! in unit of sigma_h^2     x=log(r/r0_cl)
	   !type(s1d_type)::fma_star ! in unit of nh r0_cl^3        x=log(r/r0_cl)
	   !type(s1d_type)::frho_star  ! in unit of 1Msun n_h    x=log(r/r0_cl)
	   type(s1d_type)::alpha_r, jc, rc, fr_phi!, frmax
       type(s1d_type)::jc_sample_erange
       type(s1d_type)::barp, frc_x
       type(s1d_type)::surface_den, cum_sur_den
       type(s1d_ird_type)::barp_ir, dlxb_ir
       type(s1d_type)::faniso
       !type(s1d_type)::fcri_emri
       !type(s1d_ird_type)::jc_ir

       integer df_coe_bins, dstr_bins_r,dstr_bins_e,jbin_type,ebin_type
       logical e_iregular
       real(8) emin, emax, jmin, jmax, mtot, v0, n0, r0_cl, logrmin, logrmax, phis0
       real(8) x_boundary!, logrmin, logrmax
       contains
    !   procedure:: set_binmass=>set_mass_bin_mass
       procedure:: set_diffuse_mspec
       procedure::init=>init_diffuse_mspec
    !    procedure::get_vm
   end type
   private::init_diffuse_mspec!, get_fxj0
   private::write_info_mass_bin,read_info_mass_bin!,print_norm_dms
contains
   subroutine init_diffuse_mspec(dm, n)
        implicit none
        integer n
        class(diffuse_mspec)::dm
        if(allocated(dm%mb))then
            deallocate(dm%mb)
        end if
        allocate(dm%mb(n))
        dm%n=n
    end subroutine
    subroutine init_stellar_obj_rtables(dm)
        implicit none
        class(diffuse_mspec)::dm
        integer i
        do i=1, dm%n
            call dm%mb(i)%init_r(dm%df_coe_bins, dm%dstr_bins_r ,dm%logrmin, &
            dm%logrmax, dm%jmin,dm%jmax, dm%mtot, dm%v0, dm%n0, dm%r0_cl,dm%jbin_type)
        !print*, dm%logrmin, dm%logrmax
        !call dm%mb(i)%dt%init(dstr_bins,dstr_bins,emin, emax, jmin, jmax,coeff_sts_type_dc)
        end do
        !print*, "3"
		call dm%all%init_r(dm%df_coe_bins, dm%dstr_bins_r , dm%logrmin, dm%logrmax, &
        dm%jmin,dm%jmax, dm%mtot, dm%v0, dm%n0, dm%r0_cl,dm%jbin_type)   
    end subroutine
    subroutine init_diffuse_mspec_rtables(dm)
        implicit none
        class(diffuse_mspec)::dm
        real(8) logrmin,logrmax
        integer dstr_bins, df_coe_bins

        logrmin=dm%logrmin; logrmax=dm%logrmax
        dstr_bins=dm%dstr_bins_r
        df_coe_bins=dm%df_coe_bins

		!call dm%fphi_star%init(logrmin,logrmax,dstr_bins,sts_type_dstr)
        call dm%fphi_tot%init( logrmin,logrmax,dstr_bins,sts_type_dstr)
		!call dm%fma_star%init( logrmin,logrmax,dstr_bins,sts_type_dstr)
		!call dm%frho_star%init(logrmin,logrmax,dstr_bins,sts_type_dstr)
        
		call dm%alpha_r%init( logrmin,logrmax,dstr_bins,sts_type_dstr)
		!call dm%beta_m%init( logrmin,logrmax,dstr_bins,sts_type_dstr)
		!call dm%fphi_star%set_range()
        call dm%fphi_tot%set_range()
		!call dm%fma_star%set_range()
		!call dm%frho_star%set_range()
		call dm%alpha_r%set_range()
		!call dm%beta_m%set_range()

        !logrmin=dm%logrmin; logrmax=dm%logrmax
        !call dm%n_collid_11%init(logrmin, logrmax, dm%dstr_bins, sts_type_dstr)
        ! call dm%n_collid_gw%init(logrmin, logrmax, dstr_bins, sts_type_dstr)
        ! call dm%n_collid_gw_simu%init(logrmin, logrmax, dstr_bins, sts_type_dstr)
        ! call dm%favg_mass%init(logrmin, logrmax, dstr_bins, sts_type_dstr)
        ! call dm%favg_mass2%init(logrmin, logrmax, dstr_bins, sts_type_dstr)
        call dm%fvm%init(dstr_bins, dstr_bins, logrmin, logrmax, dm%jmin, dm%jmax, sts_type_dstr)
        call dm%fvm%set_range()
        ! call dm%favg_mass%set_range()
        ! call dm%favg_mass2%set_range()
        !call dm%n_collid_11%set_range()
        ! call dm%n_collid_gw%set_range()
        ! call dm%n_collid_gw_simu%set_range()
    end subroutine
    subroutine init_dms_dc(dm)
        implicit none
        class(diffuse_mspec)::dm
        integer i
        call dm%dc0%init(dm%df_coe_bins,  dm%emin, dm%emax, dm%jmin,dm%jmax,&
            dm%ebin_type,dm%jbin_type)
        !call dm%dc0_bk%init(df_coe_bins,  emin, emax, jmin,jmax)
        !print*, "1"
        do i=1, dm%n
            !print*, "i=",i
            call dm%mb(i)%dc%init(dm%df_coe_bins,  dm%emin, dm%emax, dm%jmin,dm%jmax,&
                dm%ebin_type,dm%jbin_type)
            !print*, "df_coe_bins=", df_coe_bins
        end do

    end subroutine
    subroutine init_diffuse_mspec_etables(dm)
        implicit none
        class(diffuse_mspec)::dm
        real(8) logemin,logemax,jmin, jmax,tmin, tmax, logrmin, logrmax
        real(8) smin,smax
        integer jb_type,eb_type
        integer i

        jmin=dm%jmin; jmax=dm%jmax
        jb_type=dm%jbin_type
        eb_type=dm%ebin_type
        !print*, "logemin=",logemin, logemax
        !print*, "jmin,jmax=",jmin,jmax
        
        if(jb_type.eq.Jbin_type_log)then
            tmin=log10(jmin);tmax=log10(jmax)            
        end if
        if(eb_type.eq.ebin_type_log)then
            smin=log10(dm%emin);smax=log10(dm%emax)
        else
            smin=dm%emin;smax=dm%emax
        end if
        call dm%ra%init(dm%df_coe_bins,dm%df_coe_bins,smin,smax,tmin,tmax,sts_type_dstr)
        call dm%rp%init(dm%df_coe_bins,dm%df_coe_bins,smin,smax,tmin,tmax,sts_type_dstr)
        call dm%pd%init(dm%df_coe_bins,dm%df_coe_bins,smin,smax,tmin,tmax,sts_type_dstr)
        
        !call dm%jc%init(smin,smax,dm%df_coe_bins,sts_type_dstr)
        !if(dm%e_iregular)then
        !    call dm%barp_ir%init(smin,smax,dm%dstr_bins,sts_type_dstr)
        !else
        !    call dm%barp%init(smin,smax,dm%dstr_bins,sts_type_dstr)
        !    call dm%barp%set_range()
        !end if

		call dm%ra%set_range()
		call dm%rp%set_range()
		call dm%pd%set_range()
        !call dm%frmax%set_range()
        call dm%jc%init(smin,smax,dm%df_coe_bins,sts_type_dstr)
        call dm%rc%init(smin,smax,dm%df_coe_bins,sts_type_dstr)
        !==========================
        ! note that fr_phi use the bin size the same as phi
        call dm%fr_phi%init(smin,smax,dm%dstr_bins_r,sts_type_dstr)
        !==========================
		call dm%jc%set_range()
        call dm%rc%set_range()
        call dm%fr_phi%set_range()



        !print*, "4"
        !print*, "1"
        !print*, "ctl%nbins=",ctl%m_bins, ctl%bin_mass(1)
        !print*, "si=", ctl%asymptot_ini(1:8,1),ctl%bin_mass_m1(1)
        !print*, "2"

    end subroutine

    subroutine init_diffuse_mspec_gxtables(dm)
        implicit none
        class(diffuse_mspec)::dm
        real(8) jmin, jmax,tmin, tmax, logrmin, logrmax
        integer jb_type
        integer i
        
        !logemin=dm%logemin; logemax=dm%logemax
        jmin=dm%jmin; jmax=dm%jmax
        jb_type=dm%jbin_type

        do i=1, dm%n
            call dm%mb(i)%init_e(dm%df_coe_bins, dm%dstr_bins_e, dm%emin, dm%emax,&
             jmin, jmax, dm%mtot, dm%v0, dm%n0, dm%r0_cl, dm%ebin_type, jb_type)
        !print*, dm%logrmin, dm%logrmax
        !call dm%mb(i)%dt%init(dstr_bins,dstr_bins,emin, emax, jmin, jmax,coeff_sts_type_dc)
        end do
        !print*, "3"
		call dm%all%init_e(dm%df_coe_bins, dm%dstr_bins_e, dm%emin, dm%emax,  &
        jmin, jmax, dm%mtot, dm%v0, dm%n0, dm%r0_cl, dm%ebin_type, jb_type)       

    end subroutine
    
    subroutine set_diffuse_mspec(dm, df_coe_bins, dstr_bins_r, dstr_bins_e,  logrmin, logrmax, jmin, jmax, mtot, &
         v0, n0, r0_cl, xb, idx_ref,eb_type,jb_type,e_iregular)
        implicit none
        class(diffuse_mspec)::dm
        integer df_coe_bins, dstr_bins_r, dstr_bins_e, i, idx_ref,jb_type, eb_type
        real(8) jmin, jmax, mtot, v0, n0, r0_cl
        real(8) xb, tmin, tmax,logrmin, logrmax
        logical::e_iregular
        dm%e_iregular=e_iregular

        dm%idx_ref=idx_ref
        dm%logrmin=logrmin; dm%logrmax=logrmax
        !dm%logemin=logemin; dm%logemax=logemax
        dm%df_coe_bins=df_coe_bins; dm%dstr_bins_r=dstr_bins_r; dm%dstr_bins_e=dstr_bins_e; 
        
        dm%jmin=jmin; dm%jmax=jmax; dm%mtot=mtot; 
        dm%v0=v0; dm%n0=n0; dm%r0_cl=r0_cl
        !logrmin=log10(0.5d0*r0_cl/(10**logemax)); logrmax=log10(0.5d0*r0_cl/(10**logemin))
		!rrhmin=log10(0.5d0)-emax;rrhmax=log10(0.5d0)-emin
        !dm%logrmin=logrmin;dm%logrmax=logrmax
        dm%x_boundary=xb
        dm%jbin_type=jb_type; dm%ebin_type=eb_type
        
        !dm%favg_mass%nsam=1
        !dm%favg_mass2%nsam=1
    end subroutine

    subroutine set_mass_bin_mass_given(dm,masses, m1,m2, asym, n)
        implicit none
        class(diffuse_mspec) dm
        integer i, j, n
        real(8) masses(n), m1(n),m2(n), asym(n_tot_comp+1, n)
		!print*, "n=",n, n_tot_comp
        do i=1, n
			!print*, "i=",i
            dm%mb(i)%mc=masses(i)
            dm%mb(i)%m1=m1(i)
            dm%mb(i)%m2=m2(i)
            dm%mb(i)%all%asymp=asym(1, i)
			do j=1, n_tot_comp
				!print*, "j=",asym(j+1,i), asym(1,i)
				dm%mb(i)%dsp(j)%p%asymp=asym(j+1,i)*asym(1,i)
			end do
        end do
    end subroutine
	! subroutine get_fanso(dm)
	! 	implicit none
	! 	type(diffuse_mspec)::dm
	! 	real(8) e0
	! 	integer i
	! 	do i=1, dm%n
	! 		call mb_get_fanso(dm%mb(i))
	! 	end do
	! end subroutine
	subroutine get_fslope(dm,source)
		implicit none
		type(diffuse_mspec)::dm
		integer i,source
		do i=1, dm%n
            !print*, "i=",i
			call mb_get_fslope(dm%mb(i),source)
		end do

        call mb_get_fslope(dm%all,source)
	end subroutine
    subroutine set_mass_bin_mass_log(dm, mmin, mmax) ! in log 10 scale
        implicit none
        class(diffuse_mspec)::dm
        integer i, n
        real(8) mmin, mmax
        call set_range(dm%mb(:)%mc,dm%n,mmin,mmax,1)
        do i=1, dm%n
            dm%mb(i)%m1=10**(real(i-1)/real(dm%n)*(log10(mmax)-log10(mmin))+log10(mmin))
        end do
        do i=2, dm%n+1
            dm%mb(i-1)%m2=10**(real(i-1)/real(dm%n)*(log10(mmax)-log10(mmin))+log10(mmin))
        end do
    end subroutine
    !subroutine get_dc0_bk(dm)
    !    implicit none
    !    class(diffuse_mspec)::dm
    !    integer i, j
    !    real(8) kappa,sigma32, n0
    !    type(coeff_type)::cej
    !    external::fgx_bk
!
    !    kappa=(4*pi)**2*log(dm%mtot/1d0)  ! mass = 1 msun
    !    sigma32=(2*pi*dm%v0**2)**(-3/2d0)
    !    n0=dm%n0
    !    associate(dc=>dm%dc0_bk)
    !        do i=1, dc%nbin
    !            do j=1, dc%nbin
    !                call get_coeff_ej_sigma0_only(cej,10**dc%s2_de_110%xcenter(i), &
    !                    dc%s2_de_110%ycenter(j), dm%mtot, fgx_bk)
    !                dc%s2_de_110%fxy(i,j)=0d0
    !                dc%s2_de_0%fxy(i,j)=cej%e_0 *sigma32*n0*kappa
    !                if(cej%ee<0d0) cej%ee=abs(cej%ee)
    !                dc%s2_dee%fxy(i,j)=cej%ee*sigma32*n0*kappa
    !                dc%s2_dj_111%fxy(i,j)=0d0
    !                dc%s2_dj_rest%fxy(i,j)=cej%j_rest *sigma32*n0*kappa
    !                if(cej%jj<0d0) cej%jj=abs(cej%jj)
    !                dc%s2_djj%fxy(i,j)=cej%jj*sigma32*n0*kappa
    !                dc%s2_dej%fxy(i,j)=cej%ej*sigma32*n0*kappa                    
    !            end do
    !        end do
    !    end associate
    !end subroutine
    subroutine get_dc0(dm)
        implicit none
        integer n, nbin
        class(diffuse_mspec)::dm
    ! type(diffuse_coeffient_type):: dfs(n), dfs_tot
        integer i, j, k
        dm%dc0%s2_de_0%fxy=0
        dm%dc0%s2_dee%fxy=0
        dm%dc0%s2_dj_rest%fxy=0
        dm%dc0%s2_djj%fxy=0
        dm%dc0%s2_dej%fxy=0
        
        do i=1, dm%df_coe_bins
            do j=1, dm%df_coe_bins
                do k=1, dm%n
                    !dm%dc0%s2_de_110%fxy(i,j)=dm%dc0%s2_de_110%fxy(i,j)+dm%mb(k)%dc%s2_de_110%fxy(i,j)
                    dm%dc0%s2_de_0%fxy(i,j)=dm%dc0%s2_de_0%fxy(i,j)+dm%mb(k)%dc%s2_de_0%fxy(i,j)
                    dm%dc0%s2_dee%fxy(i,j)=dm%dc0%s2_dee%fxy(i,j)+dm%mb(k)%dc%s2_dee%fxy(i,j)
                    !dm%dc0%s2_dj_111%fxy(i,j)=dm%dc0%s2_dj_111%fxy(i,j)+dm%mb(k)%dc%s2_dj_111%fxy(i,j)
                    dm%dc0%s2_dj_rest%fxy(i,j)=dm%dc0%s2_dj_rest%fxy(i,j)+dm%mb(k)%dc%s2_dj_rest%fxy(i,j)
                    dm%dc0%s2_djj%fxy(i,j)=dm%dc0%s2_djj%fxy(i,j)+dm%mb(k)%dc%s2_djj%fxy(i,j)
                    dm%dc0%s2_dej%fxy(i,j)=dm%dc0%s2_dej%fxy(i,j)+dm%mb(k)%dc%s2_dej%fxy(i,j)
                end do
                !dm%dc0%s2_de_0%fxy(i,j)=dm%dc0%s2_de_0%fxy(i,j)+dm%dc0_bk%s2_de_0%fxy(i,j)
                !dm%dc0%s2_dee%fxy(i,j)=dm%dc0%s2_dee%fxy(i,j)+dm%dc0_bk%s2_dee%fxy(i,j)
                !!======
                !!dm%dc0%s2_de_110%fxy(i,j)=dm%dc0%s2_de_110%fxy(i,j)+dm%dc0_bk%s2_de_110%fxy(i,j)
                !!dm%dc0%s2_dj_111%fxy(i,j)=dm%dc0%s2_dj_111%fxy(i,j)+dm%dc0_bk%s2_dj_110%fxy(i,j)
                !!======
                !dm%dc0%s2_dj_rest%fxy(i,j)=dm%dc0%s2_dj_rest%fxy(i,j)+dm%dc0_bk%s2_dj_rest%fxy(i,j)
                !dm%dc0%s2_djj%fxy(i,j)=dm%dc0%s2_djj%fxy(i,j)+dm%dc0_bk%s2_djj%fxy(i,j)
                !dm%dc0%s2_dej%fxy(i,j)=dm%dc0%s2_dej%fxy(i,j)+dm%dc0_bk%s2_dej%fxy(i,j)
            end do
        end do        
    end subroutine     
     
end module



subroutine get_n_from_particle(mstar,n_star,&
    m1,  m2, idx, nsam)
    implicit none
    integer n_star
    real(8):: mstar(n_star)
    integer i, nsam, idx(n_star)
    real(8) jmax, m1,m2
    nsam=0
    if(n_star.eq.0) return
    do i=1, n_star
        !print*, ps_arr%sp(i)%m, m1, m2
        if(mstar(i).ge.m1.and.mstar(i).le.m2)then
            nsam=nsam+1
            idx(nsam)=i
        end if
    end do
    !print*, "nsam, n_star=",nsam, n_star
end subroutine


subroutine copy_dms(dm_source, dm_target)
    use md_dms
    implicit none
    type(diffuse_mspec)::dm_source
    type(diffuse_mspec),target:: dm_target
    integer i,j
    dm_target=dm_source
    do i=1, dm_target%n
        !do j=1, n_tot_comp
        associate(mb=>dm_target%mb(i))
            mb%dsp(1)%p=>mb%star
            mb%dsp(2)%p=>mb%sbh
            mb%dsp(3)%p=>mb%ns
            mb%dsp(4)%p=>mb%wd
            mb%dsp(5)%p=>mb%bd
            mb%dsp(6)%p=>mb%rg
            mb%dsp(7)%p=>mb%dark_matter
            mb%dsp(8)%p=>mb%nakedHe
            ! mb%dsp(9)%p=>mb%bstar
            ! mb%dsp(10)%p=>mb%bbh
            
        end associate
        associate(mb=>dm_target%all)
            mb%dsp(1)%p=>mb%star
            mb%dsp(2)%p=>mb%sbh
            mb%dsp(3)%p=>mb%ns
            mb%dsp(4)%p=>mb%wd
            mb%dsp(5)%p=>mb%bd
            mb%dsp(6)%p=>mb%rg
            mb%dsp(7)%p=>mb%dark_matter
            mb%dsp(8)%p=>mb%nakedHe
            ! mb%dsp(9)%p=>mb%bstar
            ! mb%dsp(10)%p=>mb%bbh
        end associate
        !end do
    end do
    

end subroutine