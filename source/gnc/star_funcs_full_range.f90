subroutine get_alpha_full_range(fphi,logr, alpha_out,spp)
	use com_sts_type
	use model_basic
	use md_coeff
	use md_star_pot
	implicit none
	type(s1d_type)::fphi
	real(8) logr, alpha_out,phi_out
	type(star_pot_para)::spp
	integer i
	
	call get_phi_star_full_range(spp,logr, phi_out)
	select case(ctl%ebin_type)
	case(ebin_type_log)
		alpha_out=10**logr*10**phi_out/spp%mbh_dmless 
	end select
end subroutine
subroutine get_phi_star_full_range(spp, logr, phi_out)
	!get full range of phi(star) assuming dehnen potential outside 
	use com_sts_type
	use model_basic
	use md_star_pot
	!use md_coeff
	implicit none
	real(8) logr, phi_out,phi_out_tmp,phi_out_0
	real(8) phi1,phi2
	!type(s1d_type)::phi_star
	type(star_pot_para)::spp
	integer n,i
	associate(phi_star=>spp%fphi_star)
		n=phi_star%nbin
		if(logr<=phi_star%xb(n).and.logr>=phi_star%xb(1))then
		!if(logr<=sample_logrmax.and.logr>=sample_logrmin)then
		!if(logr<=phi_star%xmax.and.logr>=phi_star%xmin)then
			call phi_star%get_value_s(logr,phi_out)
		!elseif(logr<phi_star%xb(1))then
		!elseif(logr<sample_logrmin)then
		elseif(logr>phi_star%xb(n).and.logr<phi_star%xmax)then
			phi1=phi_star%fx(n)
			phi2=log10((spp%spt_rho_rmin*4*pi*10**(phi_star%xmin*3)/3d0+spp%phi_r1r2_s2)/10**phi_star%xmax)
			phi_out=(phi2-phi1)/(phi_star%xmax-phi_star%xb(n))*(logr-phi_star%xb(n))+phi1
		elseif(logr<phi_star%xb(1).and.logr>phi_star%xmin)then
			phi1=log10(spp%phi_r1r2_s+4*pi*spp%spt_rho_rmin*10**(phi_star%xmin*2)/3d0)
			phi2=phi_star%fx(1)
			phi_out=(phi2-phi1)/(phi_star%xb(1)-phi_star%xmin)*(logr-phi_star%xmin)+phi1
		elseif(logr<=phi_star%xmin)then
			!call get_phi_dehnen(logr, phi_out_tmp)
			!call get_phi_dehnen(phi_star%xb(1),phi_out_0)
			!phi_out=log10(10**phi_out_tmp-10**phi_out_0+10**phi_star%fx(1))

			phi_out=log10(spp%phi_r1r2_s+4*pi*spp%spt_rho_rmin*(3*10**(phi_star%xmin*2)-10**(logr*2))/6d0)
			!print*, "logr<phi_star%xmin",logr, phi_star%fx(1), spt_rho_rmin, phi_out
		!elseif(logr>phi_star%xb(1))then
		!elseif(logr>sample_logrmax)then
		elseif(logr>=phi_star%xmax)then
			
			phi_out=(spp%spt_rho_rmin*4*pi*10**(phi_star%xmin*3)/3d0+spp%phi_r1r2_s2)/10**logr
			select case(ctl%ebin_type)
			case(ebin_type_log)
				phi_out=log10(phi_out)
			end select
		end if
	end associate
end subroutine
 
subroutine get_rho_full_range_log(frho,spt_rho_rmin,logr, rho_out)
	use com_sts_type
	use model_basic,only:ctl,ini_den_model_dehnen,ini_den_model_plummer
	implicit none
	type(s1d_type)::frho
	real(8) spt_rho_rmin,logr, rho_out, rho1,rho2 

	if(logr<=frho%xb(frho%nbin).and.logr>=frho%xb(1))then
		!here we should use linear interpolation
		!as if using spline interplotation,  as frho is in log scale
		! if it suddenly has very sharp transitions to small values such as -100,
		! it will cause large interpolation error

		call frho%get_value_l(logr,rho_out)
		!if(rho_out<0) rho_out=0
		rho_out=10**rho_out
	elseif(logr>frho%xmin.and.logr<frho%xb(1))then
		rho2=frho%fx(1)
		rho1=log10(spt_rho_rmin)
		rho_out=(rho2-rho1)/(frho%xb(1)-frho%xmin)*(logr-frho%xmin)+rho1
		rho_out=10**rho_out
	elseif(logr<frho%xmax.and.logr>frho%xb(frho%nbin))then
		rho2=frho%fx(frho%nbin)-4d0
		rho1=frho%fx(frho%nbin)
		rho_out=(rho2-rho1)/(frho%xmax-frho%xb(frho%nbin))*(logr-frho%xb(frho%nbin))+rho1
		rho_out=10**rho_out
	elseif(logr>=frho%xmax)then
		rho_out=0
		
	elseif(logr<=frho%xmin)then
		rho_out=spt_rho_rmin
	end if
end subroutine

subroutine get_rho_full_range(frho,spt_rho_rmin,logr, rho_out)
	use com_sts_type
	use model_basic,only:ctl,ini_den_model_dehnen,ini_den_model_plummer
	implicit none
	type(s1d_type)::frho
	real(8) spt_rho_rmin,logr, rho_out, rho1,rho2

	if(logr<=frho%xb(frho%nbin).and.logr>=frho%xb(1))then
		call frho%get_value_s(logr,rho_out)
		if(rho_out<0) rho_out=0
	elseif(logr>frho%xmin.and.logr<frho%xb(1))then
		rho2=frho%fx(1)
		rho1=spt_rho_rmin
		rho_out=(rho2-rho1)/(frho%xb(1)-frho%xmin)*(logr-frho%xmin)+rho1
	elseif(logr<frho%xmax.and.logr>frho%xb(frho%nbin))then
		rho2=0
		rho1=frho%fx(frho%nbin)
		rho_out=(rho2-rho1)/(frho%xmax-frho%xb(frho%nbin))*(logr-frho%xb(frho%nbin))+rho1
	elseif(logr>=frho%xmax)then
		rho_out=0d0
	elseif(logr<=frho%xmin)then
		rho_out=spt_rho_rmin
	end if
end subroutine

subroutine get_rho_full_range_spp(spp,logr, rho_out)
	use com_sts_type
	use model_basic 
	use md_star_pot
	implicit none
	!type(s1d_type)::frho
	real(8) logr, rho_out
	type(star_pot_para)::spp
	call get_rho_full_range(spp%frho_star,spp%spt_rho_rmin,logr,rho_out) 
end subroutine 
subroutine get_plummer_den(mtot,ra_scale, logr, rho_d)
	use com_main_gw
	implicit none
	real(8) logr, rho_d, r,ra_scale
	real(8) gamma, mtot
	real(8),parameter::c=3d0/4d0/pi
	!gamma=ctl%denhenmodel_gamma
	!mtot=ctl%denhenmodel_mtot
	!ra_scale=ctl%plummer_model_ra_critical
	r=10**logr
	rho_d=c*mtot*ra_scale**2/(r**2+ra_scale**2)**2.5d0
	!print*, "mtot, logr,rho=", mtot,logr,rho_d,r, ra_scale
end subroutine
subroutine get_rho_dehnen_outside(logr, rho_out)
    use model_basic,only:ctl
    implicit none
    real(8) fma_out,logr, rho_out,ysp
    integer j
!    rho_out=0
!    do j=1, ctl%m_bins
!        call get_dehnen_den(ctl%denhenmodel_mtot*ctl%asymptot_ini(1,j)*ctl%bin_mass(j), logr, ysp)
!        rho_out=rho_out+ysp
!    end do	
end subroutine 
subroutine get_ini_dehnen_dens(rho, Mtot, ra_scale,gamma)
	! get the dehnen density profile for initial condition Dehnen 1993, MNRAS, 265, 250
	! output rho%fx in unit of nh r0_cl^{-3}
	use com_sts_type
	use constant
	use model_basic
	implicit none
	type(s1d_type)::rho  ! with units of ctl%n0 (AU^-3) 
							! xb in units of log r/r0_cl
	integer i
	real(8) mtot ! total mass in units of m0_cl
	real(8) ra_scale ! scale radius in units of r0_cl
	real(8) gamma
	real(8) r
	!print*, "rho%nbin=", rho%nbin,gamma
	do i=1, rho%nbin
		r=10**rho%xb(i)
		call get_ini_dehnen_dens_func(rho%fx(i),r,mtot,ra_scale,gamma)
		!rho%fx(i)=(3-gamma)*mtot/4d0/pi*ra_scale/(r**gamma*(r+ra_scale)**(4-gamma))
	end do
end subroutine

subroutine get_ini_dehnen_dens_func(rho,r, Mtot, ra_scale,gamma)
	use constant
	implicit none
	real(8) mtot ! total mass in units of m0_cl
	real(8) ra_scale ! scale radius in units of r0_cl
	real(8) gamma
	real(8) r
	real(8) rho
	rho=(3-gamma)*mtot/4d0/pi*ra_scale/(r**gamma*(r+ra_scale)**(4-gamma))
end subroutine

subroutine get_ini_dehnen_fna(fna, ntot, ra_scale,gamma)
	use com_sts_type
	use constant
	use model_basic
	implicit none
	type(s1d_type)::fna  ! with units of ctl%n0 (AU^-3) 
							! xb in units of log r/r0_cl
	integer i
	real(8) ntot ! total number in units of m0_cl/msun
	real(8) ra_scale ! scale radius in units of r0_cl
	real(8) gamma
	real(8) r
	!print*, "rho%nbin=", rho%nbin,gamma
	do i=1, fna%nbin
		r=10**fna%xb(i)
		fna%fx(i)=ntot*(r/(r+ra_scale))**(3-gamma)
	end do
end subroutine

subroutine get_ini_plummer_fna(fna, ntot, ra_scale)
	! output rho%fx in unit of nh r0_cl^{-3}
	use com_sts_type
	use constant
	use model_basic
	implicit none
	type(s1d_type)::fna
	integer i
	real(8) ntot ! total number in units of m0_cl/1msun
	real(8) ra_scale ! scale radius in units of r0_cl
	real(8) gamma
	real(8) r
	!print*, "rho%nbin=", rho%nbin,gamma
	do i=1, fna%nbin
		call get_plummer_fna(ntot,ra_scale,fna%xb(i),fna%fx(i))
	end do
end subroutine


subroutine get_plummer_fna(ntot, ra_scale,logr, fna)
	use com_main_gw
	implicit none
	real(8) logr, fna, r
	real(8) gamma, ntot, ra_scale
	!real(8),parameter::c=3d0/4d0/pi
	!gamma=ctl%denhenmodel_gamma
	!mtot=ctl%denhenmodel_mtot
	!ra_scale=ctl%plummer_model_ra_critical
	r=10**logr
	fna=ntot*r**3/(r**2+ra_scale**2)**1.5d0
	!print*, "mtot, logr,rho=", mtot,logr,rho_d,r, ra_scale
end subroutine

subroutine get_ini_plummer_dens(rho, Mtot, ra_scale)
	! output rho%fx in unit of nh r0_cl^{-3}
	use com_sts_type
	use constant
	use model_basic
	implicit none
	type(s1d_type)::rho  ! with units of ctl%n0 (AU^-3) 
							! xb in units of log r/r0_cl
	integer i
	real(8) mtot ! total mass in units of m0_cl
	real(8) ra_scale ! scale radius in units of r0_cl
	real(8) gamma
	real(8) r
	!print*, "rho%nbin=", rho%nbin,gamma
	do i=1, rho%nbin
		call get_plummer_den(mtot,ra_scale,rho%xb(i),rho%fx(i))
	end do
end subroutine
subroutine get_beta_full_range(spp,logr, beta_out)
	use com_sts_type
	use model_basic 
	use md_star_pot
	implicit none
	type(s1d_type)::fma
	real(8) logr, beta_out
	real(8) fma_out,fma2
	type(star_pot_para)::spp

	associate(fma=>spp%fma_star)
		if(logr<=fma%xb(fma%nbin).and.logr>=fma%xb(1))then
		!if(logr<=sample_logrmax.and.logr>=sample_logrmin)then
		!if(logr<=fma%xmax.and.logr>=fma%xmin)then
			call fma%get_value_s(logr, fma_out)
			beta_out=fma_out 
			if(beta_out<0) beta_out=0
		elseif(logr>fma%xb(fma%nbin).and.logr<fma%xmax)then
			beta_out=(spp%M_r_within_max-fma%fx(fma%nbin))/(fma%xmax-fma%xb(fma%nbin))*&
				(logr-fma%xb(fma%nbin))+fma%fx(fma%nbin)
		elseif(logr<fma%xb(1).and.logr>fma%xmin)then		
			fma2=4*pi/3d0*spp%spt_rho_rmin*10**(fma%xmin*3)
			beta_out=(fma%fx(1)-fma2)/(fma%xb(1)-fma%xmin)*(logr-fma%xmin)+fma2
		!elseif(logr>fma%xb(fma%nbin))then
		!elseif(logr>sample_logrmax)then
		elseif(logr>=fma%xmax)then
			beta_out=spp%M_r_within_max
		elseif(logr<=fma%xmin)then
			beta_out=4*pi/3d0*spp%spt_rho_rmin*10**(logr*3)
		end if
	end associate
end subroutine
subroutine get_gx_full_range(gx,xi_i,enx,gx_out)
    use com_main_gw
    implicit none
    type(s1d_type)::gx
    real(8) enx,gx_out, xi_i
	!logical,save::first=.true.
	
    if(enx>=gx%xmin.and.enx<=gx%xmax)then
        call gx%get_value_s(enx,gx_out)
    else
	!	select case(ctl%ini_den_model)
	!	case(ini_den_model_dehnen)
    !    	call get_gx_dehnen_outside(logenx,xi_i,mbh_dmless,&
    !        	ctl%denhenmodel_gamma,ctl%denhenmodel_ra_critical,gx_out)
	!	case(ini_den_model_plummer)
	!		call get_gx_plummer_outside(logenx,xi_i,&
    !        	ctl%plummer_model_ra_critical,gx_out)
	!	case default
	!		print*, "error in gx full range", ctl%ini_den_model
	!		stop
	!	end select
		gx_out=0
    end if
	!write(*,fmt="(A20,2E30.23)"), "logenx=",logenx, gx_out
end subroutine
subroutine get_gx_full_range_ir(gx,enx,gx_out)
    use com_sts_type
    use model_basic
    implicit none
    type(s1d_ird_type)::gx
    real(8) enx,gx_out
	!logical,save::first=.true.
	!if(fresh_common_gx)then
	!	common_gx=gx
	!	call get_none_zero_s1d(gx,common_gx)
	!	fresh_common_gx=.false.
	!	call gx%print("gx")
	!	call common_gx%print("common_gx")
	!	read(*,*)
	!end if
	!print*, "enx=",enx
	
    if(enx>=gx%xmin.and.enx<=gx%xmax)then
        call gx%get_value_s(enx,gx_out)
		!call gx%print("gx")
		!print*, "enx,gx_out=",enx,gx_out
		!read(*,*)
		!gx_out=10**gx_out
		
		!call gx%get_value_s(enx,gx_out)
		!if(gx_out<0) gx_out=0
    else
	!	select case(ctl%ini_den_model)
	!	case(ini_den_model_dehnen)
    !    	call get_gx_dehnen_outside(logenx,xi_i,mbh_dmless,&
    !        	ctl%denhenmodel_gamma,ctl%denhenmodel_ra_critical,gx_out)
	!	case(ini_den_model_plummer)
	!		call get_gx_plummer_outside(logenx,xi_i,&
    !        	ctl%plummer_model_ra_critical,gx_out)
	!	case default
	!		print*, "error in gx full range", ctl%ini_den_model
	!		stop
	!	end select
		gx_out=0
		!gx_out=-100
    end if

end subroutine
subroutine get_gx_full_range_ir_log(gx,enx,gx_out)
    use com_sts_type
    use model_basic
    implicit none
    type(s1d_ird_type)::gx
    real(8) enx,gx_out
	!logical,save::first=.true.
	!if(fresh_common_gx)then
	!	common_gx=gx
	!	call get_none_zero_s1d(gx,common_gx)
	!	fresh_common_gx=.false.
	!	call gx%print("gx")
	!	call common_gx%print("common_gx")
	!	read(*,*)
	!end if
	!print*, "enx=",enx
	
    if(enx>=gx%xmin.and.enx<=gx%xmax)then
        call gx%get_value_s(enx,gx_out)

		!call gx%print("gx")
		!print*, "enx,gx_out=",enx,gx_out
		!read(*,*)
		gx_out=10**gx_out
		
		!call gx%get_value_s(enx,gx_out)
    else
	!	select case(ctl%ini_den_model)
	!	case(ini_den_model_dehnen)
    !    	call get_gx_dehnen_outside(logenx,xi_i,mbh_dmless,&
    !        	ctl%denhenmodel_gamma,ctl%denhenmodel_ra_critical,gx_out)
	!	case(ini_den_model_plummer)
	!		call get_gx_plummer_outside(logenx,xi_i,&
    !        	ctl%plummer_model_ra_critical,gx_out)
	!	case default
	!		print*, "error in gx full range", ctl%ini_den_model
	!		stop
	!	end select
		gx_out=0
    end if

end subroutine
subroutine get_beta_dehnen_outside(logr, beta_out)
    use model_basic,only:ctl
    implicit none
    real(8) fma_out,logr, beta_out,ysp
    integer i

!   fma_out=0
!   do i=1, ctl%m_bins
!       call get_dehnen_fma(ctl%denhenmodel_mtot*ctl%asymptot_ini(1,i)*ctl%bin_mass(i), logr, ysp)
!       !print*, "i,logr,ysp=",i,logr,ysp
!       fma_out=fma_out+ysp
!   end do	
!   beta_out=fma_out 
end subroutine

subroutine get_beta_plummer_outside(mtot,ra,logr, beta_out)
    use model_basic,only:ctl
    implicit none
    real(8) fma_out,logr, beta_out,ysp,mtot,ra
    integer i

    fma_out=0
    do i=1, ctl%m_bins
        call get_plummer_fma(mtot*ctl%asymptot_ini(1,i)*ctl%bin_mass(i),ra, logr, ysp)
        !print*, "i,logr,ysp=",i,logr,ysp
        fma_out=fma_out+ysp
    end do	
    beta_out=fma_out 
end subroutine

!subroutine get_dehnen_fma(mtot,logr,fmaout)
!	use constant
!	use model_basic,only:ctl
!	implicit none
!	real(8) mtot, r, logr, fmaout
!	r=10**logr
!	fmaout=mtot*(r/(r+ctl%denhenmodel_ra_critical))**(3-ctl%denhenmodel_gamma)
!end subroutine
subroutine get_plummer_fma(mtot,ra,logr,fmaout)
	use constant
	use model_basic,only:ctl
	implicit none
	real(8) mtot, r, logr, fmaout,ra
	r=10**logr
	fmaout=mtot*r**3/(r**2+ra**2)**1.5d0
end subroutine
subroutine get_dehnen_pot(phi_star, r, mtot, ra, gamma)
	implicit none
	real(8) phi_star, mtot, ra, gamma, r

	if(gamma.ne.2d0)then
		phi_star=mtot/ra/(2-gamma)*(1-(r/(r+ra))**(2-gamma))
	else
		phi_star=-mtot/ra*(log(r/(r+ra)))
	end if
end subroutine

subroutine get_plummer_pot(phi_star, r, mtot, ra)
	implicit none
	real(8) phi_star, mtot, ra, gamma, r

	phi_star=mtot/(r**2+ra**2)**0.5

end subroutine
subroutine get_ini_plummer_pot(pot,plummer)
	! get the plummer potential for initial condition
	! in unit of Mbh/r0_cl =sigma_h^2
	use com_sts_type
	use constant
	use model_basic,only:ctl,type_mplummer
	implicit none
	type(s1d_type)::pot
	type(type_mplummer)::plummer
	real(8) mtot ! total mass in units of MBH
	real(8) ra_scale ! scale radius in units of r0_cl
	real(8) gamma
	real(8) r
	integer i
	do i=1, pot%nbin
		r=10**pot%xb(i)
		call get_plummer_pot(pot%fx(i), r,  &
		plummer%mtot, plummer%ra_crit)
	end do
end subroutine

subroutine get_ini_dehnen_pot(pot,dehnen)
	! get the dehnen potential for initial condition Dehnen 1993, MNRAS, 265, 250
	! the sign of potential is inversed compared to Dehnen
	! in unit of Mbh/r0_cl =sigma_h^2
	use com_sts_type
	use constant
	use model_basic,only:ctl,type_mdehnen
	implicit none
	type(s1d_type)::pot
	type(type_mdehnen)::dehnen
	real(8) mtot ! total mass in units of MBH
	real(8) ra_scale ! scale radius in units of r0_cl
	real(8) gamma
	real(8) r
	integer i
	do i=1, pot%nbin
		r=10**pot%xb(i)
		call get_dehnen_pot(pot%fx(i), r,  &
		dehnen%mtot, dehnen%ra_crit, dehnen%gamma)
	end do
end subroutine
 