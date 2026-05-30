 
subroutine get_spp_fma(spp)
	use com_main_gw
	use md_star_pot
	implicit none 
	!type(diffuse_mspec)::dm
	type(star_pot_para)::spp
	integer i
	!logical::patch_outside
	!type(sts_fc_type)::frho_tot
	
	spp%fma_star%fx=0
	!if(spp%fma_star%nbin.ne.spp%frho_star%nbin)then
	!	print*, "error! dms%fma_star%nbin.ne.dms%all%all%fmden%nbin", spp%fma_star%nbin, spp%frho_star%nbin
	!	stop
	!end if
	call get_fma_tot(spp%fma_star,sample_logrmin,spp)
	!call spp%fma_star%print("fma_star")
	!stop
	!call get_beta_full_range(dms%fma_star,dms%fma_star%xmax, &
	!	M_r_within_max)
	call get_fma_tot_one(sample_logrmax,sample_logrmin,spp%M_r_within_max,spp)


	call get_fna_tot_one(spp)
	if(rid.eq.0)then
		print*, "N_within,M_within=",spp%N_r_within_max, spp%m_r_within_max
	end if
end subroutine
subroutine get_fna_tot_one(spp)
	use com_main_gw
	use md_star_pot
	implicit none
	type(star_pot_para)::spp
	integer i
	real(8) n_cum, nr
	n_cum=0
	do i=1, dms%n
		call dms%mb(i)%all%fna%get_value_s(sample_logrmax,nr)
		n_cum=n_cum+nr
		!print*, "i=",i
		!if(rid.eq.0)then
		!	call dms%mb(i)%all%fna%print('fna')
		!	print*, "nr=", nr, n_cum
		!end if
	end do

	spp%N_r_within_max=n_cum
	
	!print*, "n_cum=", n_cum, spp_new%M_r_within_max
	!read(*,*)
end subroutine
subroutine get_fma_tot_one(xb,rmin,fx, spp)
	use my_intgl
	use constant
	use com_sts_type
	!use model_basic,only:sample_logrmin
	use md_star_pot
	use model_basic,only:fr_funcs_int_acc_a,fr_funcs_int_acc_r
	implicit none
	!type(s1d_type)::frho
	type(star_pot_para)::spp
	integer idid
	real(8) rmin, yout,radius,xb, fx
	integer i
	!logical::patch_outside

	!rmin=frho%xb(1)-2; 
	!rmin=sample_logrmin

	yout=0
	call my_integral_acc(rmin,xb,yout,fr_funcs_int_acc_a,fr_funcs_int_acc_r,FCN, idid)
	fx=yout*4*pi*log(10d0)+10**(rmin*3)/3d0*spp%spt_rho_rmin*4*pi
	!if(yout<0)print*, "yout=",yout
	!print*, "xb,fx=",xb,fx,yout, 10**(rmin*3)/3d0*spt_rho_rmin*4*pi
contains 
	subroutine FCN(N,X,Y,F,IPAR,RPAR)
		implicit none
		integer n, ipar(100)
		real(8) x, y(n), f(n), rpar(100),ysp, rho_out

		call get_rho_full_range_spp(spp, x,rho_out)

		F(1)=rho_out* (10**X)**3
		!print*, "1:X,F=",X, F(1)
	end subroutine
end subroutine
subroutine get_fma_tot(fma,rmin, spp)
	use my_intgl
	use constant
	use com_sts_type
	use md_star_pot
	implicit none
	type(s1d_type)::fma
	type(star_pot_para)::spp
	integer idid
	real(8) rmin, yout,radius
	integer i
	
	do i=1, fma%nbin
		call get_fma_tot_one(fma%xb(i),rmin,fma%fx(i),spp)
		if(fma%fx(i)<0)fma%fx(i)=0
	enddo

end subroutine
 
subroutine get_spp_starpt_one(xb, fx, spp)
	use com_main_gw
	use md_star_pot
	implicit none
	!type(diffuse_mspec)::dm
	!type(s1d_type)::frho_star
	type(star_pot_para)::spp
	!call get_starpt(fphi_star%fx,fphi_star%nbin, frho_star%fx,frho_star%nbin, fphi_star%xb&
	!	, ctl%n0,ctl%v0, r0_cl)
	real(8) rmin, rmax, xb, fx
	real(8) logradius,radius,mtot
	integer i,idid,ier
	!logical::patch_outside

	call check_spp_data(spp,ier)
	if(ier<0)then
		print*, "error in get_dms_starpt_one"
		stop
	end if

	rmin=sample_logrmin; rmax=sample_logrmax
	!spt_rho_rmin=frho_star%fx(1)

	fx=0
	logradius=xb
	radius=10**xb
	!print*, "star starpot", spp%spt_rho_rmin
	
	call my_integral_acc(rmin,rmax,fx,fr_funcs_int_acc_a,fr_funcs_int_acc_r,FCN, idid)
	fx=fx*4*pi*log(10d0)+spp%spt_rho_rmin*10**(sample_logrmin*3)/3d0/radius*4*pi
	!print*, "rmin,rmax=",rmin,rmax, fx
contains 
	subroutine FCN(N,X,Y,F,IPAR,RPAR)
		implicit none
		integer n, ipar(100), j
		real(8) x, y(n), f(n), rpar(100),ysp

		call get_rho_full_range_spp(spp,x,ysp)

		if(X<logradius)then
			F(1)=ysp * (10**X)**3/radius
		else
			F(1)=ysp* (10**X)**2
		end if
		!print*, "X,F=",X, F(1)
	end subroutine
end subroutine
subroutine get_spp_starpt(spp)
	use com_main_gw
	use md_star_pot
	implicit none
	!type(diffuse_mspec)::dm
	type(star_pot_para)::spp
	!type(s1d_type)::fphi_star, frho_star
	!call get_starpt(fphi_star%fx,fphi_star%nbin, frho_star%fx,frho_star%nbin, fphi_star%xb&
	!	, ctl%n0,ctl%v0, r0_cl)
	real(8) logradius,radius,mtot,fout
	integer i,idid,ier
	
	call check_spp_data(spp,ier)
	if(ier<0) then
		print*, "error in get_spp_starpt"
		stop
	end if

	associate(fphi_star=>spp%fphi_star)
		do i=fphi_star%nbin, 1, -1
			call get_spp_starpt_one(fphi_star%xb(i), fphi_star%fx(i), spp)
			if(i<fphi_star%nbin)then  ! to make sure an increase function of phi
				if(fphi_star%fx(i)<fphi_star%fx(i+1))then
					fphi_star%fx(i)=fphi_star%fx(i+1)+1d-11
				end if
			end if
			!print*, "logr, phi=", logradius(i), phi_star(i)
		enddo
		select case(ctl%ebin_type)
		case(ebin_type_log)
			fphi_star%fx=log10(fphi_star%fx) !-phi_star(nbin)
		end select 
	end associate

	!call spp%fphi_star%print("1 phi")
	!do i=1, spp%fphi_star%nbin
	!	call get_spp_starpt_one_2(spp%fphi_star%xb(i),spp%fphi_star%fx(i),spp)
	!end do
	!call spp%fphi_star%print("2 phi")
	!read(*,*)
end subroutine
subroutine get_spt_phi_constants(spp)
	use md_star_pot
	use model_basic
	use com_sts_type
	use mpi_comu,only:rid
	implicit none
	!type(s1d_type)::frho_star
	type(star_pot_para)::spp
	real(8) fx,rmin,rmax, real_rmin
	integer idid
	!if(rid.eq.0)then
		!print*, "get constants:",frho_star%xmin,sample_logrmin,frho_star%xmax,sample_logrmax
	!end if
	!rmin=sample_logrmin;rmax=sample_logrmax
	rmin=spp%frho_star%xmin; rmax=spp%frho_star%xmax
	fx=0
	call my_integral_acc(rmin,rmax,fx,fr_funcs_int_acc_a,fr_funcs_int_acc_r,FCN, idid)
	spp%phi_r1r2_s=fx*4*pi*log(10d0)
	!real_rmin=10**rmin
	!if(spp%phi_r1r2_s+spp%spt_rho_rmin*(10**(rmin*2)*3-real_rmin**2)*4*pi/6d0<10**sample_logemax&
	!	.and.spp%mbh_dmless.eq.0)then
	!	if(rid.eq.0)then
	!		print*, "spp%phi_r1r2_s, cor,sample_emax=",spp%phi_r1r2_s, &
	!		10**sample_logemax-(10**(rmin*2)*3-real_rmin**2)*4*pi/6d0*spp%spt_rho_rmin,sample_emax
	!	end if
	!	spp%phi_r1r2_s=10**sample_logemax-(10**(rmin*2)*3-real_rmin**2)*4*pi/6d0*spp%spt_rho_rmin
	!end if
	fx=0
	!call spp%frho_star%print('frho')
	!call dms%all%all%barge_ir%print("barge")
	!print*, "rho_min=",spp%spt_rho_rmin
	call my_integral_acc(rmin,rmax,fx,fr_funcs_int_acc_a,fr_funcs_int_acc_r,FCN2, idid)
	spp%phi_r1r2_s2=fx*4*pi*log(10d0)

	spp%phi_star0=spp%phi_r1r2_s+0.5d0*10**(rmin*2)*spp%spt_rho_rmin
	!print*, "s2=",spp%phi_r1r2_s2
	!read(*,*)
	!if ((spt_rho_rmin*4*pi*10**(dms%fphi_star%xmin*3)/3d0+phi_r1r2_s2)/10**dms%fphi_star%xb(dms%fphi_star%nbin)&
	!	<dms%fphi_star%fx(dms%fphi_star%nbin))then
	!		print*, (spt_rho_rmin*4*pi*10**(dms%fphi_star%xmin*3)/3d0+phi_r1r2_s2)&
	!		/10**10**dms%fphi_star%xb(dms%fphi_star%nbin), dms%fphi_star%fx(dms%fphi_star%nbin)
	!		stop
	!end if
	!if(rid.eq.0)then
	!	print*, "phi_r1r2_s2=",phi_r1r2_s2
	!end if
	!print*, "phi_r1r2_s2=",spp%phi_r1r2_s2, ctl%plummer_model_mtot*(10**(rmax*3)/(10**(rmax*2)+1d0)**1.5&
	!	-10**(rmin*3)/(10**(rmin*2)+1d0)**1.5)
contains 
	subroutine FCN(N,X,Y,F,IPAR,RPAR)
		implicit none
		integer n, ipar(100), j
		real(8) x, y(n), f(n), rpar(100),ysp

		call get_rho_full_range_spp(spp,x,ysp)
		F(1)=ysp* (10**X)**2
	end subroutine
	subroutine FCN2(N,X,Y,F,IPAR,RPAR)
		implicit none
		integer n, ipar(100), j
		real(8) x, y(n), f(n), rpar(100),ysp

		call get_rho_full_range_spp(spp,x,ysp)
		F(1)=ysp* (10**X)**3
		!F(1)=ysp* (X)**2
	end subroutine
end subroutine
subroutine get_dms_alpha(dm,spp)
	use md_dms
	use md_star_pot
	implicit none
	!real(8) mbh
	type(diffuse_mspec)::dm
	type(star_pot_para)::spp

	dm%alpha_r%fx=10**spp%fphi_star%fx*(10**spp%fphi_star%xb)/spp%mbh_dmless
endsubroutine


subroutine get_jc_dmless(jcdm, rc, spp)
	use com_sts_type
	use md_star_pot
	implicit none
	type(s1d_type)::jcdm, rc
	type(star_pot_para)::spp
	real(8) r_c, jc_dmless
	integer i
	if(jcdm%nbin.ne.rc%nbin) then
		print*, "error! jdcm%nbin /= rc%nbin"
		stop
	end if
	do i=1, jcdm%nbin
		jcdm%fx(i)=jc_dmless(rc%fx(i),spp)
	end do
end subroutine

     
subroutine get_dehnen_rc(rc_theory, gamma, ra, mtot)
	use com_sts_type
	use constant
	implicit none
	real(8) gamma, ra, mtot,x
	type(s1d_type)::rc_theory
	integer i
	if(gamma.eq.1d0)then
		!print*, "finish it"
		!stop
		do i=1, rc_theory%nbin
			x=10**rc_theory%xb(i)			
			rc_theory%fx(i)=-ra+mtot/4d0/x+((-ra+mtot/4d0/x)**2-(ra**2-mtot/x*ra))**0.5
		end do
	end if
end subroutine
subroutine get_dehnen_jc(jc_theory, rc_theory,gamma, ra, mtot)
	use com_sts_type
	use constant
	implicit none
	real(8) gamma, ra, mtot,rc
	type(s1d_type)::rc_theory,jc_theory
	integer i
	if(gamma.eq.1d0)then
		!print*, "finish it"
		!stop
		do i=1, jc_theory%nbin
			rc=rc_theory%fx(i)			
			jc_theory%fx(i)=(rc*mtot*(rc/(rc+ra))**2)**0.5
		end do
	end if
end subroutine
subroutine get_dehnen_fe(fe_theory, gamma, ra, mtot)
	use com_sts_type
	use constant
	use model_basic,only:ctl
	use md_coeff
	implicit none
	type(s1d_type)::fe_theory
	real(8) gamma, ra, mtot,eta, x
	integer i

	eta=ra/mtot
	if(gamma.eq.0d0)then
		do i=1, fe_theory%nbin			
			select case(ctl%ebin_type)
			case(ebin_type_log)
				x=eta*10**fe_theory%xb(i)
			case(ebin_type_lin)
				x=eta*fe_theory%xb(i)
			end select
			fe_theory%fx(i)=2**0.5*3*eta**0.5/pi**1.5d0/ra**2*((2*x)**0.5*(3-4*x)/(1-2*x)-3*asinh((2*x/(1-2*x))**0.5) )
		end do
	end if
	if(gamma.eq.1d0)then
		!print*, "finish it"
		!stop
		do i=1, fe_theory%nbin			
			select case(ctl%ebin_type)
			case(ebin_type_log)
				x=eta*10**fe_theory%xb(i)
			case(ebin_type_lin)
				x=eta*fe_theory%xb(i)
			end select
			fe_theory%fx(i)=2**(-2d0)*eta**0.5/pi**1.5d0/ra**2/(1-x)**2.5*(3*asin(x**0.5)-(x*(1-x))**0.5*(3+2*x-24*x**2+16*x**3))
		end do
	end if
end subroutine
subroutine get_plummer_fe(fe_theory, ra, mtot)
	use com_sts_type
	use constant
	use model_basic,only:ctl
	use md_coeff
	implicit none
	type(s1d_type)::fe_theory
	real(8)  ra, mtot, x, n0, sigma0
	real(8) xmin
	integer i
	
	!n0=mtot/ra**3
	!sigma0=(mtot/ra)**0.5
	
	do i=1, fe_theory%nbin			
		select case(ctl%ebin_type)
		case(ebin_type_log)
			xmin=10**fe_theory%xb(fe_theory%nbin)
			x=10**fe_theory%xb(i)
		case(ebin_type_lin)
			xmin=fe_theory%xb(fe_theory%nbin)
			x=fe_theory%xb(i)
		end select
		fe_theory%fx(i)=24d0*2**0.5*ra**2*x**3.5d0/7d0/pi**3/mtot**4*(2*pi)**1.5d0
	end do

end subroutine
subroutine get_plummer_fma_s1d(fma_theory,mtot,ra)
	use com_sts_type
	use constant
	use model_basic,only:ctl
	implicit none
	type(s1d_type)::fma_theory
	real(8)  ra, mtot!, r
	integer i
	!mtot=ctl%plummer_model_mtot
	!ra=ctl%plummer_model_ra_critical
	do i=1, fma_theory%nbin
		!r=10**fe_theory%xb(i)
		call get_beta_plummer_outside(mtot,ra,fma_theory%xb(i),fma_theory%fx(i))
	end do

end subroutine
   
subroutine get_nx(nx,gx,rmax,  mbhin,spp)
	use com_sts_type
	use constant
	use my_intgl
	use model_basic,only:ctl
	use md_coeff,only:ebin_type_log,ebin_type_lin
	use md_star_pot
	implicit none
	type(s1d_type)::nx, gx,  fphi_star, barp
	type(s1d_type):: rmax
	type(star_pot_para)::spp

	integer i, idid
	real(8) logx,fout, logrmin, logr, mbhin
	select case(ctl%ebin_type)
	case(ebin_type_log)
		barp=nx
		call get_barp_xy(barp%xb, barp%fx,barp%nbin,rmax,mbhin,spp)
		nx%xb=barp%xb
		do i=1, nx%nbin
			logx=nx%xb(i)		
			!print*, "logrmin, max=", logrmin, rmax%fx(i)
			fout=0
			nx%fx(i)=2**1.5d0*pi**(0.5d0)*gx%fx(i)*10**nx%xb(i)*log(10d0)*barp%fx(i)! to log10 bin
			!print*, "nx%fx(i)=", nx%xb(i), nx%fx(i), gx%fx(i), fout
			!read(*,*)
			if(nx%fx(i)>0)then
				nx%fx(i)=log10(nx%fx(i))
			else
				nx%fx(i)=-100
			end if
		end do
		! call gx%print("gx")
		! call nx%print("nx")
		! read(*,*)
	case(ebin_type_lin)
		call get_barp_xy(log10(barp%xb), barp%fx,barp%nbin,rmax,mbhin,spp)
		do i=1, nx%nbin
			fout=0
			nx%fx(i)=2**1.5d0*pi**(0.5d0)*gx%fx(i)*barp%fx(i)
			!print*, "nx%fx(i)=", nx%xb(i), nx%fx(i), gx%fx(i), fout
			!read(*,*)
		end do
	end select

	
end subroutine
 
subroutine get_pd_dmless_mpi(spp,jc_dm_less,rp_dm,ra_dm,pd_dm_less)
	use com_sts_type
	use constant
	use model_basic,only:ctl
	use com_main_gw,only:dms
	use md_coeff
	use md_star_pot
	use mpi_comu
	implicit none
	type(star_pot_para)::spp
	type(s2d_type)::pd_dm_less, rp_dm, ra_dm
	type(s1d_type)::jc_dm_less
	real(8) p_EJ_dmless, enx, jm, jc, rp, ra,p_EJ_dmless_fast
	integer i,j,ierr
	real(8) t1, t2
	integer ibg,ied,nblock
	!call cpu_time(t1)
	nblock=pd_dm_less%nx/ctl%ntasks
	ibg=rid*nblock+1
	ied=(rid+1)*nblock

	do i=1, pd_dm_less%nx
		do j=ibg, ied
			select case(ctl%ebin_type)
			case(ebin_type_log)
				enx=10**pd_dm_less%xcenter(i)
			case(ebin_type_lin)
				enx=pd_dm_less%xcenter(i)
			end select
			jm=10**pd_dm_less%ycenter(j)
			jc=jc_dm_less%fx(i)
			rp=rp_dm%fxy(i,j)
			ra=ra_dm%fxy(i,j)
			!print*, "i,j,enx, jm,jc,ra,rp=",i,j,enx, jm,jc, ra,rp
			pd_dm_less%fxy(i,j)=p_EJ_dmless_fast(spp,enx,jm,jc,rp,ra)
			
		end do
	end do
	call mpi_barrier(mpi_comm_world,ierr)
	call collect_data_mpi_x(pd_dm_less%fxy,pd_dm_less%nx,ibg,ied,pd_dm_less%nx/ctl%ntasks,ctl%ntasks)
	!call cpu_time(t2)
	!print*, "mpi,period,t=",t2-t1,rid

end subroutine


subroutine get_pd_dmless(spp,jc_dm_less,rp_dm,ra_dm,pd_dm_less)
	use com_sts_type
	use com_sts_type
	use constant
	use model_basic,only:ctl
	use com_main_gw,only:dms
	use md_coeff
	use md_star_pot
	implicit none
	type(star_pot_para)::spp
	type(s2d_type)::pd_dm_less, rp_dm, ra_dm
	type(s1d_type)::jc_dm_less
	real(8) p_EJ_dmless, enx, jm, jc, rp, ra,p_EJ_dmless_fast
	integer i,j
	real(8) t1, t2
	call cpu_time(t1)
	do i=1, pd_dm_less%nx
		! print*, "i=",i
		do j=1, pd_dm_less%ny
			!enx=10**pd_dm_less%xcenter(pd_dm_less%nx)
			!j=42
			select case(ctl%ebin_type)
			case(ebin_type_log)
				enx=10**pd_dm_less%xcenter(i)
			case(ebin_type_lin)
				enx=pd_dm_less%xcenter(i)
			end select
			jm=10**pd_dm_less%ycenter(j)
			jc=jc_dm_less%fx(i)
			rp=rp_dm%fxy(i,j)
			ra=ra_dm%fxy(i,j)
			!print*, "i,j,enx, jm,jc,ra,rp=",i,j,enx, jm,jc, ra,rp
			pd_dm_less%fxy(i,j)=p_EJ_dmless_fast(spp,enx,jm,jc,rp,ra)
			!print*, "pd_dm_less%fxy(i,j)=", pd_dm_less%fxy(i,j)
			!pd_dm_less%fxy(i,j)=p_EJ_dmless_fast(spp,enx,jm,jc,rp,ra)
			!print*, "pd_dm_less%fxy(i,j)=", pd_dm_less%fxy(i,j)
			!read(*,*)
			!if(ctl%debug.ge.2)then
			!	print*, "i,j,enx, jm,jc,ra,rp=",i,j,enx, jm,jc, ra,rp
			!	print*, "pd_dm_less%fxy,pk=", pd_dm_less%fxy(i,j), 2*pi/(2*enx)**1.5
			!	!call dms%pd%print("pd")
			!	!read(*,*)
			!end if
	
	!		stop
		end do
	end do
	call cpu_time(t2)
	print*, "org period t=:", t2-t1, " s"
	!stop
end subroutine


subroutine get_rp_ra_dm(spp, jc_dm, rc, rmax, rp,ra )
	use com_main_gw
	implicit none
	type(s1d_type)::phi_star, jc_dm, rc, rmax
	type(s2d_type)::rp,ra
	type(star_pot_para)::spp
	real(8)r_root
	integer i,j, k
	real(8) enx, jm, jcdm_rc, logr, r,logrmax_tmp
	real(8) t1, t2
	call cpu_time(t1)
	do i=1, rp%nx
		select case(ctl%ebin_type)
		case(ebin_type_log)
			enx=10**rp%xcenter(i)
		case(ebin_type_lin)
			enx=rp%xcenter(i)
		end select
		!call cpu_time(t1)
		!do j=1,100
		!	print*, "j=",j
		call get_rmax_accurate(spp,rmax,rp%xcenter(i),logrmax_tmp)
		!end do
		!call cpu_time(t2)
		!print*, "rmax:", (t2-t1)/100d0
		!call cpu_time(t1)
		do j=1, rp%ny			
			jm=10**rp%ycenter(j)
			call get_rpra_dmless(spp, enx, jm, jc_dm%fx(i), log10(rc%fx(i)), logrmax_tmp, &
				rp%fxy(i,j),ra%fxy(i,j))
			!call get_rpra_dmless_old(spp, enx, jm, jc_dm%fx(i), log10(rc%fx(i)), logrmax_tmp, &
			!	rp%fxy(i,j),ra%fxy(i,j))
			!print*, "enx, jm, rp, ra=", enx, jm, rp%fxy(i,j),ra%fxy(i,j)
		end do
		!call cpu_time(t2)
		!print*, "rpra:", (t2-t1)/real(rp%ny)
		!read(*,*)
	end do
	call cpu_time(t2)
	print*, "org:t2-t1", t2-t1,rid 
end subroutine


subroutine get_rp_ra_dm_mpi(spp, jc_dm, rc, rmax, rp,ra )
	use com_main_gw
	implicit none
	type(s1d_type)::phi_star, jc_dm, rc, rmax
	type(s2d_type)::rp,ra
	type(star_pot_para)::spp
	real(8)r_root
	integer i,j, k
	real(8) enx, jm, jcdm_rc, logr, r,logrmax_tmp
	real(8) t1, t2
	integer ibg,ied,ierr,nblock
	!call cpu_time(t1)
	nblock=rp%nx/ctl%ntasks
	ibg=rid*nblock+1
	ied=(rid+1)*nblock
	do i=1, rp%nx 
		select case(ctl%ebin_type)
		case(ebin_type_log)
			enx=10**rp%xcenter(i)
		case(ebin_type_lin)
			enx=rp%xcenter(i)
		end select
		call get_rmax_accurate(spp,rmax,rp%xcenter(i),logrmax_tmp)
		!end do
		!call cpu_time(t2)
		!print*, "rmax:", (t2-t1)/100d0
		!call cpu_time(t1)
		do j=ibg, ied		
			jm=10**rp%ycenter(j)
			call get_rpra_dmless(spp, enx, jm, jc_dm%fx(i), log10(rc%fx(i)), logrmax_tmp, &
				rp%fxy(i,j),ra%fxy(i,j))
		end do
	end do
	if(rp%nx.ne.rp%ny) then
		print*, "error! rpnx.ne.rpny"
		stop
	end if
	call mpi_barrier(mpi_comm_world,ierr)

	call collect_data_mpi_x(rp%fxy,rp%nx,ibg,ied,rp%nx/ctl%ntasks,ctl%ntasks)
	call collect_data_mpi_x(ra%fxy,ra%nx,ibg,ied,ra%nx/ctl%ntasks,ctl%ntasks)
	!do i=1, ctl%ntasks
	!	if(i.eq.rid+1)then
	!		call rp%print("f:rp")
	!	end if
	!	call mpi_barrier(mpi_comm_world,ierr)
	!end do
	!call mpi_barrier(mpi_comm_world,ierr)
	!stop
	!print*, dms%rp%fxy(50,20:80),rid
    !print*, dms%rp%fxy(10,20:80),rid
	
	!call cpu_time(t2)
	!print*, "mpi:t2-t1", t2-t1,rid
end subroutine



subroutine get_frphi(phi_star,  rmax)
	use com_main_gw
	!use coms_sts_type
	!use model_basic:only:dms,ctl,ebin_type_log,ebin_type_lin
	implicit none 
	type(s1d_type)::phi_tot_sorted, rmax, phi_star, frho_star
	real(8) phi_out, fx(phi_star%nbin)
	integer i, nbin	
	rmax%xmin=sample_logemin
	rmax%xmax=sample_logemax
	do i=1, rmax%nbin
		if(spp_new%mbh_dmless.ne.0)then
			fx(i)=log10(10**phi_star%fx(i)+spp_new%mbh_dmless/10**phi_star%xb(i))
		else
			fx(i)=phi_star%fx(i)
		end if
	end do

	do i=1,rmax%nbin
		select case(ctl%ebin_type)
		case(ebin_type_log)
			rmax%xb(i)=fx(rmax%nbin-i+1)
			rmax%fx(i)=phi_star%xb(rmax%nbin-i+1)
		case(ebin_type_lin)
!			call get_rmax_accurate(phi_star,phi_tot_sorted,log10(rmax%xb(i)),rmax%fx(i))
			print*, "get_rmax:finish the code"
			stop
		end select
	end do
	do i=2, rmax%nbin
		if(rmax%fx(i)>rmax%fx(i-1))then
			if(rid.eq.0)then
				print*, "error! rmax is not declining!"
				print*, "i=",i
				call rmax%print("rmax")
				call phi_star%print("phi_star")
                stop    
			end if
		end if
	end do
	if(ctl%chattery.ge.3.and.rid.eq.0)then
        call rmax%print("fr_phi")
    end if
	!print*, "2"
	!call rmax%print("rmax")
	!call phi_tot_sorted%print("phi_tot_sorted")
	!read(*,*)
end subroutine

 
subroutine get_rmax_r_le_r1(spp, logex,logrmax)
	use model_basic
	use md_star_pot
	implicit none
	type(star_pot_para)::spp
	real(8) logex,logrmax
	real(8) px, qx, r1
	real(8) delta, term
	r1=10**sample_logrmin
	px=-(2*pi*spp%spt_rho_rmin*r1*r1+spp%phi_r1r2_s-10**logex)/(2*pi*spp%spt_rho_rmin)
	qx=-3*spp%mbh_dmless/4d0/pi/spp%spt_rho_rmin
	delta=qx**2+px**3
	term=(-qx+delta**0.5)**(1/3d0)+(-qx-delta**0.5)**(1/3d0)
	if(term>0)then
		logrmax=log10(term)
	else
		print*, "error! term<0", term
		stop
	end if
	!print*, "mbh/ex=",spp%mbh_dmless/10**logex
	!print*, "r1,cs=",r1, spp%phi_r1r2_s
	!print*, "delta, px,qx=",delta, px, qx
	!print*, "mbh,logex,logrmax=",spp%mbh_dmless,logex,logrmax, 10**logrmax
	!read(*,*)
end subroutine
subroutine get_rmax_accurate(spp,  fr_phi, logex, rmax)
	use model_basic
	use md_star_pot
	implicit none
	type(s1d_type):: fr_phi!,frho_star
	type(star_pot_para)::spp
	real(8) rmax, logex, ex
	real(8) rmax0, rtbis_yacc 
	real(8) xl,xh
	real(8) par(50), phi_star, width!,phi_star_tmp
	integer ier,niter
	integer::nmax_iter
	real(8)::require_acc
	!logical slient
	!if(ctl%debug.eq.1)then
	!	call fr_phi%print("fr_phi")
	!end if

	nmax_iter=10

	!call fr_phi%print("fr_phi")
	!read(*,*)
	ex=10**logex
	!select case(ctl%ebin_type)
	!case(ebin_type_log)
		if(logex.ge.log10emax_factor.and.spp%mbh_dmless.eq.0) then
			rmax=ctl%log10rmin_factor
			!print*, "get_rmax_accurate:ex, emax_factor=",ex,emax_factor
			return
		end if
		
		!if(logex.le.log10emin_factor) then
		!	rmax=ctl%log10rmax_factor
			!print*, "get_rmax_accurate:ex, emax_factor=",ex,emax_factor
		!	return
		!end if
		!if(logex.le.log10emin_factor)		
		!logex=8d0
		!print*, "true_log10emax_factor=",true_log10emax_factor
		if(logex.gt.sample_logemax.and.spp%mbh_dmless.ne.0.and.spp%spt_rho_rmin>0)then
			if(10**logex-2*pi*spp%spt_rho_rmin*10**(sample_logrmin*2)-spp%phi_r1r2_s>0)then
				call get_rmax_r_le_r1(spp,logex,rmax)
				return
			end if			
		end if
		call fr_phi%get_value_l(logex,rmax0)
		!if(logex>1e-4)then
		!	print*, "rmax0=",rmax0
		!end if
	!case(ebin_type_lin)

	!	if(ex>emax_factor) then
	!		!print*, "get_rmax_accurate:ex, emax_factor=",ex,emax_factor
	!		rmax=dms%logrmin
	!		return
	!	end if
	!	call fr_phi%get_value_l(ex,rmax0)
	!case default
	!	print*, "error in rmax 0", ctl%ebin_type
	!	stop
	!end select
	require_acc=1d-13

	if(spp%mbh_dmless.ne.0)then
		!call get_phi_star_full_range(fphi_star,logex,phi_star_tmp)
		!if(10**phi_star_tmp<mbh_dmless/10**rmax0*1d-5)then
		!	rmax=log10(mbh_dmless/(ex-10**phi_star_tmp))
		!	print*, "rmax=",rmax
		!	!return 
		!	ctl%debug=1
		!end if
		require_acc=1d-11
	end if

	!print*, "logex,rmax0=",logex, rmax0
	
	!if(mbh_dmless.eq.0.and.ex>fphi_star%fx(fphi_star%nbin))then
	!	rmax=fphi_star%xb(1)
	!	return
	!end if
	!print*, "start"
	width=0.1d0
	!if(mbh_dmless.ne.0)then
100		xl=rmax0-width
		xh=rmax0+width
	!else
	!	xl=max(rmax0-width,dms%logrmin)
	!	xh=min(rmax0+width,dms%logrmax)
	!	if(func(xl,par).eq.0)then
	!		rmax=xl
	!		return
	!	end if
	!	if(func(xh,par).eq.0)then
	!		rmax=xh
	!		return
	!	end if
	!end if
	
	!if(logex>1e-4)then
	!	print*, "xl,xh,rmax0=",xl,xh,rmax0, func(xl,par),func(xh,par)
	!end if
	rmax=rtbis_yacc(func,xl,xh,require_acc,par,niter,500,ier,.true.) 

	if(ier.eq.1)then
		!print*, "nmax_iter, rmax0, xl, xh=",nmax_iter, rmax0, xl, xh
		!call error_handle()
		if(nmax_iter>=0)then
			
			width=width/5d0
			require_acc=require_acc*2
			nmax_iter=nmax_iter-1
			rmax0=rmax

			!print*, "xl,xh,width,nmax_iter=",xl,xh,width,nmax_iter,rmax
			!print*, "fl,fh=",func(xl,par),func(xh,par)
			!print*, "xmid,dx, require_acc=",rmax0,par(50), require_acc
			!print*, "fx,fx2=",func(rmax0,par), func(rmax0+par(50),par), func(rmax0-par(50),par)
			!read(*,*)
			if(abs(par(50)/rmax0)>1d-15)then
				goto 100
			else
				return
			end if
		else
			print*, "get_rmax_accurate: reach maximum iteration,required acc not achieved", &
			"width, require_accr=", width, require_acc
			print*, "f(xl),f(xh)=",func(xl,par),func(xh,par)
			call error_handle()
			stop
		end if
	elseif(ier.eq.2)then
		if(nmax_iter>0)then
			width=width*5d0	
			nmax_iter=nmax_iter-1
			goto 100
		else
			print*, "get_rmax_accurate: reach maximum iteration, roots not bracked", &
			"width, require_accr=", width, require_acc
			print*, "f(xl),f(xh)=",func(xl,par),func(xh,par)
			call error_handle()
			stop
		end if

	elseif(ier.ne.0)then
		print*, "get_rmax_accurate problem"
		call error_handle()
		stop
	end if
	!print*, "ier=",ier
	!stop

	!if(ctl%debug.eq.1)then
	!	print*, "rmax=",rmax
	!	read(*,*)
	!end if
	!print*, "ex,rmax=",ex,logex,rmax,func(rmax,par)
	!block
	!	real(8) phi_out
	!	call get_phi_star_full_range(fphi_star,rmax,phi_out)
	!	print*, "phi_out=",phi_out,phi_out-logex,require_acc
	!end block
	!stop
	!read(*,*)
contains 
	real(8) function func(x,par)
		implicit none
		real(8) x, f, df, par(50), beta, phi_star_tmp

		!call get_starpt_at_r(frho_star,x,  phi_star) 
		call get_phi_star_full_range(spp,x,phi_star)
		!print*, "x, phi_star=",x,phi_star
		select case(ctl%ebin_type)
		case(ebin_type_log)
			phi_star_tmp=10**phi_star
		case(ebin_type_lin)
			phi_star_tmp=phi_star
		case default
			print*, "error in rmax"
			stop
		end select
		if(spp%mbh_dmless.eq.0)then
			func=ex-phi_star_tmp
			!print*, "func:x,ex,phi_star_tmp=",x,ex,phi_star_tmp
		else
			func=ex-phi_star_tmp-spp%mbh_dmless/10**x
		endif
		!call get_beta_full_range(fma, x, beta)
		!df=beta/10**(x)*log(10d0)
		!print*, "f,ex,x,phi=",func,ex,x,phi_star
	end function
	subroutine error_handle()
		implicit none
		type(s1d_type)::test
		integer i
		call test%init(xl,xh,300,sts_type_grid)
		call test%set_range()
		do i=1, 300
			test%fx(i)=func(test%xb(i),par)
			!read(*,*)
		end do
		print*, "rmax0,width=",rmax0,width
		print*, "best=",func(par(50),par)

		print*, "logex,xl, xh=", logex, xl, xh
		call test%print("test")
	    call fr_phi%print("fr_phi")
		call spp%fphi_star%print("fphi_star")
		call fr_phi%get_value_l(logex,rmax0)
		print*, "fr_phi(logex)=rmax0=",rmax0
		stop
	end subroutine
end subroutine
