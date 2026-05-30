 
subroutine get_rpra_dmless_fast(spp,enx,jm,Jc_dm,logrc,logrmax, rp,ra)
	use ieee_arithmetic
	use com_sts_type
    use model_basic,only:ctl,dms,jmin_value,sample_logrmin
	use md_coeff
	use md_star_pot
	use constant
	use MPI_comu,only:rid
	implicit none
	!logical require_accurate
	integer::nb,i,minlocate(1)
	!type(s1d_type)::phi_star
	type(star_pot_para)::spp
	real(8)  enx, jm, rp, ra, logrmin,logrmax
	real(8) par(50),r1,root_finding_1,r2
	real(8) jc_dm, logrc, ac, ec,e0,pxmin,vr_tmp
	real(8) phi_star_tmp,phi0, edge0, edge1, jph2,rtbis_yacc, rho_tmp,mu
	integer ier,niter
	logical,parameter::slient=.true.
	!real(8) xb(phi%nbin), yb(phi%nbin)

!	respect to the Keplerian case
!	ac=1d0/2d0/enx
!	ec=(1-(jm*jc_dm)**2)**0.5
	!if(jm>=1d0-1d-5) then 
	!	rp=10**logrc; ra=10**logrc
	!	return
	!end if
	if(jm>=1d0-1d-5) then 
		call get_rho_full_range_spp(spp,logrc,rho_tmp)
		mu=(1-jm**2)**0.5/(1+4*pi*rho_tmp*10**(logrc*4)/Jc_dm**2)**0.5
		!print*, "the:rp,ra=",(1-mu)*10**logrc,(1+mu)*10**logrc
		rp=(1-mu)*10**logrc; ra=(1+mu)*10**logrc
		!rp=10**logrc; ra=10**logrc
		return
	end if

	if(logrmax<logrc)then
		print*, "get_rpra_dmless:error! logrmax<logrc,rid",logrmax,logrc,rid
		print*, "enx,jm,jc_dm=",enx,jm, jc_dm
		stop
	end if
	!logrmin=ctl%log10rmin_factor+log10(jmin_value)*2.5
	jph2=(jm*jc_dm)**2
	
	if(func(sample_logrmin,par)<0.or.enx<spp%phi_r1r2_s+2*pi*spp%spt_rho_rmin*10**(sample_logrmin*2))then
		logrmin=logrc+log10(jmin_value)*2-2
		edge0=func(logrmin,par)
		if(edge0<1d-8.and.edge0.ge.-1d-12) then
			r1=logrmin
			goto 100
		end if
		edge1=func(logrc,par)
		if(edge1<1d-8.and.edge1.ge.-1d-12)then
			r1=logrc
			goto 100
		end if
		if(func(logrmin,par)*func(logrc,par).gt.0)then
			print*, "in get_rpra_dmless_fast"
			print*, "no root between rmin and rc"
			call error_handle()
		end if
		!print*, "logrc, logrmax=",logrc, logrmax+0.1
		!print*, func(logrc,par),func(logrmax+0.1,par)
		r1=root_finding_1(func,logrmin,logrc,1d-4,10,par,ier,slient)
		!print*, "r1=",r1
		!r1=root_finding_1(func,logrmin,logrc,1d-4,10,par,ier,slient)
		!print*, "r1=",r1
		!r1=rpra_finding_fast(func,logrmin,logrc,ier)
		
		call root_finding_error()

	!print*, "2"
100		rp=10**r1
		if(ieee_is_nan(rp))then
			print*, "error! rp is NaN", logrmin, logrc, r1, rp
			
			stop
		end if
	else
		if(spp%spt_rho_rmin>0d0)then
			call get_rp_from_qp(spp, enx, jph2, rp)
		else
			e0=spp%phi_r1r2_s-enx
			rp=-spp%mbh_dmless/(2d0*e0)*&
			(1-(1+(jm*jc_dm)**2*2*e0/spp%mbh_dmless**2)**0.5d0)
		end if
	end if


	edge0=func(logrc,par)
	if(edge0<1d-8.and.edge0>=-1d-12) then
		r2=logrc
		goto 200
	end if
	edge1=func(logrmax,par)
	if(edge1<1d-8.and.edge1>=-1d-12)then
		r2=logrmax
		goto 200
	end if
	if(func(logrc,par)*func(logrmax+0.5,par).gt.0)then
		print*, "in get_rpra_dmless_fast"
		print*, "no root between rc and rmax"
		call error_handle()
	end if

	!r2=rtbis_yacc(func,logrc,logrmax+0.5,1d-8,par,ier,slient)
	r2=root_finding_1(func,logrc,logrmax+0.5,1d-4,10,par,ier,slient)
	!print*, "r2=",r2
	!r2=root_finding_1(func,logrc,logrmax+0.5,1d-4,10,par,ier,slient)
	!print*, "r2=",r2
	!read(*,*)
	!r2=rpra_finding_fast(func,logrc,logrmax+0.5,ier)
	if(ier.eq.1)then
		
		print*, "r2 root not found closest=", r2, func(r2,par), logrc, logrmax+0.5
		call error_handle()
	end if

200	ra=10**r2
	if(rp>ra)then
		print*, "get_rpra_dmless_fast:error! rp>ra",rp,ra
		print*, "assuming within r?",func(sample_logrmin,par)<0&
		.or.enx<spp%phi_r1r2_s+2*pi*spp%spt_rho_rmin*10**(sample_logrmin*2)
		stop
	end if
	!print*, "rp,ra=",rp,ra, func(log10(rp),par),func(log10(ra),par), jc_dm, jm
	!read(*,*)
	!if(func(r1,par)>0.001.or.func(r2,par)>0.001)then
	!	print*, "????"
	!	print*, "rp,ra=",rp,ra, func(log10(rp),par),func(log10(ra),par), jc_dm, jm
	!	read(*,*)
	!end if
	!read(*,*)
contains
	subroutine root_finding_error()
		if(ier.eq.1)then
			if(.not.slient)then
				select case(ctl%ebin_type)
				case(ebin_type_log)
					print*, "func(phi_star%xb(1),par)>0, 10**phi_star%fx(1)*2<enx, logrc>pxmin"
					print*, func(spp%fphi_star%xb(1),par)>0, 10**spp%fphi_star%fx(1)*2<enx, logrc>pxmin
					print*, "r1 root not found closest=", r1, func(r1,par)
				case(ebin_type_lin)
					print*, "func(phi_star%xb(1),par)>0, phi_star%fx(1)*2<enx, logrc>pxmin"
					print*, func(spp%fphi_star%xb(1),par)>0, spp%fphi_star%fx(1)*2<enx, logrc>pxmin
					print*, "r1 root not found closest=", r1, func(r1,par)
				end select
				call error_handle()
			else
				select case(ctl%ebin_type)
				case(ebin_type_log)
					print*, "func(phi_star%xb(1),par)>0, 10**phi_star%fx(1)*2<enx, logrc>pxmin"
					print*, func(spp%fphi_star%xb(1),par)>0, 10**spp%fphi_star%fx(1)*2<enx, logrc>pxmin
					print*, "r1 root not found closest,r1,logrmin, func(r1),ex,jm=", r1, logrmin, func(r1,par), enx,jm
				case(ebin_type_lin)
					print*, "func(phi_star%xb(1),par)>0, 10**phi_star%fx(1)*2<enx, logrc>pxmin"
					print*, func(spp%fphi_star%xb(1),par)>0, spp%fphi_star%fx(1)*2<enx, logrc>pxmin
					print*, "r1 root not found closest,r1,logrmin, func(r1),ex,jm=", r1, logrmin,func(r1,par), enx,jm
				end select
				call error_handle()
			end if
		!elseif(ier.eq.2)then
		!	print*, "root_finding_error:??",ier
		!	call error_handle()
		!	stop
		end if
	end subroutine
	real(8) function func(x,par)
		implicit none
		real(8) x, par(50)

		!call get_r_root(phi_star,rc,enx,1d0, jc_dm,func)
		!print*, "func=", func
		!stop
        call get_r_root(spp,x, enx, jph2, func)
	end function
	subroutine error_handle()
		implicit none
		type(s1d_type)::test
		print*, "logex,ex,jm=",log10(enx),enx,jm
		print*, "logmin, logrc=", logrmin, logrc
		print*, "fun(min), fun(rc)=", func(logrmin, par), func(logrc,par)
		print*, "rc, max=", logrc, logrmax
		print*, "fun(rc), fun(max)=", func(logrc, par), func(logrmax,par)
		!call dms%fr_phi%print("rmax")
		!call dms%rc%print("rc")
		call spp%fphi_star%print("fphi")
		call test%init(logrmin,logrmax,100,sts_type_grid)
		call test%set_range()
		do i=1, test%nbin
			test%fx(i)=func(test%xb(i), par)
			!print*, test%xb(i),par(1)
		end do
		call test%print("test")
		stop
	end subroutine
end subroutine
!subroutine rpra_finding_fast(func,spp,rmin,rmax,ier)
!	use com_main_gw
!	implicit none
!	type(star_pot_para)::spp
!	real(8),external::func
!	real(8) rmin,rmax
!	real(8) trmin,trmax
!	integer ier, nbin
!	integer ibg, ied, table_type, ibgmin,iedmax
!	real(8)  line_slope, line_c , y1, y2, x1, x2
!	real(8) f1,f2,par(50)
!	trmin=spp%fphi_star%xmin
!	trmax=spp%fphi_star%xmax
!	nbin=spp%fphi_star%nbin
!	table_type=spp%fphi_star%bin_type
!	call return_idx(rmin,trmin,trmax,nbin,ibgmin,table_type)
!	call return_idx(rmax,trmin,trmax,nbin,iedmax,table_type)
!	if(ibgmin.eq.iedmax)then
!		print*, "rmin,rmax=",rmin,rmax
!		print*, "ibg,ied=",ibg,ied
!		call spp%fphi_star%print("fphi_star")
!		read(*,*)
!	end if
!	f1=func(rmin,par)
!	f2=func(rmax,par)
!	if(f1*f2>0)then
!		print*, "error! roots not bracked", rmin,rmax, f1, f2
!		stop
!	end if
!	y1=10**spp%fphi_star%fx(ibg)
!	y2=10**spp%fphi_star%fx(ied)
!	x1=1d0/spp%fphi_star%xb(ibg)
!	x2=1d0/spp%fphi_star%xb(ied)
!	line_slope=(y1-y2)/(x1-x2)
!	line_c=y1-line_slope*x1
!
!end subroutine
subroutine get_rpra_dmless(spp,enx,jm,Jc_dm,logrc,logrmax, rp,ra)
	use ieee_arithmetic
	use com_sts_type
    use model_basic,only:ctl,dms,jmin_value,sample_logrmin
	use md_coeff
	use md_star_pot
	use constant
	use MPI_comu,only:rid
	implicit none
	integer::nb,i,minlocate(1)
	!type(s1d_type)::phi_star
	type(star_pot_para)::spp
	real(8)  enx, jm, rp, ra, logrmin,logrmax
	real(8) par(50),r1,root_finding_1,r2
	real(8) jc_dm, logrc, ac, ec,e0,pxmin,vr_tmp
	real(8) phi_star_tmp,phi0, edge0, edge1, jph2, mu, rho_tmp
	integer ier
	logical,parameter::slient=.true.
	!real(8) xb(phi%nbin), yb(phi%nbin)

!	respect to the Keplerian case
!	ac=1d0/2d0/enx
!	ec=(1-(jm*jc_dm)**2)**0.5
	if(jm>=1d0-1d-5) then 
		call get_rho_full_range_spp(spp,logrc,rho_tmp)
		mu=(1-jm**2)**0.5/(1+4*pi*rho_tmp*10**(logrc*4)/Jc_dm**2)**0.5
		!print*, "the:rp,ra=",(1-mu)*10**logrc,(1+mu)*10**logrc
		rp=(1-mu)*10**logrc; ra=(1+mu)*10**logrc
		!rp=10**logrc; ra=10**logrc
		return
	end if
	if(logrmax<logrc)then
		print*, "get_rpra_dmless:error! logrmax<logrc,rid",logrmax,logrc,rid
		print*, "enx,jm,jc_dm=",enx,jm, jc_dm
		stop
	end if
	!logrmin=ctl%log10rmin_factor+log10(jmin_value)*2.5
	jph2=(jm*jc_dm)**2


	pxmin=spp%fphi_star%xmin
	select case(ctl%ebin_type)
	case(ebin_type_log)
		!call get_phi_star_full_range(dms%fphi_star,pxmin,phi0)
		phi0=spp%phi_r1r2_s+4*pi*spp%spt_rho_rmin*10**(spp%fphi_star%xmin*2)/3d0
	case(ebin_type_lin)
		!phi0=phi_star%fx(1)
	end select

	if(phi0<0.001*enx.and.&
		(logrc<pxmin.or.(func(pxmin,par)>0.and.logrc>pxmin)).and.spp%mbh_dmless>0&
		.and. enx<phi0+spp%mbh_dmless**2/2d0/(jm*jc_dm)**2)then
		!print*, "phi_star%xb(1),fx(1)=",phi_star%xb(1), phi_star%fx(1)
		!print*, "enx,jm=",enx,jm
		!print*, "use fast"
		e0=phi0-enx
		
		rp=-spp%mbh_dmless/(2d0*e0)*&
		(1-(1+(jm*jc_dm)**2*2*e0/spp%mbh_dmless**2)**0.5d0)
		!print*, "logrc, pxmin,func(pxmin,par)=",logrc, pxmin,func(pxmin,par)
		!print*, "ex,jm,rp,mbh_dmless=", enx,jm,rp,spp%mbh_dmless
		!print*, "expect KP rp=", spp%mbh_dmless/(2*enx)*(1-(1-jm**2)**0.5)
		!print*, "phi0,pxmin=",phi0,10**pxmin
		!logrmin=logrc+log10(jmin_value)*2-2
		!r1=root_finding_1(func,logrmin,logrc,1d-10,10,par,ier,slient)
		!print*, "root finding=",r1,10**r1
		!read(*,*)
		if(rp<=0)then
			print*, 'error! rp<=0,rp, mbh_dmless', rp, spp%mbh_dmless
			print*, "phi0,enx,pxmin=",phi0,enx,pxmin
			print*, "e0,jm,jc_dm=",e0,jm,jc_dm
			call spp%fphi_star%print("phi_star")
			stop
		end if
		if((1+(jm*jc_dm)**2*2*e0/spp%mbh_dmless**2)<0)then
			print*, 'error! (1+(jm*jc_dm)**2*2*e0/mbh_dmless**2)<=0,rp, mbh_dmless=', &
			(1+(jm*jc_dm)**2*2*e0/spp%mbh_dmless**2), rp, spp%mbh_dmless
			print*, "phi0,enx,pxmin=",phi0,enx,pxmin, logrc
			print*, "func(pxmin,par)=",func(pxmin,par)
			print*, "e0,jm,jc_dm=",e0,jm,jc_dm
			call spp%fphi_star%print("phi_star")
			stop
		end if
	
	else
			logrmin=logrc+log10(jmin_value)*2-2
			edge0=func(logrmin,par)
			if(edge0<1d-8.and.edge0.ge.-1d-12) then
				r1=logrmin
				goto 100
			end if
			edge1=func(logrc,par)
			if(edge1<1d-8.and.edge1.ge.-1d-12)then
				r1=logrc
				goto 100
			end if
			if(func(logrmin,par)*func(logrc,par).gt.0)then
				print*, "in get_rpra_dmless"
				print*, "no root between rmin and rc"
				call error_handle()
			end if
			!print*, "logrc, logrmax=",logrc, logrmax+0.1
			!print*, func(logrc,par),func(logrmax+0.1,par)

			!r1=rtbis_yacc(func,logrmin,logrc,1d-8,par,ier,slient)
			r1=root_finding_1(func,logrmin,logrc,1d-8,10,par,ier,slient)
			!print*, "enx,r1,jm=",enx,10**r1, jm, func(r1,par)
			if(ier.eq.1)then
				if(.not.slient)then
					select case(ctl%ebin_type)
					case(ebin_type_log)
						print*, "func(phi_star%xb(1),par)>0, 10**phi_star%fx(1)*2<enx, logrc>pxmin"
						print*, func(spp%fphi_star%xb(1),par)>0, 10**spp%fphi_star%fx(1)*2<enx, logrc>pxmin
						print*, "r1 root not found closest=", r1, func(r1,par)
					case(ebin_type_lin)
						print*, "func(phi_star%xb(1),par)>0, phi_star%fx(1)*2<enx, logrc>pxmin"
						print*, func(spp%fphi_star%xb(1),par)>0, spp%fphi_star%fx(1)*2<enx, logrc>pxmin
						print*, "r1 root not found closest=", r1, func(r1,par)
					end select
					call error_handle()
				else
					select case(ctl%ebin_type)
					case(ebin_type_log)
						print*, "func(phi_star%xb(1),par)>0, 10**phi_star%fx(1)*2<enx, logrc>pxmin"
						print*, func(spp%fphi_star%xb(1),par)>0, 10**spp%fphi_star%fx(1)*2<enx, logrc>pxmin
						print*, "r1 root not found closest,r1,logrmin, func(r1),ex,jm=", r1, logrmin, func(r1,par), enx,jm
					case(ebin_type_lin)
						print*, "func(phi_star%xb(1),par)>0, 10**phi_star%fx(1)*2<enx, logrc>pxmin"
						print*, func(spp%fphi_star%xb(1),par)>0, spp%fphi_star%fx(1)*2<enx, logrc>pxmin
						print*, "r1 root not found closest,r1,logrmin, func(r1),ex,jm=", r1, logrmin,func(r1,par), enx,jm
					end select
					call error_handle()
				end if
			end if

		!print*, "2"
100			rp=10**r1
			if(ieee_is_nan(rp))then
				print*, "error! rp is NaN", logrmin, logrc, r1, rp
				
				stop
			end if
		
	end if	
	!block
	!	print*, func(log10(rp),par)
	!	print*, logrc<pxmin,func(pxmin,par)>0,10**phi_star%fx(1)*3<enx
	!	call get_phi_star_full_range(phi_star,log10(rp),phi_star_tmp)
	!		vr_tmp=2*(10**phi_star_tmp+mbh_dmless/rp-enx)-(jm*jc_dm)**2/rp**2
	!		if(abs(vr_tmp)>0.01)then
	!			print*, "get_rpra_dmless"
	!			print*, "vr(rp), rp, ex, jum =",vr_tmp, rp, enx,jm, logrc,pxmin
	!			print*, "e0=",e0,r1
	!			read(*,*)
	!		end if
	!end block

	edge0=func(logrc,par)
	if(edge0<1d-8.and.edge0>=-1d-12) then
		r2=logrc
		goto 200
	end if
	edge1=func(logrmax,par)
	if(edge1<1d-8.and.edge1>=-1d-12)then
		r2=logrmax
		goto 200
	end if
	if(func(logrc,par)*func(logrmax+0.5,par).gt.0)then
		print*, "in get_rpra_dmless"
		print*, "no root between rc and rmax"
		call error_handle()
	end if

	!r2=rtbis_yacc(func,logrc,logrmax+0.5,1d-8,par,ier,slient)
	r2=root_finding_1(func,logrc,logrmax+0.5,1d-8,10,par,ier,slient)
	if(ier.eq.1)then
		
		print*, "r2 root not found closest=", r2, func(r2,par)
		call error_handle()
	end if

200	ra=10**r2
	if(rp>ra)then
		print*, "get_rpra_dmless:error! rp>ra",rp,ra
		print*, "assuming within r?",(logrc<pxmin.or.(func(pxmin,par)>0.and.logrc>pxmin))&
				.and.spp%mbh_dmless>0,  enx<phi0+spp%mbh_dmless**2/2d0/(jm*jc_dm)**2
		print*, "phi0=",phi0
		stop
	end if

	!if(jm>1-1d-4)then
	!	print*, "rp,ra=",rp,ra, jm
		
		
	!	read(*,*)
	!end if
	!print*, "rp,ra=",rp,ra, func(log10(rp),par),func(log10(ra),par), jc_dm, jm
	!read(*,*)
	!if(func(r1,par)>0.001.or.func(r2,par)>0.001)then
	!	print*, "????"
	!	print*, "rp,ra=",rp,ra, func(log10(rp),par),func(log10(ra),par), jc_dm, jm
	!	read(*,*)
	!end if
	!read(*,*)
contains
	real(8) function func(x,par)
		implicit none
		real(8) x, par(50)

		!call get_r_root(phi_star,rc,enx,1d0, jc_dm,func)
		!print*, "func=", func
		!stop
        call get_r_root(spp,x, enx, jph2, func)
	end function
	subroutine error_handle()
		implicit none
		type(s1d_type)::test
		print*, "logex,ex,jm=",log10(enx),enx,jm
		print*, "logmin, logrc=", logrmin, logrc
		print*, "fun(min), fun(rc)=", func(logrmin, par), func(logrc,par)
		print*, "rc, max=", logrc, logrmax
		print*, "fun(rc), fun(max)=", func(logrc, par), func(logrmax,par)
		!call dms%fr_phi%print("rmax")
		!call dms%rc%print("rc")
		call spp%fphi_star%print("fphi")
		call test%init(logrmin,logrmax,100,sts_type_grid)
		call test%set_range()
		do i=1, test%nbin
			test%fx(i)=func(test%xb(i), par)
			!print*, test%xb(i),par(1)
		end do
		call test%print("test")
		stop
	end subroutine
end subroutine
subroutine get_rp_from_qp(spp, enx, jph2, rp)
	use constant
	use md_star_pot
	use model_basic,only:sample_logrmin,ctl
	implicit none
	real(8) jph2, rp
	type(star_pot_para)::spp
	real(8) ct,dt,et,ex,enx
	real(8) rho0, r1,e0
	real(8) r_root(4), i_root(4), pa(5)
	real(8) delta, Ix, Jx,wx, gx, Rx
	real(8) term1, term2,tx,sgng, tx_2, term3, term4, fx, term6,term7
	real(8) term5,kx
	integer nsol,i
	!logical::debug=.false.

	rho0=spp%spt_rho_rmin

	r1=10**sample_logrmin
	ct=-(2*pi*rho0*r1**2+spp%phi_r1r2_s-enx)/(4*pi*rho0)
	dt=-3*spp%mbh_dmless/8d0/pi/rho0
	et=3*jph2/4d0/pi/rho0 
	Ix=et+3*ct*ct
	Jx=ct*et-ct*ct*ct-dt*dt
	delta=Ix*Ix*Ix-27*Jx*Jx
	!print*, "ct,dt,et=",ct,dt,et,delta
	!do i=1, nsol
	!	if(r_root(i)<r1)then
	!rp=r_root(i)
	!read(*,*)
	if(delta>0)then
		print*, "stop:delta>0", delta
		print*, "ct,dt,et=",ct,dt,et
		print*, "enx=",enx
		stop
	end if

	term4=Ix*Ix*Ix/27d0/Jx/Jx
	delta=et/(ct*ct); wx=dt*dt/(ct**3)
	if(abs(term4)>1d-5)then
		if(term4<0.99999)then
			term1=-Jx*(1d0-(1-term4)**0.5d0)
			term2=-Jx*(1d0+(1-term4)**0.5d0)
			tx=0.5d0*(term1**(1d0/3d0)+term2**(1d0/3d0))-ct
			tx_2=tx**0.5
		else
		!	tx=(-Jx)**(1/3d0)/2d0*(2-2d0/9d0*(1-term4)-20d0/243d0*(1-term4)**2)-ct
			!print*, "tx=",tx, 0.5d0*(-Jx)**(1d0/3d0)*(((1d0-(1-term4)**0.5d0))**(1d0/3d0)&
		!		+((1d0+(1-term4)**0.5d0))**(1d0/3d0))-ct
			term6=wx-delta
			term7=(2*wx+wx*wx-3*delta+2d0/3d0*delta*delta-2*delta*wx-delta**3/27d0)/(1+wx-delta)**2
			tx=ct*(term6/3d0-term6*term6/9d0+term6**3*5d0/81d0 &
				-(term7/9d0+10d0/243d0*term7*term7+154d0/6561d0*term7**3)*(1+wx-delta)**(1/3d0))
			!print*, "tx=",tx, 0.5d0*(-Jx)**(1d0/3d0)*(((1d0-(1-term4)**0.5d0))**(1d0/3d0)&
			!	+((1d0+(1-term4)**0.5d0))**(1d0/3d0))-ct, term6
			tx_2=tx**0.5
		end if
	else
		term1=-Jx*(term4/2d0+term4*term4/8d0+term4**3/16d0+term4**4*5d0/128d0)
		term2=-Jx*(1d0+(1-term4)**0.5d0)
		tx=0.5d0*(term1**(1d0/3d0)+term2**(1d0/3d0))-ct
		tx_2=tx**0.5
		!print*, "term1:2=",term1, -Jx*(1d0-(1-term4)**0.5d0)
	end if

	fx=tx/ct
	ex=delta/(2*fx+3)**2
	term5=ex
	if(abs(term5)<1d-5)then
		term3=(3+2*fx)/fx*(term5/2d0+term5*term5/8d0+term5**3/16d0)
		!print*, "term3=",term3, dt/tx**1.5+3*ct/tx+2!, &
		!(3+2*fx)/fx*(1-(1+ex/fx)**0.5), ex/fx
		!read(*,*)
		!ctl%debug=1
	else
		term3=(3+2*fx)/fx*(1-(1-term5)**0.5) !dt/tx**1.5+3*ct/tx+2
!		print*, "term3:3, delta/wx, gx/fx=",term3, abs(delta/wx), abs(gx/fx)
	end if
	if(abs(term3)>1d-5)then
		if(term3<1)then
			rp=tx_2*(1-(1-term3)**0.5)
		else
			print*, "error! term3>1",term3
			goto 100
		end if
	else
		rp=tx_2*(term3/2d0+term3*term3/8d0+term3**3/16d0)
		!print*, "get_rp_from_qp:rp2=",rp, term3
	!	read(*,*)
	end if 
	if(rp<=0.or.ctl%debug.eq.1.or.isnan(rp))then		
100		print*, "rp,enx,jph2=",rp, enx,jph2
		print*, "tx_2=",tx_2
		print*, "(dt+3*ct*tx_2+2*tx**1.5d0)=",(dt+3*ct*tx_2+2*tx**1.5d0)
		print*, "tx,term3,fx,ex=",tx, term3,fx,ex
		print*, "ct,dt,et=",ct,dt,et
		print*, "delta, wx=", delta, wx
		print*, "term6,term7=",term6,term7
		print*, "term1,term2,term4,term5=",term1,term2,term4,term5
		read(*,*)
		ctl%debug=0
		stop
		!debug=.false.
	end if
	return
	
	!	end if
	!end do
	!print*, "nreal_sol=",nsol
	!print*, "r_root=",r_root
	!print*, "i_root=",i_root
	!read(*,*)

end subroutine
subroutine get_r_root(spp, logr, enx, jph2, func)
    use com_sts_type
    use model_basic,only:ctl
	use md_coeff
	use md_star_pot
    implicit none
    type(s1d_type)::phi_star
	type(star_pot_para)::spp
    real(8) phi_star_tmp, logr, func, r, enx, jph2
	
	!if(x>phi%xb(phi%nbin))
    !print*, phi_star%xmin, logr, phi_star%xmax

    call get_phi_star_full_range(spp,logr, phi_star_tmp)
	r=10**logr
	select case(ctl%ebin_type)
	case(ebin_type_log)
	    func=2*(10**phi_star_tmp+spp%mbh_dmless/r-enx)-jph2/r**2
	case(ebin_type_lin)
		func=2*(phi_star_tmp+spp%mbh_dmless/r-enx)-jph2/r**2
	end select
    
	
        
    !print*, "x, phi, phi_tot,func=",enx,10**phi_star_tmp,10**phi_star_tmp+1d0/r, func, &
	!	2*(10**phi_star_tmp+1d0/r-enx)*r**2

end subroutine
  

real(8) function p_EJ_dmless_fast(spp,enx,j,jc_dmless,rp,ra)
	use model_basic
	use md_star_pot
	implicit none
	real(8) enx,j
	type(star_pot_para)::spp
	!type(s1d_type)::phi_star,fma
	!type(s1d_type)::frho
	real(8) ra,rp,j_dmless2,jc_dmless
	integer i,idid, nbin
	integer idxmin, idxmax
	real(8) phi_star_tmp,yout, up0, up1, umin, umax, logrmin, logrmax, logrminmid,logrmaxmid
	real(8) logr0,logr1,res_tmp, yp1, yp0,xp1, xp0,logrmid
	integer flag_end,idx1,idx0
	real(8) line_c, line_slope, px, qx,ax, fu1, fu0, term1, term2, term0
	
	if(rp.eq.ra)then
		call get_phi_star_full_range(spp,log10(rp),phi_star_tmp)
		p_EJ_dmless_fast=2*pi*ra/(2*(10**phi_star_tmp+spp%mbh_dmless/rp-enx))**0.5
		return
	end if
	nbin=common_aux%nbin

	call get_aux_function_for_period_pi2(common_aux,spp,enx,j,jc_dmless,rp,ra)
	res_tmp=0
	do i=1, nbin-1
		yp1=common_aux%fx(i+1)
		yp0=common_aux%fx(i)
		xp1=common_aux%xb(i+1)
		xp0=common_aux%xb(i)
		res_tmp=res_tmp+(yp1+yp0)*(xp1-xp0)
	end do
	p_EJ_dmless_fast=res_tmp

end function


real(8) function p_EJ_dmless(spp,enx,j,jc_dmless,rp,ra)
	use model_basic
	use md_star_pot
	implicit none
	real(8) enx,j
	type(star_pot_para)::spp
	!type(s1d_type)::phi_star,fma
	!type(s1d_type)::frho
	real(8) ra,rp,j_dmless2,jc_dmless
	integer i,idid
	real(8)   phi_star_tmp,yout
	
	if(rp.eq.ra)then
		call get_phi_star_full_range(spp,log10(rp),phi_star_tmp)
		p_EJ_dmless=2*pi*ra/(2*(10**phi_star_tmp+spp%mbh_dmless/rp-enx))**0.5
		return
	end if
	call get_aux_function_for_period_pi2(common_aux,spp,enx,j,jc_dmless,rp,ra)
	!print*,"enx, j, jc, rp, ra=", enx, j,jc_dmless, rp, ra
	!if(ctl%debug.ge.2)then
	!	call aux%print("aux")
	!end if
	!read(*,*)

	yout=0

	call my_integral_acc(0d0,pi/2d0,yout,pd_int_acc_a,pd_int_acc_r, fcn,idid)
	p_EJ_dmless=yout*2
	if(idid<0)then
		print*, "error in p_EJ_dmless"
		call common_aux%print("common_aux")
		stop
	end if
	if(p_EJ_dmless<0)then
		print*, "error! p_EJ_dmless<0??", p_EJ_dmless
		call common_aux%print("common_aux")
		print*, "enx, jm, jc_dmless, rp, ra=",enx, j, jc_dmless, rp, ra
		stop
	end if
	!print*, "p_EJ_dm=",p_EJ_dmless
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
			call spp%fphi_star%print("phi_star")
			call spp%frho_star%print("frho")
			call spp%fma_star%print("fma")
			stop
		end if
	end subroutine
end function
     

real(8) function jc_dmless(rc_EJ,spp)
	use com_sts_type
	use md_star_pot
	implicit none	
	real(8)  rc_EJ
	real(8) mfunc, r0_cl
	type(star_pot_para)::spp
	!type(s1d_type)::fma
	!rc_EJ=r_c(phi,fma,E)
	!call fbeta%get_value_l(log10(rc_EJ),mfunc)
	call get_beta_full_range(spp,log10(rc_EJ),mfunc)
	!call get_value_at_x_fc(fma, log10(rc_EJ), yout, 1)
	!mfunc=yout(1)
	!dpdr=-mbh/rc_EJ**2-mfunc/rc_EJ**2
	jc_dmless=((spp%mbh_dmless+mfunc)*rc_EJ)**0.5d0
end function
real(8) function jc_EJ(jc_dm,mbhin,r0_cl)
	use com_sts_type
	implicit none	
	real(8) mbhin, jc_dm
	real(8)  r0_cl
	
	jc_EJ=jc_dm*(r0_cl*mbhin)**0.5d0
end function
subroutine get_xrc(xrc,spp)
	use com_sts_type
	use md_star_pot
	use md_coeff
	use model_basic,only:ctl
	implicit none
	type(s1d_type)::xrc, rcx
	type(star_pot_para)::spp
	real(8) phi_tmp,beta_tmp,r
	integer i

	call xrc%init(ctl%log10rmin_factor,ctl%log10rmax_factor,ctl%dstr_bins_r,sts_type_dstr)
	!call xrc%set_range()
	call rcx%init(ctl%log10rmin_factor,ctl%log10rmax_factor,ctl%dstr_bins_r,sts_type_dstr)
	call rcx%set_range()

	do i=1, xrc%nbin
		call get_phi_star_full_range(spp,rcx%xb(i),phi_tmp)
		call get_beta_full_range(spp,rcx%xb(i),beta_tmp)
		r=10**rcx%xb(i)
		rcx%fx(i)=log10((spp%mbh_dmless+2*10**phi_tmp*r-beta_tmp)/2d0/r)
	end do
	xrc%xb(1:xrc%nbin)=rcx%fx(rcx%nbin:1:-1)
	xrc%xmax=rcx%fx(1)
	xrc%xmin=rcx%fx(rcx%nbin)
	xrc%fx(1:xrc%nbin)=rcx%xb(rcx%nbin:1:-1)

end subroutine
subroutine get_x_given_jc(rcmin,rcmax,logx,jc,spp,dm)
	use md_star_pot
	use md_dms
	use model_basic,only:ctl
	implicit none
	real(8) logx, jc
	type(star_pot_para)::spp
	type(diffuse_mspec)::dm
	real(8) rtbis_yacc, logrc,par(50), phi_tmp, beta_tmp
	integer ier,niter
	real(8) jcmincrit, rcmin,rcmax
	!rmin=10**dm%logrmin
	!jcmincrit=(rmin*(rmin**3*4*pi/3d0*spp%spt_rho_rmin+spp%mbh_dmless))**0.5d0
	!print*, "rmin,jc,jcmincrit=",rmin,jc,jcmincrit
	!if(jc>jcmincrit.or.spp_new%mbh.eq.0)then
	!print*, "jc=",jc
		logrc=rtbis_yacc(func,rcmin,rcmax,&
		1d-13, par,niter,1000,ier,.true.)
	!else
	!	call solve_rc_given_jc_ana(logrc,jc,spp)
	!	logrc=
	!end if
	
	call get_phi_star_full_range(spp,logrc,phi_tmp)
	call get_beta_full_range(spp,logrc,beta_tmp)
	logx=log10((spp%mbh_dmless+2*10**phi_tmp*10**logrc-beta_tmp)&
		/2d0/10**logrc)
	!print*, "logrc,jc, logx=",logrc,jc, logx
	!read(*,*)
contains
	real(8) function func(x,par)
		implicit none
		real(8) x, par(50)
		call get_beta_full_range(spp,x,beta_tmp)
		func=1-(spp%mbh_dmless+beta_tmp)*10**x/jc**2
		!print*, "x, func=",x,func
	end function
end subroutine
subroutine solve_rc_given_jc_ana(logrc,jc,spp)
	use md_star_pot
	use constant
	implicit none
	type(star_pot_para)::spp
	real(8) logrc, jc
	real(8) d,e, delta
	d=3*spp%mbh_dmless/4d0/pi/spp%spt_rho_rmin; 
	e=-3*jc**2/4d0/pi/spp%spt_rho_rmin
	delta=e**3-27*d**4
	!print*, "deta=",delta
	!read(*,*)
end subroutine
subroutine get_jc_minmax(xmin,xmax,rcmin,rcmax,jcmin,jcmax,spp)
	use md_star_pot
	implicit none
	real(8) jcmin,jcmax,xmin,xmax,rcmin,rcmax
	real(8) r_c_iter,jc_dmless
	integer ier
	type(star_pot_para)::spp

	rcmax=r_c_iter(spp,10**xmin,ier)
    
    jcmax=jc_dmless(rcmax,spp)
    !print*, "xmin,rcmax=",xmin,rcmax,jcmax

    rcmin=r_c_iter(spp,10**xmax,ier)
    jcmin=jc_dmless(rcmin,spp)
	rcmin=log10(rcmin); rcmax=log10(rcmax)
    !print*, "xmax,rcmin,jcmin=",xmax,rcmin,jcmin
	
end subroutine
subroutine get_rc(rc,  spp)
	use com_sts_type
	use md_coeff
	use model_basic,only:ctl
	use md_star_pot
	use mpi_comu,only:rid
	implicit none
	type(s1d_type)::rc!,xrc
	type(star_pot_para)::spp
	real(8) r_c,r_c_iter
	integer i,ier
	!real(8) t1, t2
	!call cpu_time(t1)
	!call alpha%init_intrp()
	!call beta%init_intrp()
	
	do i=1, rc%nbin
		select case(ctl%ebin_type)  
		case(ebin_type_log)
			!rc%fx(i)=r_c(spp,10**rc%xb(i),ier)
			!if(rid.eq.0)then
			!	print*, "i,x=", i,10**rc%xb(i)
			!end if
			rc%fx(i)=r_c_iter(spp,10**rc%xb(i),ier)
			
		case(ebin_type_lin)
			rc%fx(i)=r_c(spp,rc%xb(i),ier)
		end select
		!if(ier.ne.0)then
		!	print*, "in get_rc,i=", i
			!stop
		!end if
	end do
	!stop
	!call cpu_time(t2)
	!print*, "t2-t1", t2-t1
end subroutine

real(8) function r_c_iter_re(spp,ex,ier) ! radius of circular orbit
	use com_sts_type
	use md_coeff
	use model_basic,only:ctl
	use md_star_pot
	use constant
	implicit none	
	type(star_pot_para)::spp
	!type(s1d_type)::fphi, fma
	real(8) ex,par(50),rmin,rmax, logrc, rc_tmp
	real(8) edge_value_min,edge_value_max,r1, r2
	integer max_iteration,iter,ier
	real(8) require_acc,r_c_given_rmin_rmax, rc2

	!if((emax_factor-ex)<1d-10.and.mbh_dmless.eq.0)then
	!	ier=0
	!	r_c=(6*(emax_factor-ex)/(4*pi*spt_rho_rmin))**0.5d0
	!	return
	!end if
	require_acc=1d-12

	max_iteration=8000

	rmin=spp%fphi_star%xmin
	rmax=spp%fphi_star%xmax

	edge_value_max=func_rc(rmax)
	edge_value_min=func_rc(rmin)
	if(abs(edge_value_max)<1d-15)then
		r_c_iter_re=10**rmax
		return
	end if
	if(abs(edge_value_min)<1d-16)then
		r_c_iter_re=10**rmin
		! ex is very close to emax (if mbh=0)
		!print*, "abs(edge_value_min)<1d-16, r_c=", r_c, edge_value_min, ex, rmin
		return
	end if
	
	r1=rmin
	iter=0
100 r2=func(r1)
	!print*, "r1,r2=",r1, r2,r2-r1
	iter=iter+1
	if(iter<max_iteration.and.abs(r2-r1)>require_acc)then
		
		r1=r2
		goto 100
	else
		!if(iter.eq.max_iteration.or. abs(r2-r1)>require_acc)then
			!print*, "error! iter, r1,r2=",iter, r1,r2
			!read(*,*)
		!end if
		IF(abs(R2-R1)>1e-5)THEN
			print*, "iter=",iter
			print*, "rc_iter: error!: r1,r2,dr,ex=",r1,r2,r2-r1,ex
			print*, "r2,func(r2)=",r2,func(r2), func(func(r2))
			stop
		end if
	end if
	r_c_iter_re=10**r1
	!rc2=r_c_given_rmin_rmax(spp,ex,spp%fphi_star%xmin,spp%fphi_star%xmax,ier)
	!if(abs(rc2-r_c_iter)>1d-6)then
	!	print*, "rc2,rc=", rc2, r_c_iter, rc2-r_c_iter
	!	read(*,*)
	!end if
	
	!read(*,*)
	!print*, "ex, rc=", ex, r_c_iter
	!stop

contains
	real(8) function func(x)
		implicit none
		real(8) x,phifun,r,betafun, rcv
		!call get_alpha_full_range(fphi,x,alphafun)
		call get_phi_star_full_range(spp,x,phifun)
		call get_beta_full_range(spp,x,betafun)
		!r=10**x
		rcv=(spp%mbh_dmless+2*10**phifun*10**x-betafun)/(2d0*ex)
		!print*, rcv
		!read(*,*)
		func=log10(rcv)

	end function
	real(8) function func_rc(x)
		implicit none
		real(8) x,phifun,r,betafun
		!call get_alpha_full_range(fphi,x,alphafun)
		call get_phi_star_full_range(spp,x,phifun)
		call get_beta_full_range(spp,x,betafun)
		r=10**x
		select case(ctl%ebin_type)
		case(ebin_type_log)
			func_rc=2*ex-(spp%mbh_dmless/r+2*10**phifun-betafun/r)
		case(ebin_type_lin)
			func_rc=2*ex-(spp%mbh_dmless/r+2*phifun-betafun/r)
		case default
			print*, "error in rc"
			stop
		end select
		!print*, "ex, x, alphafun,betafun=", ex, x, alphafun,betafun
	end function
end function

real(8) function r_c_iter(spp,ex,ier) ! radius of circular orbit
	use com_sts_type
	use md_coeff
	use model_basic,only:ctl,dms
	use md_star_pot
	use constant
	implicit none	
	type(star_pot_para)::spp
	!type(s1d_type)::xrc
	!type(s1d_type)::fphi, fma
	real(8) ex,emid,r_c_iter_re,r_c_given_rmin_rmax 
	real(8) r_c_given_r0
	integer ier
	!real(8) t1,t2
	integer rc_method 
	rc_method=4

	select case(rc_method)
	case(2)
		r_c_iter=r_c_given_rmin_rmax(spp, ex, spp%fphi_star%xmin,spp%fphi_star%xmax, ier)
	case(3)
		r_c_iter=r_c_iter_re(spp,ex,ier)
	case(4)
		r_c_iter=r_c_given_r0(spp,ex,dms%frc_x,ier)
		!print*, "r_c_iter=", r_c_iter
		!read(*,*)
	end select
		!call cpu_time(t2)
		!print*, "r_c_iter=",r_c_iter,t2-t1
		!read(*,*)
	!end if
end function

 
real(8) function r_c_given_rmin_rmax(spp,ex,rmin_in,rmax_in, ier)
	use com_sts_type
	use md_coeff
	use model_basic,only:ctl
	use md_star_pot
	use constant
	implicit none	
	type(star_pot_para)::spp
	!type(s1d_type)::fphi, fma
	real(8) ex,par(50),rmin,rmax, logrc,rmax_in,rmin_in
	integer ier
	real(8) width,edge_value_min,edge_value_max,root_finding_1
	integer max_iteration, debug
	debug=0
	max_iteration=10
	!print*, "ex=",ex, mbh, r_c	
	!rmin=0.02*mbh/(2.*ex); rmax=50.*mbh/(2*ex)
	width=2
100	rmin=rmin_in-width; rmax=rmax_in+width
	!print*, "rmin,rmax=",rmin,rmax
	!print*, "ex=",ex
	edge_value_max=func(rmax,par)
	edge_value_min=func(rmin,par)
	if(abs(edge_value_max)<1d-15)then
		r_c_given_rmin_rmax=10**rmax
		return
	end if
	if(abs(edge_value_min)<1d-16)then
		r_c_given_rmin_rmax=10**rmin
		! ex is very close to emax (if mbh=0)
		!print*, "abs(edge_value_min)<1d-16, r_c=", r_c, edge_value_min, ex, rmin
		return
	end if
	!if(ex>1d0)then
	!	print*, "rmin,rmax=",rmin,rmax, edge_value_max,edge_value_min
	!	debug=1
	!end if
	if(edge_value_max*edge_value_min>0)then
		if(max_iteration>0)then
			width=width*2
			max_iteration=max_iteration-1
			goto 100
		else
			print*, "rc error!"
			print*, "ex,ex>emax_factor,emax_factor=",&
				ex,ex>emax_factor,emax_factor
			print*, "rmin,rmax=",rmin,rmax
			print*, func(rmin,par), func(rmax,par)
			call spp%fphi_star%print("fphi")
			call spp%fma_star%print("fma")
			stop
		end if
	end if
	
	logrc=root_finding_1(func,rmin,rmax,1d-12,20,par,ier,.true.)
	!logrc= rtbis_yacc(func,rmin,rmax,1d-15,par,ier,.true.)
	select case(ier)
	case(1)
		print*, "r_c: ex=",ex
		call spp%fphi_star%print("fphi")
		call spp%fma_star%print("fma")
		stop
	case(2)
		!print*, "logrc,dt=",logrc,func(logrc,par)
	end select
	r_c_given_rmin_rmax=10**logrc

	!if(ex>10**sample_logemax)then
	!	print*, "r_c=",r_c, (1.5d0*10**(sample_logrmin*2)-0.75*(ex-phi_r1r2_s)/pi/spt_rho_rmin)
	!	stop
	!end if
	contains
	real(8) function func(x,par)
		implicit none
		real(8) x,phifun,r,par(50),betafun
		!call get_alpha_full_range(fphi,x,alphafun)
		call get_phi_star_full_range(spp,x,phifun)
		call get_beta_full_range(spp,x,betafun)
		r=10**x
		select case(ctl%ebin_type)
		case(ebin_type_log)
			func=2*ex-(spp%mbh_dmless/r+2*10**phifun-betafun/r)
		case(ebin_type_lin)
			func=2*ex-(spp%mbh_dmless/r+2*phifun-betafun/r)
		case default
			print*, "error in rc"
			stop
		end select
		if(debug.eq.1)then
			print*, "ex, r, alphafun,betafun=", ex, r, phifun,betafun
		end if
	end function
end function
subroutine get_rc_r_le_r1(spp, ex,rc)
	use model_basic
	use md_star_pot
	implicit none
	type(star_pot_para)::spp
	real(8) ex,rc
	real(8) px, qx, r1
	real(8) delta, term,rho0

	r1=10**sample_logrmin
	rho0=spp%spt_rho_rmin
	if(rho0.ne.0d0)then
		px=-(2*pi*rho0*r1*r1+spp%phi_r1r2_s-ex)/(4*pi*rho0)
		qx=-3*spp%mbh_dmless/(16d0*pi*rho0)
		!delta=qx**2+px**3
		term=px**3/qx**2
		if(term>0.0001)then
			rc=(-qx)**(1/3d0)*((1+(1+term)**0.5)**(1/3d0)+(1-(1+term)**0.5)**(1/3d0))
		else
			rc=(-qx)**(1/3d0)*((1+(1+term)**0.5)**(1/3d0)&
			+(-term/2d0+term**2/8d0-term**3/16d0)**(1/3d0))
		end if
	else
		rc=-spp%mbh_dmless/2d0/(spp%phi_r1r2_s-ex)
	end if
	!rc=(-qx)+delta**0.5)**(1/3d0)+(-qx-delta**0.5)**(1/3d0)
	
	!print*, "rc,mbh/2ex=",rc,spp%mbh_dmless/ex/2d0
	!read(*,*)
	!print*, "r1,cs=",r1, spp%phi_r1r2_s
	!print*, "delta, px,qx=",delta, px, qx
	!print*, "mbh,logex,logrmax=",spp%mbh_dmless,logex,logrmax, 10**logrmax
	!read(*,*)
end subroutine
real(8) function r_c_given_r0(spp,ex,xrc, ier)
	use com_sts_type
	use md_coeff
	use model_basic,only:ctl,sample_logemax,sample_logrmin
	use md_star_pot
	use constant
	implicit none	
	type(star_pot_para)::spp
	type(s1d_type)::xrc
	!type(s1d_type)::fphi, fma
	real(8) ex,par(50), logrc,r0
	integer ier,niter
	real(8) width,xl,xh,require_acc, rtbis_yacc
	integer nmax_iter, debug
	!if(abs(ex-3104.8909)<0.001)then
	!	debug=1
	!	print*, "debug=1"
	!else
		debug=0
	!end if

	nmax_iter=12
	!print*, "ex=",ex, mbh, r_c	
	!rmin=0.02*mbh/(2.*ex); rmax=50.*mbh/(2*ex)

	!call xrc%print("xrc")

	call xrc%get_value_l(log10(ex),r0)
	!print*, "ex,r0=",ex,r0,log10(ex),10**r0
	require_acc=1d-13
	
	width=0.1d0

	!ex=1d8
	if(spp%mbh_dmless.ne.0)then
		require_acc=1d-11
		! if(ex>spp%phi)
	else
100		if(ex>10**sample_logemax)then
			if(ex-2*pi*spp%spt_rho_rmin*10**(sample_logrmin*2)&
				-spp%phi_r1r2_s>0)then
				!if(debug.eq.1)then
				!	print*, "1:", ex, sample_logemax, spp%spt_rho_rmin
				!end if
				
				call get_rc_r_le_r1(spp,ex,r_c_given_r0)
				!print*, "rc_given_r0=",r_c_given_r0
				return
			end if
		end if
	end if
	!if(mbh_dmless.ne.0)then
	xl=r0-width
	xh=r0+width
	logrc=rtbis_yacc(func,xl,xh,require_acc,par,niter,2000,ier,.true.)
		
	select case(ier)
	case(1)
		!if(debug.eq.1)then
		!	print*, "1:nmax_iter, r0, xl, xh=",nmax_iter, r0, xl, xh, require_acc
		!end if
		!call error_handle()
		if(nmax_iter>=0)then
			width=width/5d0
			require_acc=require_acc*2
			nmax_iter=nmax_iter-1
			r0=logrc
			goto 100
		else
			print*, "r_c_r0 no enough iteration: ex,width,xl,xh=",ex, width,xl,xh
			print*, "r0,nmax_iter=",r0, nmax_iter
			call spp%fphi_star%print("fphi")
			call spp%fma_star%print("fma")
			print*, "fun(xl),fun(xh),fun(rc)=",func(xl,par),func(xh,par),func(logrc,par)
			stop
		end if
		
	case(2)
		if(nmax_iter>0)then
			width=width*5d0	
			nmax_iter=nmax_iter-1
			if(debug.eq.1) then
				print*, "2:nmax_iter=",nmax_iter,  require_acc
			end if
			goto 100
		else
			print*, "r_c_given_r0:roots not bracket, ex, xl, xh:", ex, xl, xh
			print*, "r0=",r0
			call xrc%print("xrc")
			ier=-99
			return
		end if
		!print*, "logrc,dt=",logrc,func(logrc,par)
	end select
	r_c_given_r0=10**logrc
contains
	real(8) function func(x,par)
		implicit none
		real(8) x,phifun,r,par(50),betafun
		!call get_alpha_full_range(fphi,x,alphafun)
		call get_phi_star_full_range(spp,x,phifun)
		call get_beta_full_range(spp,x,betafun)
		r=10**x
		!select case(ctl%ebin_type)
		!case(ebin_type_log)
			func=2*ex-(spp%mbh_dmless/r+2*10**phifun-betafun/r)
		!case(ebin_type_lin)
		!	func=2*ex-(spp%mbh_dmless/r+2*phifun-betafun/r)
		!case default
		!	print*, "error in rc"
		!	stop
		!end select
		if(debug.eq.1)then
			print*, "ex, r, alphafun,betafun=", ex, r, phifun,betafun,func
		end if
	end function
end function

real(8) function r_c(spp,ex,ier) ! radius of circular orbit
	use md_star_pot
	implicit none
	type(star_pot_para)::spp
	real(8) ex,r_c_given_rmin_rmax
	integer ier
	r_c=r_c_given_rmin_rmax(spp, ex, spp%fphi_star%xmin,spp%fphi_star%xmax, ier)
end function
subroutine get_nx_tmp(nx_tmp, nx_org)
	use com_sts_type
	implicit none
	type(s1d_type)::nx_tmp, nx_org
	integer i, nn
	nn=0
	do i=1, nx_org%nbin
		if(nx_org%fx(i).ne.-100d0)then
			nn=nn+1
		end if
	end do
	call nx_tmp%init(nx_org%xmin,nx_org%xmax,nn,nx_org%bin_type)
	nn=0
	do i=1, nx_org%nbin
		if(nx_org%fx(i).ne.-100d0)then
			nn=nn+1
			nx_tmp%xb(nn)=nx_org%xb(i)
			nx_tmp%fx(nn)=nx_org%fx(i)
		end if
	end do
	nx_tmp%xmin=nx_tmp%xb(1)-nx_tmp%xstep/2d0
	nx_tmp%xmax=nx_tmp%xb(nn)+nx_tmp%xstep/2d0
end subroutine

subroutine get_bin_number()
	use com_main_gw
	implicit none
	integer i, n, idid,clone_factor, midx
	real(8) Ntot(ctl%m_bins), emin, emax,weight_tot
	real(8),external::get_clone_deep
	type(s1d_type)::nx_tmp

	do i=1, dms%n
		call dms%mb(i)%all%fna%get_value_s(sample_logrmax,ctl%bin_total_number(i))
		!call dms%mb(i)%all%fna%print("fna")
	end do
	!do i=1, dms%n
	!	ctl%asymptot_ini(1,i)=ctl%bin_total_number(i)/spp_new%N_r_within_max
		!print*, "i,asymptot_ini=",i,ctl%asymptot_ini(1,i),ctl%bin_total_number(i),spp_new%N_r_within_max
	!end do
	! print*, ctl%bin_total_number(1:dms%n)
	! read(*,*)
	!weight_tot=0
	do i=1, ctl%m_bins
		n=ctl%ini_nx_log(i)%nbin
		
		clone_factor=ctl%clone_factor(i)
		!weight_tot=weight_tot+ctl%asymptot_ini(1,i)*ctl%bin_mass(i)
	!call ctl%ini_nx_tot%init_intrp()
	!print*, "clone_factor=",clone_factor,clone_e0_factor
		Ntot(i)=0
		
		call get_nx_tmp(nx_tmp, ctl%ini_nx_log(i))
		emin=nx_tmp%xmin
		emax=nx_tmp%xmax

		! call ctl%ini_nx_log(i)%print("nx i")
		! call nx_tmp%print("nx_tmp")
		! print*, "i=",i
		call my_integral_none(emin, emax, Ntot(i), FCN, idid)
		! print*, "Ntot=", Ntot(i)
		!  read(*,*)
	end do
	! print*, "111111"
	! stop
	do i=1, ctl%m_bins
		ctl%bin_mass_particle_number(i)=Ntot(i)*ctl%n0*r0_cl**3/ctl%ini_weight_n(i)!/weight_tot!*ctl%asymptot_ini(1,i)
		!print*, "i, bin_num, weightn, ntot=", i, ctl%bin_mass_particle_number(i),ctl%ini_weight_n(i),Ntot(i)
		!read(*,*)
		ctl%ini_nper_bin(i)%fx=ctl%ini_nper_bin(i)%fx*ctl%n0*r0_cl**3/ctl%ini_weight_n(i)*ctl%ntasks!**0.5
		!if(rid.eq.0)then
			!print*, "ntasks=",ctl%ntasks
			!call ctl%ini_nper_bin(i)%print("nper bin")
			!read(*,*)
		!end if
		
		! print*, "rid,i, n(i), nr(i),mtot,iniw=",rid,i, &
		! 	sum(ctl%ini_nper_bin(i)%fx), ctl%bin_mass_particle_number(i),ctl%n0*r0_cl**3, ctl%ini_weight_n(i)
	
	end do
	! print*, "22222"
	!print*, "ctl%n0, r0_cl=", ctl%n0, r0_cl,ctl%n0*r0_cl**3
	!print*, "ctl%clone_scheme=",ctl%clone_scheme
	!print*, "ntot=", ntot, ctl%bin_mass_particle_number(1)
	!stop
	
contains 
	subroutine FCN(N,X,Y,F,IPAR,RPAR)
		implicit none
		integer n, ipar(100)
		real(8) x, y(n), f(n), rpar(100),ysp
		integer nlvl

		call nx_tmp%get_value_s(X,ysp)
		if(ctl%clone_scheme.ge.1)then
			nlvl=int(get_clone_deep(10**X,log_clone_bd_sep,  clone_e0_factor))
			F(1)=10**ysp*dble(clone_factor)**nlvl
			! print*, "x,nlvl,clone_e0_factor=",x,nlvl,clone_e0_factor
		else
			F(1)=10**ysp
		end if
		! print*, "F=",x, F(1)
	end subroutine
end subroutine

subroutine check_aux(aux,spp,ex, jum,jc, rp,ra)
	use model_basic
	use, intrinsic::ieee_arithmetic
	use md_star_pot
	implicit none
	type(s1d_type)::aux
	!type(s1d_type)::phi_star, fma
	type(star_pot_para)::spp
	real(8) rp, ra, beta_rp, beta_ra, jum, jc, r,ex,sintheta,costheta
	real(8) phi_star_tmp, vr_tmp
	integer i
	
	do i=1, aux%nbin-1
		if(aux%fx(i)>aux%fx(i+1))then
			print*, "error detected"
			call aux%print("aux")

			sintheta=sin(aux%xb(i))
			costheta=cos(aux%xb(i))
			r=rp+(ra-rp)*sintheta*sintheta
		
			call get_phi_star_full_range(spp,log10(r),phi_star_tmp)
			vr_tmp=2*(10**phi_star_tmp+spp_new%mbh_dmless/r-ex)-(jum*jc)**2/r**2

			print*, "aux%xb(i),r,2*(ra-rp)*sintheta*costheta, vr_tmp, (vr_tmp)**(-0.5d0)"
			print*, aux%xb(i),r, 2*(ra-rp)*sintheta*costheta, vr_tmp, (vr_tmp)**(-0.5d0)
			
			call get_phi_star_full_range(spp,log10(ra),phi_star_tmp)
			vr_tmp=2*(10**phi_star_tmp+spp_new%mbh_dmless/ra-ex)-(jum*jc)**2/ra**2
			print*, "vr(ra)=",vr_tmp

			call get_phi_star_full_range(spp,log10(rp),phi_star_tmp)
			vr_tmp=2*(10**phi_star_tmp+spp_new%mbh_dmless/rp-ex)-(jum*jc)**2/rp**2
			print*, "vr(rp)=",vr_tmp
			stop
		end if
	end do
end subroutine
subroutine get_aux_function_for_period(aux,phi_star,fma,ex, jum,jc, rp,ra)
	use model_basic
	use md_star_pot
	use, intrinsic::ieee_arithmetic
	implicit none
	type(s1d_type)::aux
	type(s1d_type)::phi_star, fma
	real(8) rp, ra, beta_rp, beta_ra, jum, jc, r,ex,sintheta2
	real(8) phi_star_tmp, vr_tmp
	integer i
	integer,parameter:: nbins=200
	call aux%init(0d0,1d0,nbins,sts_type_grid)
	call aux%set_range()
	!call beta%print("beta")
	call get_beta_full_range(spp_new,log10(rp), beta_rp)
	call get_beta_full_range(spp_new,log10(ra), beta_ra)
	!print*, "rp, ra, beta_rp, beta_ra=", rp, ra, beta_rp, beta_ra
	!read(*,*)
	aux%fx(1)=2**0.5*(ra-rp)**0.5*(-(spp_new%mbh_dmless+beta_rp)/rp**2+jum**2*jc**2/rp**3)**(-0.5)
	aux%fx(nbins)=2**0.5*(ra-rp)**0.5*((spp_new%mbh_dmless+beta_ra)/ra**2-jum**2*jc**2/ra**3)**(-0.5)
	!call beta%print("beta")
	!print*, "beta_rp, rp=", beta_rp, rp,log10(rp)
	!print*, "f1,fn=",aux%fx(1), aux%fx(nbins)
	if(ctl%debug.ge.2)then
		call get_phi_star_full_range(spp_new,log10(rp),phi_star_tmp)
		print*,  "aux at rp:x,f=", 2*(10**phi_star_tmp+spp_new%mbh_dmless/rp-ex)-(jum*jc)**2/rp**2
	end if
	do i=2, nbins-1
		
		sintheta2=aux%xb(i)
		!theta=0d0
		r=rp+(ra-rp)*sintheta2
		
		call get_phi_star_full_range(spp_new,log10(r),phi_star_tmp)
		vr_tmp=2*(10**phi_star_tmp+spp_new%mbh_dmless/r-ex)-(jum*jc)**2/r**2
		if(vr_tmp<=0d0)then
			print*, "error! what happened??"
			print*, "r, vr_tmp=",r, vr_tmp
			print*, "sintheta2,r, phi_star_tmp=",sintheta2,r, phi_star_tmp
			print*, "ex,jum=", ex,jum, 2*(10**phi_star_tmp+spp_new%mbh_dmless/r-ex)-(jum*jc)**2/r**2
			call fma%print("fma")
			call phi_star%print("phi_star")
			print*, "ra,rp=",ra,rp
			print*, "jc=", jc
			stop
		end if
		aux%fx(i)=2*(ra-rp)*(sintheta2*(1-sintheta2))**0.5d0*&
			(vr_tmp)**(-0.5d0)
		!print*, "x,f=",theta,r, phi_star_tmp,ex,jum, 2*(10**phi_star_tmp+1d0/r-ex)-(jum*jc)**2/r**2
		!read(*,*)
		
		if(ieee_is_nan(aux%fx(i)).or.(.not.ieee_is_finite(aux%fx(i))))then
			print*, "error!"
			print*, "sintheta2,r, phi_star_tmp=",sintheta2,r, phi_star_tmp
			print*, "ex,jum=", ex,jum, 2*(10**phi_star_tmp+spp_new%mbh_dmless/r-ex)-(jum*jc)**2/r**2
			call fma%print("fma")
			call phi_star%print("phi_star")
			print*, "ra,rp=",ra,rp
			print*, "jc=", jc
			stop
		end if
	end do
end subroutine

subroutine get_aux_function_for_period_pi2(aux,spp,ex, jum,jc, rp,ra)
	use model_basic
	use, intrinsic::ieee_arithmetic
	use md_coeff
	use md_star_pot
	implicit none
	type(s1d_type)::aux
	!type(s1d_type)::phi_star, fma, frho
	type(star_pot_para)::spp
	real(8) rp, ra, jum, jc, r,ex,sintheta,costheta
	real(8) phi_star_tmp, vr_tmp
	real(8) logrc,rmax, rpx,rax,r_c_iter
	integer ier
	integer i,j
	integer nbins
	nbins=aux%nbin

	!call aux%init(0d0,pi/2d0,nbins,sts_type_grid)

	!call aux%set_range()

	! sampling more near pericenter

	!do i=1, aux%nbin
	!	aux%xb(i)=pi/2d0*((i-1)/dble(nbins-1))**2
	!end do

	!call beta%print("beta")
    !print*, "ex,jum=",ex,jum
	!call get_aux_2bd(phi_star, frho, fma,ra,rp,jc,ex,jum,aux%fx(1),aux%fx(aux%nbin))
	call get_aux_ra(spp,rp,ra, jc,ex,jum,aux%fx(aux%nbin))
	if(aux%fx(aux%nbin)<-1d90)then
		return
	end if
	call get_aux_rp(spp,rp,ra, jc,ex,jum,aux%fx(1))

	if(ctl%debug.ge.2)then
		call get_phi_star_full_range(spp,log10(rp),phi_star_tmp)
		select case(ctl%ebin_type)
		case(ebin_type_log)
			phi_star_tmp=10**phi_star_tmp
		end select
		print*,  "aux at rp:x,f=", 2*(phi_star_tmp+spp%mbh_dmless/rp-ex)-(jum*jc)**2/rp**2
	end if
	if(jum>0.99995)then
		do i=2, nbins-1
			aux%fx(i)=(aux%fx(nbins)-aux%fx(1))*(aux%xb(i))/pi*2d0+aux%fx(1)
		end do
	else
		do i=2, nbins-1
			
			sintheta=sin(aux%xb(i))
			costheta=cos(aux%xb(i))
			r=rp+(ra-rp)*sintheta*sintheta
			
			call get_phi_star_full_range(spp,log10(r),phi_star_tmp)
			select case(ctl%ebin_type)
			case(ebin_type_log)
				phi_star_tmp=10**phi_star_tmp
			end select
			vr_tmp=2*(phi_star_tmp+spp%mbh_dmless/r-ex)-(jum*jc)**2/r**2
			!print*, "r,vr_tmp,phi_star_tmp=",r,vr_tmp, phi_star_tmp, &
			!   2*(phi_star_tmp+mbh_dmless/r-ex),(jum*jc)**2/r**2
			if(vr_tmp<=-0.01)then
				print*, "get_aux_function_for_period_pi2: error! what happened??"
				print*, "i,r, vr_tmp,mbh_dmless=",i,r, vr_tmp,spp%mbh_dmless
				print*, "x,r, phi_star_tmp=",aux%xb(i),r, phi_star_tmp	
				print*, "ex,jum=", ex,jum, 2*(phi_star_tmp+spp%mbh_dmless/r-ex),(jum*jc)**2/r**2
				call spp%fma_star%print("fma")
				call spp%fphi_star%print("phi_star")
				print*, "ra,rp=",ra,rp
				print*, "jc=", jc

				call get_phi_star_full_range(spp,log10(ra),phi_star_tmp)
				select case(ctl%ebin_type)
				case(ebin_type_log)
					phi_star_tmp=10**phi_star_tmp
				end select
				vr_tmp=2*(phi_star_tmp+spp%mbh_dmless/ra-ex)-(jum*jc)**2/ra**2
				print*, "vr(ra)=",vr_tmp

				call get_phi_star_full_range(spp,log10(rp),phi_star_tmp)
				select case(ctl%ebin_type)
				case(ebin_type_log)
					phi_star_tmp=10**phi_star_tmp
				end select
				vr_tmp=2*(phi_star_tmp+spp%mbh_dmless/rp-ex)-(jum*jc)**2/rp**2
				print*, "vr(rp)=",vr_tmp

				do j=1,i-1
					print*, "aux:i,x,r,y=",j,aux%xb(j),rp+(ra-rp)*(sin(aux%xb(j))**2),aux%fx(j)
				end do 
				!block
				!	real(8) logrc,rmax, rpx,rax,r_c_iter
				!	integer ier
					logrc=log10(r_c_iter(spp,ex,ier))
					call get_rmax_accurate(spp,dms%fr_phi,log10(ex),rmax)
					call get_rpra_dmless(spp,ex,jum,jc,logrc,rmax,rpx,rax)
					print*, "1:rpx,rpx-rp=",rpx,(rpx-rp)/rp
				!	ctl%debug=1
					call get_rpra_dmless_fast(spp,ex,jum,jc,logrc,rmax,rpx,rax)
					print*, "2:rpx,rpx-rp=",rpx,(rpx-rp)/rp
				!end block
				stop
			end if
			aux%fx(i)=2*(ra-rp)*sintheta*costheta*&
				(abs(vr_tmp))**(-0.5d0)
			!print*, "x,f=",theta,r, phi_star_tmp,ex,jum, 2*(10**phi_star_tmp+1d0/r-ex)-(jum*jc)**2/r**2
			!read(*,*)
			
			if(ieee_is_nan(aux%fx(i)).or.(.not.ieee_is_finite(aux%fx(i))))then
				print*, "get_aux_function_for_period_pi2:error! "
				print*, "x,r, phi_star_tmp=",aux%xb(i),r, phi_star_tmp
				print*, "ex,jum=", ex,jum, 2*(10**phi_star_tmp+spp%mbh_dmless/r-ex),(jum*jc)**2/r**2
				call spp%fma_star%print("fma")
				call spp%fphi_star%print("phi_star")
				print*, "ra,rp=",ra,rp
				print*, "jc=", jc
				stop
			end if
		end do
	end if
	call aux%prepare_spline()
end subroutine 

subroutine get_dphi2_dr2(frho,fma, mbhin, logr,dphidr)
	use com_sts_type
	use constant
	use md_star_pot
	!use model_basic
	implicit none
	type(s1d_type)::frho, fma
	real(8) logr, dphidr, r, rho_tmp, beta_tmp, mbhin

	!call rho%get_value_s(logr,rho_tmp)
	call get_rho_full_range_spp(spp_new,logr, rho_tmp)
	call get_beta_full_range(spp_new,logr, beta_tmp)
	r=10**logr
	dphidr=(2*(mbhin+beta_tmp)-4*pi*r**3*rho_tmp)/r**3
end subroutine

subroutine get_aux_ra(spp,rp,ra,jc,ex,jum,res)
	use com_sts_type
	use model_basic,only:ctl
	use md_star_pot
	use md_coeff
	implicit none
	type(star_pot_para)::spp
	!type(s1d_type)::fma,frho,phi_star
	real(8) res
	real(8) rp,ra, beta_rp, beta_ra, phi_star_tmp,ex
	real(8) jc, jum, Aterm_rp, Aterm_ra 
	real(8) logrp,logra, jump, r ,vr_tmp

	logra=log10(ra)
	call get_beta_full_range(spp,logra, beta_ra)
	
	jump=(jum*jc)**2
	Aterm_ra=-2*(-(spp%mbh_dmless+beta_ra)/ra**2+jump/ra**3)
	if(Aterm_ra<0)then
		print*, "get_aux_2bd:error!, Aterm_ra<0", Aterm_ra
		print*, "ex,jum, jc, ra,rp=",ex,jum,jc, ra,rp
		print*, "M_r_within_max=",spp%M_r_within_max
		call spp%frho_star%print("frho")
		call spp%fma_star%print("fma")
		call get_phi_star_full_range(spp,log10(ra),phi_star_tmp)
		select case(ctl%ebin_type)
		case(ebin_type_log)
			phi_star_tmp=10**phi_star_tmp
		end select
		vr_tmp=2*(phi_star_tmp+spp%mbh_dmless/ra-ex)-(jum*jc)**2/ra**2
		print*, "vr(ra)=",vr_tmp

		call get_phi_star_full_range(spp,log10(rp),phi_star_tmp)
		select case(ctl%ebin_type)
		case(ebin_type_log)
			phi_star_tmp=10**phi_star_tmp
		end select
		vr_tmp=2*(phi_star_tmp+spp%mbh_dmless/rp-ex)-(jum*jc)**2/rp**2
		print*, "vr(rp)=",vr_tmp
		!stop
		res=-1d99
		return
	endif
	res=2*(ra-rp)**0.5*(Aterm_ra)**(-0.5)
	!print*,"res,ra,x=",res, ra,ex,ra*2**0.5/(ex)**0.5

end subroutine

subroutine get_aux_rp(spp,rp, ra,jc,ex,jum,res)
	use com_sts_type
	use model_basic,only:ctl
	use md_coeff
	use md_star_pot
	implicit none
	!type(s1d_type)::fma,frho,phi_star
	type(star_pot_para)::spp
	real(8) res
	real(8) rp,ra, beta_rp, beta_ra, phi_star_tmp,ex
	real(8) jc, jum, Aterm_rp, Aterm_ra 
	real(8) logrp,logra, jump, r ,vr_tmp

	logrp=log10(rp)
	call get_beta_full_range(spp,logrp, beta_rp)

	jump=(jum*jc)**2
	Aterm_rp=2*(-(spp%mbh_dmless+beta_rp)/rp**2+jump/rp**3)
!	print*, "Aterm_rp, Aterm_ra=",Aterm_rp, Aterm_ra
	if(Aterm_rp<0)then
		print*, "get_aux_2bd:error!, Aterm_rp<0", Aterm_rp
		print*, "ex,jum, jc, ra,rp=",ex,jum,jc, ra,rp
		call spp%frho_star%print("frho")
		call spp%fma_star%print("fma")
		call get_phi_star_full_range(spp,log10(ra),phi_star_tmp)
		select case(ctl%ebin_type)
		case(ebin_type_log)
			vr_tmp=2*(10**phi_star_tmp+spp%mbh_dmless/ra-ex)-(jum*jc)**2/ra**2
		case(ebin_type_lin)
			vr_tmp=2*(phi_star_tmp+spp%mbh_dmless/ra-ex)-(jum*jc)**2/ra**2
		end select
		print*, "vr(ra)=",vr_tmp

		call get_phi_star_full_range(spp,log10(rp),phi_star_tmp)
		select case(ctl%ebin_type)
		case(ebin_type_log)
			vr_tmp=2*(10**phi_star_tmp+spp%mbh_dmless/rp-ex)-(jum*jc)**2/rp**2
		case(ebin_type_lin)
			vr_tmp=2*(phi_star_tmp+spp%mbh_dmless/rp-ex)-(jum*jc)**2/rp**2
		end select
		print*, "vr(rp)=",vr_tmp
		stop
	end if
	
	res=2*(ra-rp)**0.5*(Aterm_rp)**(-0.5)
	!print*,"res,rp,x=",res, rp,ex
end subroutine

subroutine get_sample_para_one_no_pd(dm,sp,spp)
	 use com_main_gw
	 implicit none
	 type(particle_sample_type)::sp
	 type(diffuse_mspec)::dm
	 type(star_pot_para)::spp
	 real(8) ex, logex, jc, jc_dmless
	 real(8) rmax,r_c, rc,jm,jc_xy,rp_xy,ra_xy
	 real(8) pd_xy,p_EJ_dmless,r_c_iter
	 integer ier

	 ex=sp%x
	 logex=log10(ex)
	 sp%en=sp%x*ctl%energy0
	 call get_rmax_accurate(spp ,  dm%fr_phi, logex,rmax)
	 !print*, "rmax=",rmax
	 !call get_rmax_accurate(dm%fphi_star,  dm%fr_phi, logex,rmax)
	 !print*, "rmax=",rmax
	 !read(*,*)
	 !rc=r_c(spp,ex,ier)
	 rc=r_c_iter(spp,ex,ier)
	 
	 if(ier.eq.-99)then
		print*, "sp%m=", sp%m
	 end if
	 jc_xy=jc_dmless(rc,spp)
	 sp%jc=jc_xy*ctl%v0*r0_cl
	 !===========================
	 sp%jm=sp%jph/sp%jc
	 call set_jm_bound(sp%jm)
	 !===========================
	 jm=sp%jm
	 call get_jm_idx(sp%jm, sample_table_idy,sample_table_rdy,sample_evjum)
	 sp%jph=jm*sp%jc

	 !if(sp%jm>0.999988.and.rc<1.62d-4.and.abs(ex-3104.8909)<0.001)then
	!	block
	!		real(8) fout1,fout2, logrc
	!		real(8) phifun,r,par(50),betafun
	!	print*, "rc,ex,jm=",rc,ex,jm
	!	logrc=log10(rc)
	!	call get_r_root(spp,logrc,ex,(jm*jc_xy)**2,fout1)
	!	
	!	!call get_alpha_full_range(fphi,x,alphafun)
	!	call get_phi_star_full_range(spp,logrc,phifun)
	!	call get_beta_full_range(spp,logrc,betafun)
!
	!	!select case(ctl%ebin_type)
	!	!case(ebin_type_log)
	!	fout2=2*ex-(spp%mbh_dmless/rc+2*10**phifun-betafun/rc)
!
	!	print*, "fout1,fout2=",fout1,fout2
!
!
	!	end block
	 !end if

	 if(jc_xy.eq.0)then
		 sp%rp=10**dms%logrmin*r0_cl
		 sp%ra=sp%rp
		 sp%period=0d0
		 print*, "jc_xy=0", sp%x, sp%jm, sp%jph
		 return
	 end if
	 call get_rpra_dmless(spp, ex, jm, jc_xy, &
				 log10(rc), rmax, rp_xy,ra_xy)
	 sp%rp=rp_xy*r0_cl       
	 sp%ra=ra_xy*r0_cl
	 !print*, "ex,jm,rp_xy,ra_xy,jc_xy=",ex,jm,rp_xy,ra_xy,jc_xy,rc
	 !pd_xy=p_EJ_dmless(spp, ex,jm,  jc_xy, rp_xy,ra_xy)
	 !sp%period=pd_xy*r0_cl/ctl%v0
	 

 end subroutine