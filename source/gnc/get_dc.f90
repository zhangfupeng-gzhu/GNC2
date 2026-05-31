
subroutine dm_get_dc_mpi(dm)
    use com_main_gw
    implicit none
    type(diffuse_mspec):: dm
    interface 
		subroutine mb_get_dc_mpi_starpt_gx_iregular(mb,spp)
			use com_main_gw
			implicit none
			type(mass_bins),target::mb
		!	type(s1d_type)::phi_star
			type(star_pot_para)::spp
		end subroutine 
    end interface
    integer i

	
	do i=1, dm%n
		if(rid.eq.0)then
			print*, "m(i)=",i,dm%mb(i)%df_coe_bins
		end if
		
		if(ctl%dc_grid_type.eq.dc_grid_irregular)then
			if(rid.eq.0)then
				print*, "dc grid type=dc_grid_irregular"
			end if
			!print*, '1'
			!call prepare_barge_ir_tables()
			!print*, "2"
			call set_dc_bins_ir(dm)
			!print*, "3"
			select case(ctl%fden_ana_est_method)
			case(fden_ana_est_method_1d_iso)
				call mb_get_dc_mpi_starpt_gx_iregular(dm%mb(i),spp_new)
			case(fden_ana_est_method_2d)
				call mb_get_dc_mpi_starpt_gx_iregular_2d(dm%mb(i),spp_new)
			end select
		endif
	end do
	 
    call get_dc0(dm)
     
end subroutine
subroutine set_dc_bins_ir(dm)
	use com_main_gw
	implicit none
	type(diffuse_mspec)::dm
	integer i
	do i=1, dm%n
		dm%mb(i)%dc%s2_de_110%xcenter=dms%dlxb_ir%xb
		dm%mb(i)%dc%s2_de_0%xcenter=dms%dlxb_ir%xb
		dm%mb(i)%dc%s2_dee%xcenter=dms%dlxb_ir%xb
		dm%mb(i)%dc%s2_dj_111%xcenter=dms%dlxb_ir%xb
		dm%mb(i)%dc%s2_dj_rest%xcenter=dms%dlxb_ir%xb
		dm%mb(i)%dc%s2_djj%xcenter=dms%dlxb_ir%xb
		dm%mb(i)%dc%s2_dej%xcenter=dms%dlxb_ir%xb
	end do
	dm%dc0%s2_de_110%xcenter=dms%dlxb_ir%xb
	dm%dc0%s2_de_0%xcenter=dms%dlxb_ir%xb
	dm%dc0%s2_dee%xcenter=dms%dlxb_ir%xb
	dm%dc0%s2_dj_111%xcenter=dms%dlxb_ir%xb
	dm%dc0%s2_dj_rest%xcenter=dms%dlxb_ir%xb
	dm%dc0%s2_djj%xcenter=dms%dlxb_ir%xb
	dm%dc0%s2_dej%xcenter=dms%dlxb_ir%xb

end subroutine 
subroutine get_coeff_ej_spt_y_gx_iregular( xb, yb,n, spp,barge,asymp, &
	e_110,e_0,ee,j_111,j_rest,jj,ej)
	use com_main_gw
	implicit none
	!type(s1d_type)::phi_star
	type(star_pot_para)::spp
	integer n,i,j
	type(s1d_ird_type),target::barge
	type(s2d_type)::gxj_ir
	real(8) asymp,energyx,jum
	integer ierr
	integer nbg, ned, nblock, ntasks
	external::fgx_mb_ir
	real(8) xb(n),yb(n)
	real(8) e_110(n,n),e_0(n,n),ee(n,n) 
	real(8) j_111(n,n),j_rest(n,n),jj(n,n),ej(n,n)
	real(8) f0
	real(8) t1, t2,t3,t4

	fc_ir_share=>barge

	fgx_g0=asymp
	!emin_dstr_factor=10**fc_share%xb(1)
	!emax_dstr_factor=10**fc_share%xb(fc_share%nbin)
	nbg=ctl%nblock_mpi_bg
	ned=ctl%nblock_mpi_ed
	nblock=ctl%nblock_size
	ntasks=ctl%ntasks 
	call cpu_time(t1) 

	do i=nbg, ned, ntasks
		!i=46
	!i=n-3;j=n-3
		energyx=xb(i)			
		!call cpu_time(t3)		
		
		!energyx=xb(n)

		call get_sigma0(energyx, fgx_mb_ir, f0)
		do j=1, n
			jum=yb(j) 
			
			call get_coeff_ffuncs_cfs_grid(energyx,jum, fgx_mb_ir,spp, &
			dms%rp%fxy(i,j),dms%ra%fxy(i,j),dms%pd%fxy(i,j),dms%jc%fx(i),&
			f0, e_110(i,j), e_0(i,j),ee(i,j), ej(i,j), jj(i,j), j_rest(i,j), j_111(i,j))  
		end do 
	end do 
	call cpu_time(t2)
	if(rid.eq.0)then
		print*, "timetakes,ibg, ied=",  t2-t1, nbg, ned
	end if
	call mpi_BARRIER(mpi_comm_world, ierr)
	call collect_data_mpi_y(e_110, n,nbg, ned, nblock, ntasks)
	 
	call mpi_BARRIER(mpi_comm_world, ierr)
	call collect_data_mpi_y(e_0, n,nbg, ned, nblock, ntasks)
	call mpi_BARRIER(mpi_comm_world, ierr)
	call collect_data_mpi_y(ee, n,nbg, ned, nblock, ntasks)
	call mpi_BARRIER(mpi_comm_world, ierr)
	call collect_data_mpi_y(ej, n,nbg, ned, nblock, ntasks)
	call mpi_BARRIER(mpi_comm_world, ierr)
	call collect_data_mpi_y(jj, n,nbg, ned, nblock, ntasks)
	call mpi_BARRIER(mpi_comm_world, ierr)
	call collect_data_mpi_y(j_rest, n,nbg, ned, nblock, ntasks)
	call mpi_BARRIER(mpi_comm_world, ierr)
	call collect_data_mpi_y(j_111, n,nbg, ned, nblock, ntasks)
end subroutine
 

subroutine get_coeff_ffuncs_cfs_grid(ex,jum, fgx,spp, rp,ra,pd,jc,&
	f0,e110,e0,ee,ej,jj, jrest,j111)
	use model_basic
	use md_star_pot
	use, intrinsic :: ieee_arithmetic
	use MPI_comu, only:rid
	implicit none
	real(8) ex,jum,e110,e0,ee,ej,jj, jrest,j111
	real(8) f0,emax
	external:: fgx
	type(s1d_type)::f1,f2, ss !, aux
	!type(s1d_type)::fphi,fma, frho
	type(star_pot_para)::spp
	real(8) rp,ra,pd,jc, lrp, lra
	integer nintn,i,idid 
	!real(8) t1, t2, t3
	!ctl%debug=1
	!print*, "emin_factor,emax_factor=",emin_factor,emax_factor
	!read(*,*)
	!call cpu_time(t1)
	if(rp.eq.ra)then
		!print*, "rp=ra=",rp,ra
		call get_coeff_ffuncs_cfs_grid_circular(ex, fgx,spp, rp,&
		f0,e110,e0,ee,ej,jj, jrest,j111)
		return
	end if
	  
	call get_aux_function_for_period_pi2(common_aux, spp,ex,jum,jc,rp,ra)
	!call common_aux%prepare_spline()

	if(jum<0.1d0)then
	nintn=f12_function_bin_size_1
	else
	nintn=f12_function_bin_size_2
	end if 
	lrp=log10(rp)
	lra=log10(ra)
	!
	call f1%init(lrp, lra, nintn, sts_type_grid)
	call f1%set_range()
	call f2%init( lrp, lra,nintn, sts_type_grid)
	call f2%set_range()
	 
	emax=10**fc_ir_share%xmax !min(emax_factor,dms%emax) 
	select case(ctl%get_dc_method)			
	case(get_dc_method_fast)
		do i=2, nintn-1
 
			call get_ffuncs_series_fast(ex,  spp_new%mbh_dmless, f1%xb(i), jum, jc, ra,rp, spp, fgx, f1%fx(i), 0.5d0,.false.)
			call get_ffuncs_series_fast(ex,  spp_new%mbh_dmless, f2%xb(i), jum, jc, ra,rp, spp, fgx, f2%fx(i), 1.5d0,.false.)
		end do
		call get_ffuncs_series_fast(ex, spp_new%mbh_dmless, lrp,  jum, jc, ra,rp, spp, fgx, f1%fx(1), 0.5d0,.true.)
		call get_ffuncs_series_fast(ex, spp_new%mbh_dmless, lra,  jum, jc, ra,rp, spp, fgx, f1%fx(nintn), 0.5d0,.true.)
		call get_ffuncs_series_fast(ex, spp_new%mbh_dmless, lrp,  jum, jc, ra,rp, spp, fgx, f2%fx(1), 1.5d0,.true.)
		call get_ffuncs_series_fast(ex, spp_new%mbh_dmless, lra,  jum, jc, ra,rp, spp, fgx, f2%fx(nintn), 1.5d0,.true.)
	case(get_dc_method_general)
		do i=2, nintn-1
			call get_ffuncs_series(ex, emax, emin_factor, spp_new%mbh_dmless, f1%xb(i), ra,rp, spp, fgx, f1%fx(i), 0.5d0)
			!call get_ffuncs_series(ex, mbh_dmless, f1%xb(i), ra,rp, fphi, fgx, f1%fx(i), 0.5d0)
			call get_ffuncs_series(ex, emax, emin_factor, spp_new%mbh_dmless, f2%xb(i), ra,rp, spp, fgx, f2%fx(i), 1.5d0)
			!call get_ffuncs_series(ex,  mbh_dmless, f2%xb(i), ra,rp, fphi, fgx, f2%fx(i), 1.5d0)
		 end do
		 
		 call get_ffuncs_series_rarp(ex, emax, emin_factor, lrp,  jum, jc, fgx, f1%fx(1), 0.5d0)
		 call get_ffuncs_series_rarp(ex, emax, emin_factor, lra,  jum, jc, fgx, f1%fx(nintn), 0.5d0)
		 !
		 call get_ffuncs_series_rarp(ex, emax, emin_factor, lrp,  jum, jc, fgx, f2%fx(1), 1.5d0)
		 call get_ffuncs_series_rarp(ex, emax, emin_factor, lra,  jum, jc, fgx, f2%fx(nintn), 1.5d0)	
 
	end select
	 

	if(ctl%chattery.ge.3)then
	print*, "ex,jum,log10(ex),log10(jum), jc, rp, ra=",ex,jum,log10(ex),log10(jum), jc, rp, ra
	call common_aux%print("aux")
	call f1%print("f1")
	call f2%print("f2")
	!read(*,*)
	end if 
	e110=0 

	call my_integral_acc(0d0,pi/2d0, e110, star_dc_int_acc_a,star_dc_int_acc_r, fcn_e110,idid)	
	if(idid<0)then
	print*, "stop, e110 problem, rid=", rid
	print*, "ex,jum,log10(ex),log10(jum), jc, rp, ra=",ex,jum,log10(ex),log10(jum), jc, rp, ra
	call common_aux%print("aux")
	call f1%print("f1")
	call f2%print("f2")
	stop
	end if
	e110=e110/pd*2
	e0=-f0
	ee=0
	 
	call my_integral_acc(1d-6,pi/2d0, ee,star_dc_int_acc_a,star_dc_int_acc_r, fcn_ee,idid)	
	if(idid<0)then
	print*, "stop, ee problem, rid=", rid
	print*, "ex,jum,log10(ex),log10(jum), jc, rp, ra=",ex,jum,log10(ex),log10(jum), jc, rp, ra
	call common_aux%print("aux")
	call f1%print("f1")
	call f2%print("f2")

	call check_aux(common_aux,spp,ex, jum,jc, rp,ra)
	call spp%fphi_star%print("fphi")
	call spp%fma_star%print("fma")
	call dms%mb(1)%all%barge%print("1:all%barge")

	stop
	end if 
	ee=abs(ee*8d0/3d0/pd)
	!print*, "ex,jum,ee,pd=", ex,jum,ee,pd
	!read(*,*)
	j111=0
	call my_integral_acc(0d0,pi/2d0, j111,star_dc_int_acc_a,star_dc_int_acc_r, fcn_j111,idid)	
	if(idid<0)then
	print*, "stop, j111 problem, rid=", rid
	print*, "ex,jum,log10(ex),log10(jum), jc, rp, ra=",ex,jum,log10(ex),log10(jum), jc, rp, ra
	call common_aux%print("aux")
	call f1%print("f1")
	call f2%print("f2")
	stop
	end if	
	j111=j111/pd
	!print*, "j111=", j111
	!stop
	jrest=0
	call my_integral_acc(0d0,pi/2d0, jrest, star_dc_int_acc_a,star_dc_int_acc_r, fcn_jrest,idid)		
	if(idid<0)then
	print*, "stop, jrest problem, rid=", rid
	print*, "ex,jum,log10(ex),log10(jum), jc, rp, ra=",ex,jum,log10(ex),log10(jum), jc, rp, ra
	call common_aux%print("aux")
	call f1%print("f1")
	call f2%print("f2")
	stop
	end if
	jrest=jrest/pd
	!print*, "jrest=", jrest
	jj=0
	!call my_integral_acc(0d0,pi/2d0, jj,1d-12,1d-9, fcn_jj,idid)		
	call my_integral_acc(0d0,pi/2d0, jj, star_dc_int_acc_a,star_dc_int_acc_r,fcn_jj,idid)		
	jj=abs(jj*2/pd)
	if(idid<0)then
	print*, "stop, jj problem, rid=", rid
	print*, "ex,jum,log10(ex),log10(jum), jc, rp, ra=",ex,jum,log10(ex),log10(jum), jc, rp, ra
	call common_aux%print("aux")
	call f1%print("f1")
	call f2%print("f2")
	stop
	end if
	!print*, "jj=", jj
	ej=0
	call my_integral_acc(0d0,pi/2d0, ej, star_dc_int_acc_a,star_dc_int_acc_r, fcn_ej,idid)		
	ej=-ej*4/pd/3d0*jum
	if(idid<0)then
	print*, "stop, ej problem, rid=", rid
	print*, "ex,jum,log10(ex),log10(jum), jc, rp, ra=",ex,jum,log10(ex),log10(jum), jc, rp, ra
	call common_aux%print("aux")
	call f1%print("f1")
	call f2%print("f2")
	stop
	end if 
contains
	subroutine fcn_e110(n, x, y, f, par, ipar)
	implicit none
	integer n, ipar(100)
	real(8) x, y(n), f(n), par(100),vr,vr2
	real(8) f1inp,r, aux_tmp
 
	r=rp+(ra-rp)*sin(x)**2
	
	call common_aux%get_value_l(x, aux_tmp)
	call f1%get_value_l(log10(r), f1inp) 
	f1inp=10**f1inp

	f(1)=f1inp*aux_tmp
 
	if(ieee_is_nan(f1inp).or.ieee_is_nan(f(1)).or.(.not.ieee_is_finite(f(1))))then
		print*, "nan detected,  f1inp=", f1inp 
		call f1%print("f1")
		call f2%print("f2")
		call common_aux%print("aux")
		print*, "ex, jm=", ex, jum
		print*, "x, aux_tmp=",x,aux_tmp
		stop
	end if
	end subroutine
	subroutine fcn_ee(n, x, y, f, par, ipar)
	use, intrinsic :: ieee_arithmetic
	implicit none
	integer n, ipar(100)
	real(8) x, y(n), f(n), par(100),vr,v,vr2
	real(8) f2inp,fphiinp,r, aux_tmp,Ev

	r=rp+(ra-rp)*sin(x)**2
	
	call common_aux%get_value_l(x, aux_tmp)
	call f2%get_value_l(log10(r), f2inp)
	!call splint_mylib(f2%xb,f2%fx,f2%y2a,f2%nbin, x, f2inp)
	f2inp=10**f2inp
	!call get_phi_star_full_range(fphi,log10(r),fphiinp)
	!fphiinp=10**fphiinp
	!fphiinp=0d0
	if(ieee_is_nan(fphiinp).or.ieee_is_nan(f2inp))then
		print*, "nan detected, fphiinp, f2inp=", fphiinp, f2inp 
	end if
	 
	if(aux_tmp.ne.0)then
		vr=(2*(ra-rp)*sin(x)*cos(x))/aux_tmp
		Ev=(vr*vr+(jum*jc)*(jum*jc)/(r*r))/(2*ex)
		!block
		!	real(8) vr
		!	!call aux%print("aux")
		!	vr=(aux_tmp/(2*(ra-rp)*sin(x)*cos(x)))**(-2)
		!	print*, "vr=",vr
			
		!read(*,*)
		!end block
		f(1)=(f0+f2inp)*Ev*aux_tmp
	else
		print*, "ee:aux_tmp=0"
		f(1)=0d0
		!stop
	end if
	 
	if(ieee_is_nan(fphiinp).or.ieee_is_nan(f2inp).or.ieee_is_nan(f(1)).or.(.not.ieee_is_finite(f(1))))then
		print*, "nan detected, fphiinp, f1inp=", fphiinp, f2inp 
		call f1%print("f1")
		call f2%print("f2")
		call common_aux%print("aux")
		print*, "ex, jm=", ex, jum
		print*, "x, aux_tmp=",x,aux_tmp
		stop
	end if
	
	end subroutine

	subroutine fcn_j111(n, x, y, f, par, ipar)
	use, intrinsic :: ieee_arithmetic
	implicit none
	integer n, ipar(100)
	real(8) x, y(n), f(n), par(100),v
	real(8) f1inp,fphiinp,vr,r, aux_tmp, Ev
	
	
	r=rp+(ra-rp)*sin(x)**2
	
	call common_aux%get_value_l(x, aux_tmp)
	call f1%get_value_l(log10(r), f1inp)

	!call splint_mylib(f1%xb,f1%fx,f1%y2a,f1%nbin, x, f1inp)
	f1inp=10**f1inp
	!call get_phi_star_full_range(fphi,log10(r),fphiinp)
	!fphiinp=10**fphiinp
	if(ieee_is_nan(fphiinp).or.ieee_is_nan(f1inp))then
		print*, "nan detected, fphiinp, f1inp=", fphiinp, f1inp 
	end if
	!if(x>pi/2d0-0.01d0)then
	!	Ev=(jum*jc)**2/(2*ex*ra**2)
	!else
	!	Ev=(fphiinp+mbh_dmless/r)/ex-1
	!end if
	vr=2*(ra-rp)*sin(x)*cos(x)/aux_tmp
	Ev=(vr*vr+(jum*jc)*(jum*jc)/(r*r))/(2*ex)
	!print*, "x, Ev=",x,Ev
	f(1)=-f1inp*jum/Ev*aux_tmp*2
	
	if(Ev<=0d0)then
		print*, "error in j111!,rid, Ev=", rid,Ev
		print*, "fphiinp, r, ex=", fphiinp, r, ex, rp, ra
		call f1%print("f1")
		call f2%print("f2")
		call common_aux%print("aux")
		print*, "ex, jm=", ex, jum
		print*, "x, aux_tmp,f(1)=",x,aux_tmp,f(1)
		print*, "x>pi/2d0-0.01d0",x>pi/2d0-0.01d0
		print*, (jum*jc)**2/(2*ex*ra**2)
		stop
	end if
	if(ieee_is_nan(fphiinp).or.ieee_is_nan(f1inp).or.ieee_is_nan(f(1)).or.(.not.ieee_is_finite(f(1))))then
		print*, "nan detected, fphiinp, f1inp=", fphiinp, f1inp 
		call f1%print("f1")
		call f2%print("f2")
		call common_aux%print("aux")
		print*, "ex, jm=", ex, jum
		print*, "x, aux_tmp=",x,aux_tmp
		stop
	end if
	end subroutine
	subroutine fcn_jrest(n, x, y, f, par, ipar)
	use, intrinsic :: ieee_arithmetic
	implicit none
	integer n, ipar(100)
	real(8) x, y(n), f(n), par(100),vr,Ev
	real(8) f1inp,fphiinp,f2inp,vr2,r, aux_tmp

	!call linear_int(f1%xb,f1%fx,f1%nbin, log10(x),f1inp)
	!if(f1inp<0) f1inp=0
	r=rp+(ra-rp)*sin(x)**2
	
	call common_aux%get_value_l(x, aux_tmp)
	call f1%get_value_l(log10(r), f1inp)

	!call splint_mylib(f1%xb,f1%fx,f1%y2a,f1%nbin, x, f1inp)
	f1inp=10**f1inp
	!call splint_mylib(f2%xb,f2%fx,f2%y2a,f2%nbin, log10(r), f2inp)
	call f2%get_value_l(log10(r), f2inp)
	f2inp=10**f2inp 
	if(ieee_is_nan(f2inp).or.ieee_is_nan(f1inp))then
		print*, "nan detected, f2inp, f1inp, fphiinp=", f2inp, f1inp
	end if
	 
	f(1)=(r**2*ex/(jum*jc**2)*(f1inp-f2inp/3d0+2d0/3d0*f0))*aux_tmp
 
	if(ieee_is_nan(f1inp).or.ieee_is_nan(f(1)).or.(.not.ieee_is_finite(f(1))))then
		print*, "nan detected,  f1inp=",  f1inp 
		call f1%print("f1")
		call f2%print("f2")
		call common_aux%print("aux")
		print*, "ex, jm=", ex, jum
		print*, "x, aux_tmp=",x,aux_tmp
		stop
	end if
	end subroutine
	subroutine fcn_jj(n, x, y, f, par, ipar)
	use, intrinsic :: ieee_arithmetic
	implicit none
	integer n, ipar(100)
	real(8) x, y(n), f(n), par(100),vr,v, Ev
	real(8) f1inp,f2inp,juu, vr2,r, aux_tmp, rjtmp

	r=rp+(ra-rp)*sin(x)**2
	
	call common_aux%get_value_l(x, aux_tmp)
	call f1%get_value_l(log10(r), f1inp)
	f1inp=10**f1inp
	call f2%get_value_l(log10(r), f2inp)
	f2inp=10**f2inp
	! call get_phi_star_full_range(spp,log10(r),fphiinp)
	! !fphiinp=10**fphiinp
	! select case(ctl%ebin_type)
	! case(ebin_type_log)
	! 	fphiinp=10**fphiinp
	! !case(ebin_type_lin)
	! end select
	if(ieee_is_nan(f2inp).or.ieee_is_nan(f1inp))then
		print*, "nan detected, f2inp, f1inp=", f2inp, f1inp
	end if
	 
	vr=2*(ra-rp)*sin(x)*cos(x)/aux_tmp
	Ev=(vr*vr+(jum*jc)*(jum*jc)/(r*r))/(2*ex)
	!print*, "ev=",ev
	!read(*,*)

	rjtmp=(r**2*ex/jc**2-jum**2/2d0/Ev)
	f(1)=(jum**2/3d0/Ev*f2inp+rjtmp*(f1inp-f2inp/3d0)&
		+2d0/3d0*r**2*ex/jc**2*f0)*aux_tmp
	if(rjtmp<-0.01)then
		print*, "???? Ev=", Ev
		print*, "r^2x/jc^2-j^2/2Ev=", r**2*ex/jc**2-jum**2/2d0/Ev
		print*, "r,ex,jum=", r,ex,jum
		! print*, "fphiinp, f1inp=", fphiinp, f1inp 
		call f1%print("f1")
		call f2%print("f2")
		call common_aux%print("aux")
		print*, "ra, rp=",ra,rp
		print*, "x, aux_tmp=",x,aux_tmp
		print*, "x>pi/2d0-0.01d0", x>pi/2d0-0.01d0
		stop
	end if


	if(ieee_is_nan(f1inp).or.ieee_is_nan(f(1)).or.(.not.ieee_is_finite(f(1))))then
		! print*, "nan detected, fphiinp, f1inp=", fphiinp, f1inp 
		call f1%print("f1")
		call f2%print("f2")
		call common_aux%print("aux")
		print*, "ex, jm=", ex, jum
		print*, "x, aux_tmp=",x,aux_tmp
		stop
	end if
	end subroutine
	subroutine fcn_ej(n, x, y, f, par, ipar)
	use, intrinsic :: ieee_arithmetic
	implicit none
	integer n, ipar(100)
	real(8) x, y(n), f(n), par(100),vr,v,vr2
	real(8) f1inp,fphiinp,f2inp,juu,r, aux_tmp

	!call linear_int(f1%xb,f1%fx,f1%nbin, log10(x),f1inp)
	!if(f1inp<0) f1inp=0
	r=rp+(ra-rp)*sin(x)**2
	
	call common_aux%get_value_l(x, aux_tmp)

	!call splint_mylib(f1%xb,f1%fx,f1%y2a,f1%nbin, log10(x), f1inp)
	call f2%get_value_l(log10(r), f2inp)
	!call splint_mylib(f2%xb,f2%fx,f2%y2a,f2%nbin, x, f2inp)
	f2inp=10**f2inp
	 
		f(1)=(f2inp+f0)*aux_tmp
	!end if
	if(ieee_is_nan(f2inp).or.ieee_is_nan(f(1)).or.(.not.ieee_is_finite(f(1))))then
		print*, "nan detected, fphiinp, f1inp=", fphiinp, f1inp 
		call f1%print("f1")
		call f2%print("f2")
		call common_aux%print("aux")
		print*, "ex, jm=", ex, jum
		print*, "x, aux_tmp=",x,aux_tmp
		stop
	end if
	end subroutine
end subroutine


subroutine get_coeff_ffuncs_cfs_grid_circular(ex,fgx,spp, rc,&
	f0,e110,e0,ee,ej,jj, jrest,j111)
use model_basic
use, intrinsic :: ieee_arithmetic
use MPI_comu, only:rid
use md_star_pot
implicit none
real(8) ex,jum,e110,e0,ee,ej,jj, jrest,j111
real(8) f0
external:: fgx
real(8)::f1,f2, ss !, aux
!type(s1d_type)::fphi,fma, frho
type(star_pot_para)::spp
real(8) rc,pd,jc, vc2, lrc,phi_out_tmp
integer nintn,i,idid 
!if(rid.eq.0)then
!	print*, "get_coeff_ffuncs_cfs_grid_circular"
!	stop
!end if
nintn=f12_function_bin_size_1

lrc=log10(rc)

call get_ffuncs_series(ex, emax_factor, emin_factor, spp_new%mbh_dmless, lrc, rc,rc, spp, fgx, f1, 0.5d0)
call get_ffuncs_series(ex, emax_factor, emin_factor, spp_new%mbh_dmless, lrc, rc,rc, spp, fgx, f2, 1.5d0)

e110=f1

e0=-f0
call get_phi_star_full_range(spp,lrc,phi_out_tmp)
vc2=2*(10**phi_out_tmp+spp_new%mbh_dmless/rc-ex)
ee=abs((vc2/ex*2d0/3d0)*(f2+f0))
!print*, "ex,jum,ee,pd=", ex,jum,ee,pd
!read(*,*)
j111=ex/vc2*f1
jrest=ex/2d0/vc2*(f1-f2/3d0+2d0/3d0*f0)
!print*, "jrest=", jrest
jj=abs(2*ex/3d0/vc2*(f2+f0))
ej=-2d0/3d0*(f2+f0)
end subroutine

subroutine get_ffuncs_series_fast(ex,  mbhin,logr, jum, jc, ra,rp, spp, func, res, ipow,edge)
	use, intrinsic :: ieee_arithmetic
	use com_sts_type
	use my_intgl
	use md_coeff
	use model_basic,only:dms,fc_share,fc_ir_share,ctl,sample_logemin,sample_logemax
	use MPI_comu,only:rid
	use md_star_pot
	implicit none
	external::func
	!type(s1d_type)::fphi
	logical edge
	type(star_pot_para)::spp
	real(8) ex,exmax,exmin,logr,fout,ratio,res,ra,rp, jum, jc
	real(8) mbhin,tmin,tmax,fx_tmp,xmin,xmax,ipow
	real(8) logx0, logx1,logxmid, logxminmid,logxmaxmid,t0,t1,phi,logxmin,logxmax
	real(8) line_slope, line_c,yp0,yp1,xp0,xp1,PX,QX, res_tmp,yout,dint 
	integer debug, nbin,flag_end
	integer i,idx0,idx1,idxmin,idxmax
	!ctl%debug=1
	if(.not.edge)then
		call get_phi_star_full_range(spp,logr, yout)
		phi=(10**yout+mbhin/10**logr)
		ratio=phi/(ex)
	else
		ratio=(jum*jc)**2/(10**(logr*2))/2d0/ex+1
	end if

	if(ratio<1)then
		print*, "get_ffuncs_series: error! ratio<1,ratio=", ratio
		select case(ctl%ebin_type)
		case(ebin_type_log)
			print*, "phi_star, r, ex=", 10**yout,10**logr, ex  
		end select
		stop
	end if
	!=======test
	!fc_ir_share%fx=10**(fc_ir_share%xb*2)
	!==========
	nbin=dms%barp_ir%nbin
	tmin=min(max((ratio-10**sample_logemax/ex)/(ratio-1)*(1d0+1d-11),0d0),1d0)
	tmax=min((ratio-10**sample_logemin/ex)/(ratio-1),1d0)
	!print*, "tmin,tmax=",tmin,tmax, 10**sample_logemax,ex
	if(tmax-tmin.le.0)then
		res=-100
		return
	end if
	logxmin=max(log10((ratio-tmax*(ratio-1))*ex),sample_logemin)
	logxmax=min(log10((ratio-tmin*(ratio-1))*ex),sample_logemax)

	if(logxmin>logxmax)then
		res=-100
		return
	end if

	call return_idx_ir(fc_ir_share%xb,nbin,logxmin,idxmin,sts_type_dstr)
	call return_idx_ir(fc_ir_share%xb,nbin,logxmax,idxmax,sts_type_dstr)
	if(idxmin<0.or.idxmax<0)then
		print*, "error!"
		print*, "logxmin,logxmax, sample_logemin,sample_logemax, nbin"
		print*, logxmin,logxmax, sample_logemin,sample_logemax, nbin
		print*, "idxmin,idxmax=",idxmin,idxmax
		print*, "tmin,tmax,tmax-tmin=",tmin,tmax,tmax-tmin
		stop
	end if
	!if(ctl%debug.eq.1)then
	!	print*, 'ex,phi=',ex,phi, log10(ex),log10(phi)
	!	print*, "logxmin,logxmax, sample_logemin,sample_logemax, nbin"
	!	print*, logxmin,logxmax, sample_logemin,sample_logemax, nbin
	!	print*, "idxmin,idxmax=",idxmin,idxmax
	!endif
	logxminmid=fc_ir_share%xb(idxmin)
	logxmaxmid=fc_ir_share%xb(idxmax)
	!print*, "logxminmid,logxmaxmid=", logxminmid,logxmaxmid
	logx0=logxmin
	res_tmp=0
	flag_end=0
100	if(logx0<logxminmid)then 
		!if(ctl%debug.eq.1)then
		!	print*, "=logx0<logxminmid", logx0, logxminmid
		!end if
		if(idxmin.ne.1)then
			idx1=idxmin
			idx0=idxmin-1
			logxmid=fc_ir_share%xb(idx1)
		else
			idx1=2
			idx0=1
			logxmid=fc_ir_share%xb(1)
		end if
		if(logxmid<logxmax)then
			logx1=logxmid
		else
			logx1=logxmax
			flag_end=1
		end if
	else
		!if(ctl%debug.eq.1)then
		!	print*, "=logx0<logxminmid", logx0, logxminmid
		!end if
		if(idxmin.ne.nbin)then
			idx1=idxmin+1
			idx0=idxmin
			logxmid=fc_ir_share%xb(idx1)
			if(logxmid<logxmax)then
				logx1=logxmid
			else
				logx1=logxmax
				flag_end=1
			end if
		else
			idx1=nbin
			idx0=nbin-1
			logxmid=fc_ir_share%xb(nbin)
			logx1=logxmax
			flag_end=1
		end if
		
	end if
	!if(ctl%debug.eq.1)then
	!	print*, "xmid,logxmaxmid",logxmid,logxmaxmid,logxmax
	!	print*, "x0,x1,flag_end=",logx0,logx1,flag_end
	!end if
	yp1=fc_ir_share%fx(idx1)
	yp0=fc_ir_share%fx(idx0)
	xp1=10**fc_ir_share%xb(idx1)
	xp0=10**fc_ir_share%xb(idx0)
	!if(ctl%debug.eq.1)then
	!	print*, "xp0,xp1,yp0,yp1=", xp0,xp1,yp0,yp1,idx0,idx1
	!end if
	line_slope=(yp1-yp0)/(xp1-xp0)
	line_c=yp0-line_slope*xp0
	!print*, "line_slope,c=",line_slope,line_c

	if(flag_end.eq.1)then
		if(10**logx1>-line_c/line_slope.and.line_slope.lt.0d0)then
			if(-line_c/line_slope>0)then
				logx1=log10(-line_c/line_slope)
			else
				print*, "error! x>-c/slope, -line_c/line_slope<0??"
				stop
			end if
		end if
		if(10**logx0>-line_c/line_slope.and.line_slope.lt.0d0)then
			if(-line_c/line_slope>0)then
				logx0=log10(-line_c/line_slope)
			else
				print*, "error! x>-c/slope, -line_c/line_slope<0??"
				stop
			end if
		end if
		
	end if
	if(idx0.eq.1)then
		if(10**logx0<-line_c/line_slope.and.line_slope.gt.0d0)then
			if(-line_c/line_slope>0)then
				logx0=log10(-line_c/line_slope)
			else
				print*, "error! x<-c/slope, -line_c/line_slope<0??"
				stop
			end if
		end if
		if(10**logx1<-line_c/line_slope.and.line_slope.gt.0d0)then
			if(-line_c/line_slope>0)then
				logx1=log10(-line_c/line_slope)
			else
				print*, "error! x<-c/slope, -line_c/line_slope<0??"
				stop
			end if
		end if
	end if
	!if(ratio-1<1d-5)then
	px=1d0/(ipow+1)*(line_slope* ratio* ex+line_c)*(ratio-1)
	qx=1d0/(ipow+2)*(ratio-1)*(ratio-1)*line_slope*ex
	t0=(ratio-10**logx1/ex)/(ratio-1)
	t1=(ratio-10**logx0/ex)/(ratio-1)
	
	if(t0<0)then
		t0=0
	end if
	if(t1<0)then
		t1=0
	end if
	
	dint=px*(t1**(ipow+1)-t0**(ipow+1))-qx*(t1**(ipow+2)-t0**(ipow+2))! 
	if(dint<0)then
		dint=0
	end if
	res_tmp=res_tmp+dint
	!if(ctl%debug.eq.1)then
	!	print*, "res_tmp=",res_tmp
	!end if
	if(flag_end.ne.1)then
		logx0=logx1
		if(idxmin<nbin)then
			idxmin=idxmin+1
			logxminmid=fc_ir_share%xb(idxmin)
		end if
		goto 100
	end if
	if(res_tmp.le.0)	then
		res=-100
	else
		res=log10(res_tmp)
	end if
	!print*, "re 
	!end block
	if(ieee_is_nan(res))then
		print*, "nan detected get_ffuncs_series_fast"
		print*, "res_tmp=",res_tmp, dint, px, qx, t0, t1
		print*, "ipow=", ipow, ratio, ex, jum, line_slope, line_c
		print*, "logx0,logx1=",logx0, logx1, flag_end
		print*, "idx0,idx1=",idx0,idx1
		print*, "logxmin,logxmax=", logxmin,logxmax
		print*, "idxmin,idxmax=",idxmin,idxmax
		stop
	end if
	!read(*,*)
end subroutine

subroutine get_ffuncs_series(ex, exmax,exmin,mbhin,logr, ra,rp, spp, func, res, ipow)
	use, intrinsic :: ieee_arithmetic
	use com_sts_type
	use my_intgl
	use md_coeff
	use model_basic,only:dms,fc_share,fc_ir_share,ctl
	use MPI_comu,only:rid
	use md_star_pot
	implicit none
	external::func
	!type(s1d_type)::fphi
	type(star_pot_para)::spp
	real(8) ex,exmax,exmin,logr,fout,ratio,res,ra,rp
	real(8) ipow,yout,mbhin,tmin,tmax,fx_tmp,xmin,xmax
	integer debug
	integer idid,i 
	
	call get_phi_star_full_range(spp,logr, yout)

	select case(ctl%ebin_type)
	case(ebin_type_log)
		ratio=(10**yout+mbhin/10**logr)/(ex) 
	end select


	if(ratio<1)then
		print*, "get_ffuncs_series: error! ratio<1,ratio=", ratio
		select case(ctl%ebin_type)
		case(ebin_type_log)
			print*, "phi_star, r, ex=", 10**yout,10**logr, ex  
		end select
		stop
	end if

	!fout=0
	tmin=min(max((ratio-exmax/ex)/(ratio-1)*(1d0+1d-9),0d0),1d0)
	!if(tmin.ne.0.and.tmin.ne.1d0)then
	!	print*, "ratio, r, ex,exmax=",ratio,10**logr, ex, exmax
	!endif
	tmax=min((ratio-exmin/ex)/(ratio-1),1d0)
	!print*, "tmin=",tmin,tmax
	
	!call my_integral_none(tmin,tmax,fout,fcn_y,idid)
	!print*, "fout=",fout
	xmin=10**fc_ir_share%xmin
	xmax=10**fc_ir_share%xmax
	!print*, "ratio,logr=", ratio,logr, tmin, tmax,xmin,xmax
	if(tmax-tmin<1d-9)then
		call func((ratio-(ratio-1)*tmin)*ex,fx_tmp,xmin,xmax)
		fout=(tmax-tmin)*fx_tmp*tmin**ipow
		if(fout.le.0)then
			res=-100
		else
			res=log10((ratio-1)*fout)
			!print*, "ratio,ex=",ratio,ex
			!print*, "rp,ra,fx_tmp=",rp,ra,fx_tmp
			!print*, "tmin,tmax=",tmin,tmax
			!print*, "res=",res
			!return
		end if
		return
	end if

	fout=0

	call my_integral_acc(tmin,tmax,fout,1d-12,1d-9,fcn_y,idid)
	!call my_integral_none(0d0,1d0,fout,fcn_y,idid)
	!print*, "fout=",fout
	!read(*,*)



	if(idid<0)then
		print*, "get_ffuncs_series:ipow,logr=", ipow, logr
		print*, "ratio,ex=", ratio,ex,rid
		print*, "tmin,tmax=", tmin, tmax
		print*, "rp,ra=",rp,ra
		!do i=1, dms%n
		!	print*, "i=",i
		!	call fc_share%print("barge")
		call fc_ir_share%print("barge")
		!end do
		call spp%fphi_star%print("get_ffuncs_series: phi")
		stop
	end if
	if(fout<-0.1d0)then
		print*, "ex, logr, ratio, fout=", ex, logr, ratio, fout
		call spp%fphi_star%print("get_ffuncs_series: phi")
		read(*,*)
	end if
	!print*, "fout=",fout
	!read(*,*)
	if(ratio.eq.1)then
		print*, "get_ffuncs_series:ratio=", ratio
		print*, "logr,ex,phi=",logr, ex, yout
		print*, dms%mb(1)%dc%s2_dee%xcenter
		read(*,*)
	end if

	if(fout<=0)then 

		res=-100
		return
	end if
	!res=0d0
	!print*, "ie=", ieee_is_finite(log10(res)), ieee_is_finite(0d0)
	!stop
	res=log10(fout*(ratio-1))
	!if(tmax-tmin<1d-9)then
	!	print*, "int:res=",res
	!	stop
	!end if
	!print*, "logr, res, fout=", logr, res, fout
	if(ieee_is_nan(res).or.(.not.ieee_is_finite(res)))then
		print*, "error!"
		print*, "fout, ratio, ipow=", fout, ratio, ipow
		stop
	end if
	!read(*,*)
contains
	subroutine fcn_y(n, x, y, f, par, ipar)
		use, intrinsic :: ieee_arithmetic
		implicit none
		integer n, ipar(100)
		real(8) x, y(n), f(n), par(100)
		real(8) fx,s
		s=ratio-(ratio-1)*x
		call func(s*ex,fx,xmin,xmax)		
		f(1)=x**ipow*fx
		!if(debug.ge.1)then
		!	print*, "f(1),x,fx,x*ex=",f(1),x,fx,s, log10(s*ex)
		!end if
		if(.not.ieee_is_finite(f(1)).or.(abs((ratio-x))>0.1d0.and.ratio<x))then
			print*, "x,x*ex,fx, ratio=",x,x*ex, fx, ratio, logr, f(1),ratio-x,ratio<x
			!call dms%mb(1)%all%barge%print("all barge")
			stop
		end if
		
		
	end subroutine
end subroutine



subroutine get_ffuncs_series_rarp(ex,exmax,exmin, logr,  jum, jc, func, res, ipow)
	use, intrinsic :: ieee_arithmetic
	use com_sts_type
	use my_intgl
	use model_basic,only:fc_share,fc_ir_share
	implicit none
	external::func
	!type(sts_fc_type)::fphi
	real(8) ex,logr,fout,ratio,res, jum, jc
	real(8) ipow,yout,tmin,tmax,exmax,exmin,xmin,xmax,fx_tmp
	integer idid
	ratio=(jum*jc)**2/(10**(logr*2))/2d0/ex+1

	if(ratio<1)then
		print*, "get_ffuncs_series_rarp: error! ratio<1,ratio=", ratio
		stop
	end if
	fout=0

	tmin=min(max((ratio-exmax/ex)/(ratio-1)*(1d0+1d-9),0d0),1d0)
	!if(tmin.ne.0.and.tmin.ne.1d0)then
	!	print*, "tmin=",tmin,ratio, exmax,ex 
	!endif
	tmax=min((ratio-exmin/ex)/(ratio-1),1d0)
	!call my_integral_none(tmin,tmax,fout,fcn_y,idid)
	xmin=10**fc_ir_share%xmin
	xmax=10**fc_ir_share%xmax
	if(tmax-tmin<1d-9)then
		call func((ratio-(ratio-1)*tmin)*ex,fx_tmp,xmin,xmax)
		fout=(tmax-tmin)*fx_tmp*tmin**ipow
		if(fout.le.0)then
			res=-100
		else
			res=log10((ratio-1)*fout)
			!print*, "ratio,ex=",ratio,ex
			!print*, "rp,ra,fx_tmp=",rp,ra,fx_tmp
			!print*, "tmin,tmax=",tmin,tmax
			!print*, "res=",res
			!return
		end if
		return
	end if


	call my_integral_acc(tmin,tmax,fout,1d-12,1d-9,fcn_y,idid)
	!call my_integral_none(0d0,1d0,fout,fcn_y,idid)
	!call my_integral_none(1d0,ratio,fout,fcn_y,idid)
	if(idid<0)then
		print*, "get_ffuncs_series_rarp: ipow,logr=", ipow, logr
		print*, "ratio=", ratio
		print*, "tmin,tmax=",tmin,tmax
		call fc_ir_share%print("barge")
		stop
	end if
	if(fout<-0.1d0)then
		print*, "ex, ratio=", ex, logr, ratio, fout
		!call fphi%print("phi")
		read(*,*)
	end if
	!print*, "ex,jm,fout=",ex,jum,fout
	!read(*,*)
	if(ratio.eq.1)then
		print*, "ratio=", ratio
		print*, "logr,ex,phi=",logr, ex, yout
		read(*,*)
	end if

	if(fout<=0)then
		res=-100
		return
	end if
	!res=0d0
	!print*, "ie=", ieee_is_finite(log10(res)), ieee_is_finite(0d0)
	!stop
	!res=log10(fout/(ratio-1)**ipow)
	res=log10(fout*(ratio-1))
	if(ieee_is_nan(res).or.(.not.ieee_is_finite(res)))then
		print*, "error!"
		print*, "fout, ratio, ipow=", fout, ratio, ipow
		stop
	end if
contains
	subroutine fcn_y(n, x, y, f, par, ipar)
		use, intrinsic :: ieee_arithmetic
		implicit none
		integer n, ipar(100)
		real(8) x, y(n), f(n), par(100)
		real(8) fx,s
		s=ratio-(ratio-1)*x
		!call func(x*ex,fx)		
		call func(s*ex,fx,xmin,xmax)		
		f(1)=x**ipow*fx
		!f(1)=(abs(ratio-x))**ipow*fx

		!print*, "f(1),x,fx,x*ex=",f(1),x,fx, x*ex
		if(.not.ieee_is_finite(f(1)).or.(abs((ratio-x))>0.1d0.and.ratio<x))then
			print*, "x,x*ex,fx, ratio=",x,x*ex, fx, ratio, logr, f(1),ratio-x,ratio<x
			!call dms%mb(1)%all%barge%print("all barge")
			stop
		end if
		
	end subroutine
end subroutine
 

subroutine mb_get_dc_mpi_starpt_gx_iregular(mb,spp)
    use com_main_gw
    implicit none
    type(mass_bins),target::mb
    real(8),allocatable:: e_110(:,:),e_0(:,:),ee(:,:)
	real(8),allocatable:: j_111(:,:),j_rest(:,:),jj(:,:),ej(:,:)
    real(8),allocatable:: ycenter(:),xcenter(:)
    integer n,i,j
    type(coeff_type)::cej
    real(8) kappa, sigma32, n0,lambda
	type(star_pot_para)::spp
	!type(s1d_ird_type)::common_gx_ir
	interface 
		subroutine get_coeff_ej_spt_y_gx_iregular(xb, yb, n, spp,barge,&
				asymp, e_110,e_0,ee,j_111,j_reset,jj,ej)
			use com_main_gw
			implicit none
			!type(s1d_type)::phi_star
			type(star_pot_para)::spp
			integer n,i,j
			type(coeff_type)::ceo_ej
			real(8) xb(n),yb(n)
			real(8) e_110(n,n),e_0(n,n),ee(n,n) 
			real(8) j_111(n,n),j_reset(n,n),jj(n,n),ej(n,n)
			type(s1d_ird_type),target::barge
			!type(s2d_type)::gxj_ir
			real(8) asymp
		end subroutine
	end interface
	!print*, "start mb get dc"
    n=mb%df_coe_bins
	!print*, "n=", n
    allocate(ycenter(n),xcenter(n))
    allocate(e_110(n,n),e_0(n,n),ee(n,n), &
			j_111(n,n),j_rest(n,n),jj(n,n),ej(n,n))
    select case(mb%dc%jbin_type) 
    case(jbin_type_log)
        ycenter=10**mb%dc%s2_de_110%ycenter    
	case default
		print*, "error! define jbin_type_lin", mb%dc%jbin_type
		stop
    end select

	!call get_ffuncs(phi_star,10**mb%dc%s2_dj%xcenter,  &
	!	f0,f1,f2,n,mb%all%barge,mb%all%asymp)
	xcenter=10**mb%dc%s2_de_110%xcenter
	!print*, "xcenter=",xcenter
	call get_lambda(lambda)
    kappa=(4*pi*mb%mc)**2*lambda!*factor
	!print*, "kappa=",kappa
    !====================test=========================
    !kappa=(4*pi*mb%mc)**2*log(mb%mtot/mb%mc)!*factor
    !=================================================

    sigma32=(2*pi*mb%v0**2)**(-3/2d0)
    n0=mb%n0
	!print*, "ctl%dejmodel=",ctl%dejmodel
    select case(ctl%Dejmodel)
	case(dejmodel_EJ) 
		if(all(mb%all%barge_ir%fx(:).eq.0d0))then
			e_110=0
			e_0=0
			ee=0
			jj=0
			j_111=0
			j_rest=0
			ej=0
		else
			call get_coeff_ej_spt_y_gx_iregular(xcenter, ycenter,n,spp,&
			mb%all%barge_ir,mb%all%asymp, e_110,e_0,ee,j_111,j_rest,jj,ej)
		end if
		

		mb%dc%s2_de_110%fxy=e_110*sigma32*n0*kappa
		mb%dc%s2_de_0%fxy=e_0 *sigma32*n0*kappa
		mb%dc%s2_dee%fxy=ee*sigma32*n0*kappa
		mb%dc%s2_dj_111%fxy=j_111 *sigma32*n0*kappa
		mb%dc%s2_dj_rest%fxy=j_rest *sigma32*n0*kappa
		mb%dc%s2_djj%fxy=jj*sigma32*n0*kappa
		mb%dc%s2_dej%fxy=ej*sigma32*n0*kappa
	
	case(dejmodel_xj)
		do i=1, n
			do j=1, n 
				mb%dc%s2_de_110%fxy(i,j)=cej%e_110*sigma32*n0*kappa
				mb%dc%s2_de_0%fxy(i,j)=cej%e_0 *sigma32*n0*kappa
				if(cej%ee<0d0) cej%ee=abs(cej%ee)
				mb%dc%s2_dee%fxy(i,j)=cej%ee*sigma32*n0*kappa
				mb%dc%s2_dj_111%fxy(i,j)=cej%j_111 *sigma32*n0*kappa
				mb%dc%s2_dj_rest%fxy(i,j)=cej%j_rest *sigma32*n0*kappa
				if(cej%jj<0d0) cej%jj=abs(cej%jj)
				mb%dc%s2_djj%fxy(i,j)=cej%jj*sigma32*n0*kappa
				mb%dc%s2_dej%fxy(i,j)=cej%ej*sigma32*n0*kappa
			end do
		end do
	
	end select
end subroutine
