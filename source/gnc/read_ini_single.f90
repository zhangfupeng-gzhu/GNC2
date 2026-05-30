
subroutine readin_model_par(fmodel)
	use com_main_gw
	use md_star_pot
	 
	use md_mbh_evl_acc
	use model_config
	implicit none
	character*(*) fmodel
	character*(200) model, inacmodel, outacmodel, ecmodel, dejmodel, model_intej
    character*(200) time_unit, byctype, bytype,str_byctype(5),taskmode
    real(8) mass(2),mr,mbh_factor, real_tmp(3),dbl_parameters(4)
    integer i,ier,int_parameters(3)
	character*(200) strall, str(12), str_massbin_mode, &
		str_jbin_bd,str_fj_bd, str_plunge_cr,str_method_int, &
		 str_ebin_type,  str_dc_grid_type, &
		str_adb_est_method, str_barge_grid_type,str_barge_evl_method, str_density_est_method, &
		str_dc_method
	character*(400) root_path
	integer nnum	
	type(type_para) tp
 !   real(8) ts_snap_dt_per_snap
    !print*, "0"
	open(unit=999,file=trim(adjustl(fmodel)),status='old')
    ! call readpar_int_sp(ctl%ntasks, 999, "#","=",ier)
	call read_one_para(999,"model_default_path", tp)
	read(unit=tp%str,fmt="(A200)") ctl%default_para_file_dir
	call GETENV("GWPATH",root_path)
	if(rid.eq.0)then
		print*, "GWPATH=",trim(adjustl(root_path))
	end if
	open(unit=1999,file=trim(adjustl(root_path))//trim(adjustl(ctl%default_para_file_dir)),status='old')
	! print*, "read_paras"
	call read_paras(1999,pa_default)
	close(unit=1999)
	!!call pa_default%print("default")
	call read_paras(999,pa_usr_set)
	!call pa_usr_set%print("usr_set")
	call override_paras(pa_default,pa_usr_set,pa_now_used)
	!call pa_now_used%print("Now")
	!print*, "start readin_mass_bins"
	close(999)

	!stop

end subroutine
subroutine read_check_point(funit,check_point_name)
	implicit none
	character*(10) str_check_point
	character*(*) check_point_name
	integer funit
100	read(funit,fmt="(A10)") str_check_point
	!print*, "str=",str_check_point
	if(str_check_point(1:1).eq."#") then
		goto 100
	end if

	if(str_check_point(1:1).ne."-")then
		print*, "error in check point 1: expecting `-' at ", trim(adjustl(check_point_name)), &
		"but get ", str_check_point
		stop
	end if
end subroutine

subroutine readin_mass_bins(fl)
	use com_main_gw
	use md_stellar_evolution
	implicit none
	character*(200) strall, str(12), str_massbin_mode,str_sg_data_mode, str_model_mode
	integer nnum, funit, i,j, num_comp, ier
	logical have_comp(10)
	character*(*) fl
	real(8) mb(20),m1(20),m2(20),mstep, mass_imf, norm, max_weight_n
	real(8) wn_close, n0,particle_n,min_particle_n, max_particle_n
	real(8) mtot, slope_clone,slope_n, mtot_bin(20),ntot_bin(20)
	integer clone_n_min,clone_n_max,i_idx
	integer model_mode
	integer,parameter::model_mode_same=1, model_mode_indvd=2

	funit=1999
	open(unit=funit,file=fl,status='old') 
	ctl%num_mdehnen=0 
	ctl%num_mplummer=0 
	
	call READPAR_STR_AUTO_sp(strall,str,12,Nnum,funit,"#",","," ") 
	read(str(1),fmt=*) ctl%m_bins
	read(str(2),fmt=*) str_massbin_mode

	call readpar_str_sp(str_sg_data_mode,funit, "#", "=")
    select case(trim(adjustl(str_sg_data_mode)))
    case("given")  
        ctl%ini_sample_sg_mode=ini_sample_mode_given
    case("mobse")  
        ctl%ini_sample_sg_mode=ini_sample_mode_mobse
		call prepare_bse_brown_data() 
     end select
	
	call readpar_str_sp(str_model_mode,funit,"#", "=")
	select case(trim(adjustl(str_model_mode)))
	case("SAME")   ! all mass bins share the same model of density profile
		model_mode=model_mode_same
	case("INDVD")  ! assign models for individual mass bins
		model_mode=model_mode_indvd
	case default
		print*, "define model input mode:", trim(adjustl(str_model_mode))
		stop
	end select

	call read_check_point(funit,"start of mass bin")
	select case (trim(adjustl(str_massbin_mode)))
	case("KROUPA") !mass bin given by kroupa function
		
		call readpar_dbl_sp(ctl%metal_z,funit,"#","=",ier)
		do i=1, ctl%m_bins
			call skip_comments(funit,"#")
			!print*, "i=",i
			read(funit,fmt=*) i_idx,ctl%bin_mass_m1(i),ctl%bin_mass_m2(i),&
			ctl%ini_weight_n(i),ctl%clone_factor(i) 
			call read_check_point(funit, "seperating mass bins")
		end do 

		call canonical_IMF_prepare_n(imf_para_nt,imf_para_nc,imf_para_nq)
		call canonical_IMF_prepare_m(imf_para_mt,imf_para_mc,imf_para_mq)

		call get_kroupa_mnfrac(ctl%bin_mass_m1(1:ctl%m_bins),ctl%bin_mass_m2(1:ctl%m_bins), &
				mtot_bin(1:ctl%m_bins),ntot_bin(1:ctl%m_bins),ctl%m_bins)
		call set_single_model(funit,mtot_bin(1:ctl%m_bins),mtot)
		ctl%asymptot_ini(1,1:ctl%m_bins)=ntot_bin(1:ctl%m_bins)
		ctl%bin_mass(1:ctl%m_bins)=mtot_bin(1:ctl%m_bins)/ntot_bin(1:ctl%m_bins)*imf_para_mq/imf_para_nq
		ctl%bin_fracmass(1:ctl%m_bins)=mtot_bin(1:ctl%m_bins)
		
		ctl%ini_mass_bin_mode=ini_mass_bin_mode_kroupa
		
		! for Kroupa mas function, the first mass bin is for brown dwarfs!
		ctl%asymptot_ini(2,2:)=1
		ctl%asymptot_ini(6,1)=1
	case("TOPH") !mass bin given by modified kroupa function
		
		call readpar_dbl_sp(ctl%metal_z,funit,"#","=",ier)
		do i=1, ctl%m_bins
			call skip_comments(funit,"#")
			!print*, "i=",i
			read(funit,fmt=*) i_idx,ctl%bin_mass_m1(i),ctl%bin_mass_m2(i),&
			ctl%ini_weight_n(i),ctl%clone_factor(i)
			!print*, "m1,m2=", ctl%bin_mass_m1(i),ctl%bin_mass_m2(i)
			call read_check_point(funit, "seperating mass bins")
		end do

		call triple_break_IMF_prepare_n(imf_para_nt,imf_para_nc,imf_para_nq,topheavy_alpha,topheavy_xb)
		call triple_break_IMF_prepare_m(imf_para_mt,imf_para_mc,imf_para_mq,topheavy_alpha,topheavy_xb)

		call get_topheavy_mnfrac(ctl%bin_mass_m1(1:ctl%m_bins),ctl%bin_mass_m2(1:ctl%m_bins), &
				mtot_bin(1:ctl%m_bins),ntot_bin(1:ctl%m_bins),ctl%m_bins)
		call set_single_model(funit,mtot_bin(1:ctl%m_bins),mtot)
		ctl%asymptot_ini(1,1:ctl%m_bins)=ntot_bin(1:ctl%m_bins)
		ctl%bin_mass(1:ctl%m_bins)=mtot_bin(1:ctl%m_bins)/ntot_bin(1:ctl%m_bins)*imf_para_mq/imf_para_nq
		ctl%bin_fracmass(1:ctl%m_bins)=mtot_bin(1:ctl%m_bins)
		
		ctl%ini_mass_bin_mode=ini_mass_bin_mode_topheavy
		
		! for Kroupa mas function, the first mass bin is for brown dwarfs!
		ctl%asymptot_ini(2,2:)=1
		ctl%asymptot_ini(6,1)=1
	case("LOGBIN")  ! mass bin given by power law function
		call READPAR_STR_AUTO_sp(strall,str,12,Nnum,funit,"#",","," ")
		!print*, nnum
		read(str(1),fmt=*) ctl%bin_mass_min
		read(str(2),fmt=*) ctl%bin_mass_max
		read(str(3),fmt=*) mass_imf
		read(str(4),fmt=*) min_particle_n
		read(str(5),fmt=*) max_particle_n
		read(str(6),fmt=*) clone_n_min
		read(str(7),fmt=*) clone_n_max

		do i=1, ctl%m_bins
			ctl%asymptot_ini(2,i)=1d0
		end do

		call set_range(mb(1:ctl%m_bins), ctl%m_bins,log10(ctl%bin_mass_min),log10(ctl%bin_mass_max),0)
		ctl%bin_mass(1:ctl%m_bins)=10**mb(1:ctl%m_bins)
		mstep=mb(2)-mb(1)
		m1(1)=log10(ctl%bin_mass_min)
		m2(ctl%m_bins)=log10(ctl%bin_mass_max)
		!print*, "mstep=",mstep
		do i=2, ctl%m_bins
			m1(i)=m1(i-1)+mstep+1d-4
		end do
		do i=ctl%m_bins-1, 1, -1
			m2(i)=m2(i+1)-mstep
		end do
		!print*, "m1=",m1(1:ctl%m_bins)
		!print*, "m2=",m2(1:ctl%m_bins)
		ctl%bin_mass_m1(1:ctl%m_bins)=10**m1(1:ctl%m_bins)
		ctl%bin_mass_m2(1:ctl%m_bins)=10**m2(1:ctl%m_bins)
		
		norm=(ctl%bin_mass_max**(mass_imf+1)-ctl%bin_mass_min**(mass_imf+1))/(mass_imf+1d0)
		
		!write(*,fmt="(A4,20F10.4)") "m1=",ctl%bin_mass_m1(1:ctl%m_bins)
		!write(*,fmt="(A4,20F10.4)") "mc=", ctl%bin_mass(1:ctl%m_bins)
		!write(*,fmt="(A4,20F10.4)") "m2=",ctl%bin_mass_m2(1:ctl%m_bins)
		!ctl%asymptot_ini(1,:)=ctl%asymptot_ini(1,:)/sum(ctl%asymptot_ini(1,1:ctl%m_bins))
		!print*, sum(ctl%asymptot_ini(1,1:ctl%m_bins))

		mtot_bin(1:ctl%m_bins)=(10**(m2(1:ctl%m_bins)*(mass_imf+1))-10**(m1(1:ctl%m_bins)*(mass_imf+1)))/(mass_imf+1)/norm
		call set_single_model(funit,mtot_bin,mtot)
		
		slope_clone=(clone_n_max-clone_n_min)/(log10(ctl%bin_mass_max)-log10(ctl%bin_mass_min))
		slope_n=(max_particle_n-min_particle_n)/(log10(ctl%bin_mass_max)-log10(ctl%bin_mass_min))
		do i=1, ctl%m_bins
			particle_n=slope_n*log10(ctl%bin_mass(i)/ctl%bin_mass_min)+min_particle_n
			!print*, "particle_n=",particle_n
			max_weight_n=mtot_bin(i)*m0_cl/ctl%bin_mass(i)/particle_n
!			print*, "ini_weight_n:mbins=",max_weight_n
			!ctl%bin_mass_particle_number(ctl%m_bins)=min_particle_n
			call find_close_number(max_weight_n, wn_close)
			ctl%ini_weight_n(i)=wn_close
			ctl%clone_factor(i)=int(slope_clone*log10(ctl%bin_mass(i)/ctl%bin_mass_min))+clone_n_min
		end do
		do i=1, ctl%m_bins
			write(*, fmt="(A20, 5E13.5, I8)") "m1,mc,m2,aymp, iwn, weight_clone=", &
			ctl%bin_mass_m1(i),ctl%bin_mass(i), &
			ctl%bin_mass_m2(i), mtot_bin(i), ctl%ini_weight_n(i), ctl%clone_factor(i)
		end do
		ctl%ini_mass_bin_mode=ini_mass_bin_mode_pow
	case("GIVEN")
		call skip_comments(funit,"#")
		
		select case(model_mode)
			case(model_mode_indvd)
				do i=1, ctl%m_bins
					call skip_comments(funit,"#")
					read(funit,fmt=*) ctl%bin_mass_m1(i),ctl%bin_mass(i),ctl%bin_mass_m2(i),&
					ctl%ini_weight_n(i),ctl%clone_factor(i)
					read(funit,fmt=*) ctl%asymptot_ini(2:7,i)
				
					!read(funit,fmt=*) str_massbin_mode
					!print*, str_massbin_mode
					!stop
					call readin_models(funit,i)
				end do
				
				do i=2,ctl%m_bins
					if(ctl%bin_mass_m1(i).eq.ctl%bin_mass_m2(i-1))then
						print*, "error! bin_mass_m1=bin_mass_m2", i
						print*, "m1=",ctl%bin_mass_m1(i), ctl%bin_mass_m2(i-1)
						stop
					end if
				end do
			case(model_mode_same)
				do i=1, ctl%m_bins
					call skip_comments(funit,"#")
					read(funit,fmt=*) ctl%bin_mass_m1(i),ctl%bin_mass(i),ctl%bin_mass_m2(i),&
					ctl%ini_weight_n(i),ctl%clone_factor(i)
					read(funit,fmt=*) ctl%asymptot_ini(2:7,i)
				end do
				call skip_comments(funit,"#")
				read(funit,fmt=*) ctl%asymptot_ini(1,1:ctl%m_bins)
				mtot_bin(1:ctl%m_bins)=ctl%asymptot_ini(1,1:ctl%m_bins)*ctl%bin_mass(:ctl%m_bins)
				call set_single_model(funit,mtot_bin(1:ctl%m_bins),mtot)
		end select
		do i=1, ctl%m_bins
			if(ctl%bin_mass_m1(i).eq.ctl%bin_mass_m2(i))then
				ctl%bin_mass_m2(i)=ctl%bin_mass_m1(i)*1.0000001
			end if
		end do
		ctl%ini_mass_bin_mode=ini_mass_bin_mode_given
	case default
		print*, "readin_mass_bins:error!"
		stop
	end select
	!do i=1, ctl%m_bins
		!print*, "i,ml,ml_in=",i,ctl%ini_model_list(i),ctl%ini_model_list_in(i)
	!end do
	ctl%asymptot_now=ctl%asymptot_ini
	ctl%bin_fracmass_now=ctl%bin_fracmass
	call check_bins()
	call update_particle_existence()
	close(unit=funit)
	

	! call get_brown_evl()
end subroutine
subroutine check_bins()
	use model_basic
	implicit none
	integer i,j
	real(8) m1, m2, m1j,m2j
	logical intervals_overlap
	do i=1, ctl%m_bins
		if(ctl%bin_mass_m1(i)>ctl%bin_mass_m2(i))then
			print*, "bin i, m1>m2", i, ctl%bin_mass_m1(i),ctl%bin_mass_m2(i)
			stop
		end if
		m1=ctl%bin_mass_m1(i)
		m2=ctl%bin_mass_m2(i)
		do j=i+1, ctl%m_bins
			m1j=ctl%bin_mass_m1(j)
			m2j=ctl%bin_mass_m2(j)
			if(intervals_overlap(m1,m2,m1j,m2j))then
				print*, "error!, bin i and j overlap", i, j, m1,m2,m1j,m2j
				stop
			end if
		end do
	end do
end subroutine
logical function intervals_overlap(a, b, c, d)
    implicit none
    real(8), intent(in) :: a, b, c, d
	real(8)  a1,b1,c1,d1
	a1=min(a,b)
	b1=max(a,b)
	c1=min(c,d)
	d1=max(c,d)
    ! 检查重叠条件
    intervals_overlap=(a1<d1).and.(c1<b1)
end function intervals_overlap
subroutine get_kroupa_mnfrac(m1, m2, mfrac, nfrac, n)
	implicit none
	integer n,i
	real(8) m1(n),m2(n),mfrac(n), nfrac(n)
	real(8) t(4),c(3), q,canonical_IMF_func_mfrac_at_bin,canonical_IMF_func_nfrac_at_bin

	call canonical_IMF_prepare_m(t,c,q)
	!print*, "t=",t
	do i=1, n
		mfrac(i)=canonical_IMF_func_mfrac_at_bin(t,c,q, m1(i),m2(i))
		!print*, "i,m1,m2, mfrac(i)=",i,m1(i), m2(i), mfrac(i)
	end do

	call canonical_IMF_prepare_n(t,c,q)
	!!print*, "t=",t
	do i=1, n
		nfrac(i)=canonical_IMF_func_nfrac_at_bin(t,c,q, m1(i),m2(i))

	end do
end subroutine

subroutine get_topheavy_mnfrac(m1, m2, mfrac, nfrac, n)
	use model_basic,only:topheavy_alpha,topheavy_xb
	implicit none
	integer n,i
	real(8) m1(n),m2(n),mfrac(n), nfrac(n)
	real(8) t(4),c(3), q,triple_break_IMF_func_mfrac_at_bin
	real(8) triple_break_IMF_func_nfrac_at_bin

	call triple_break_IMF_prepare_m(t,c,q,topheavy_alpha,topheavy_xb)
	!print*, "t=",t
	do i=1, n
		mfrac(i)=triple_break_IMF_func_mfrac_at_bin(t,c,q, m1(i),m2(i),topheavy_alpha,topheavy_xb)
		!print*, "i,m1,m2, mfrac(i)=",i,m1(i), m2(i), mfrac(i)
	end do

	call triple_break_IMF_prepare_n(t,c,q,topheavy_alpha,topheavy_xb)
	!!print*, "t=",t
	do i=1, n
		nfrac(i)=triple_break_IMF_func_nfrac_at_bin(t,c,q, m1(i),m2(i),topheavy_alpha,topheavy_xb)

	end do
end subroutine


subroutine update_particle_existence()
	use com_main_gw
	implicit none
	integer i
	
	do i=1, n_tot_comp
		if(sum(ctl%asymptot_now(i+1,:))>0)then
			ctl%exist_stellar_type(i)=1
		else
			ctl%exist_stellar_type(i)=0
		end if
		select case(i)
		case(1)
			if(ctl%exist_stellar_type(i).ge.1)then
				ctl%idxstar=1
			else
				ctl%idxstar=-1
			end if
		case(2)
			if(ctl%exist_stellar_type(i).ge.1)then
				ctl%idxsbh=1
			else
				ctl%idxsbh=-1
			end if
		case(3)
			if(ctl%exist_stellar_type(i).ge.1)then
				ctl%idxbbh=1
			else
				ctl%idxbbh=-1
			end if
		case(4)
			if(ctl%exist_stellar_type(i).ge.1)then
				ctl%idxns=1
			else
				ctl%idxns=-1
			end if
		case(5)
			if(ctl%exist_stellar_type(i).ge.1)then
				ctl%idxwd=1
			else
				ctl%idxwd=-1
			end if
		case(6)
			if(ctl%exist_stellar_type(i).ge.1)then
				ctl%idxbd=1
			else
				ctl%idxbd=-1
			end if
		case(7)
			if(ctl%exist_stellar_type(i).ge.1)then
				ctl%idxrg=1
			else
				ctl%idxrg=-1
			end if
		case(8)
			if(ctl%exist_stellar_type(i).ge.1)then
				ctl%idxdm=1
			else
				ctl%idxdm=-1
			end if
		case(9)
			if(ctl%exist_stellar_type(i).ge.1)then
				ctl%idxnHe=1
			else
				ctl%idxnHe=-1
			end if
		end select
	end do

    ctl%idx_ref=1
    ctl%mass_ref=ctl%bin_mass(ctl%idx_ref)

end subroutine 
subroutine set_single_model(funit,mfracbin,mtot)
	use com_main_gw
	implicit none
	real(8) mfracbin(ctl%m_bins),mfrac_tot,mtot
	real(8) rho0
	integer funit, i
	
	call readin_models(funit,1)
	
	
	mfrac_tot=sum(mfracbin(1:ctl%m_bins))
	!print*, "mbins=",ctl%m_bins,ctl%num_mdehnen

		!mtot_bin(i)=(10**(m2(i)*(mass_imf+1))-10**(m1(i)*(mass_imf+1)))/(mass_imf+1)/norm*mtot
		!mtot_bin(i)=ctl%asymptot_ini(1,i)*ctl%bin_mass(i)/msun
		select case(trim(adjustl(ctl%str_ini_den_model(1))))
		case("Dehnen")			
			mtot=ctl%dehnen(1)%mtot
			ctl%num_mdehnen=0
			do i=1, ctl%m_bins
				ctl%str_ini_den_model(i)=ctl%str_ini_den_model(1)
				ctl%ini_model_list(i)=ini_den_model_dehnen
				ctl%num_mdehnen=ctl%num_mdehnen+1
				!print*, "ctl%num_mdehnen=",ctl%num_mdehnen
				ctl%ini_model_list_in(i)=ctl%num_mdehnen
				
				ctl%dehnen(i)%mtot=mtot*mfracbin(i)/mfrac_tot
				ctl%dehnen(i)%ra_crit=ctl%dehnen(1)%ra_crit
				ctl%dehnen(i)%gamma=ctl%dehnen(1)%gamma
			end do
		case("Plummer")
			mtot=ctl%plummer(1)%mtot
			ctl%num_mplummer=0
			do i=1, ctl%m_bins
				ctl%str_ini_den_model(i)=ctl%str_ini_den_model(1)
				ctl%ini_model_list(i)=ini_den_model_plummer
				ctl%num_mplummer=ctl%num_mplummer+1
				ctl%ini_model_list_in(i)=ctl%num_mplummer
				
				ctl%plummer(i)%mtot=mtot*mfracbin(i)/mfrac_tot
				ctl%plummer(i)%ra_crit=ctl%plummer(1)%ra_crit
			end do
		  
		end select
	

end subroutine
subroutine readin_models(funit,i)
	use model_basic
	use md_star_pot,only:spp_new
	use mpi_comu,only:rid
	implicit none
	integer funit
	integer i,ier
	character*(200) str_ini_model
	character*(200) str_units
	!integer,parameter::unit_numeric=1, unit_msunpc=2
	!integer unit_select
	real(8) mconv
	real(8) rconv
	call readpar_str_sp(str_ini_model,funit,"#","=")
	ctl%str_ini_den_model(i)=trim(adjustl(str_ini_model))
	
	call readpar_str_sp(str_units,funit,"#","=")
	select case(adjustl(trim(str_units)))
	case("Numeric")
		!unit_select=unit_numeric
		mconv=1d0
		rconv=1d0
	case("MsunPc")
		!unit_select=unit_msunpc
		mconv=m0_cl
		rconv=r0_cl/pc
	case default
		print*, "error! unit not defined", adjustl(trim(str_units))
		stop
	end select
	!print*, "mconv,rconv=",mconv,rconv
	select case(trim(adjustl(ctl%str_ini_den_model(i))))
	case("Dehnen")
		ctl%ini_model_list(i)=ini_den_model_dehnen
		ctl%num_mdehnen=ctl%num_mdehnen+1
		ctl%ini_model_list_in(i)=ctl%num_mdehnen
		!print*, "i, Dehnen, iin=", i, ctl%num_mdehnen
		call readpar_dbl_sp(ctl%dehnen(ctl%num_mdehnen)%mtot,funit, "#", "=", ier)
		call readpar_dbl_sp(ctl%dehnen(ctl%num_mdehnen)%ra_crit,funit, "#", "=", ier)
		call readpar_dbl_sp(ctl%dehnen(ctl%num_mdehnen)%gamma,funit, "#", "=", ier)
		!ctl%mtot=ctl%mtot+ctl%denhen(ctl%num_mdehnen)%mtot*m0_cl

		if(ctl%dehnen(ctl%num_mdehnen)%gamma<0.5.and.spp_new%mbh.ne.0d0)then
			print*, "warnning!!, for denhnenmodel with MBH mass!=0 and gamma<0.5, self-consistent solution is not garentee!!",spp_new%mbh
		end if
		
		ctl%dehnen(ctl%num_mdehnen)%mtot=ctl%dehnen(ctl%num_mdehnen)%mtot/mconv
		ctl%dehnen(ctl%num_mdehnen)%ra_crit=ctl%dehnen(ctl%num_mdehnen)%ra_crit/rconv
		!print*, "mtot, ra_crit=",ctl%dehnen(ctl%num_mdehnen)%mtot,ctl%dehnen(ctl%num_mdehnen)%ra_crit
		!stop
	case("Plummer")
		ctl%ini_model_list(i)=ini_den_model_plummer
		ctl%num_mplummer=ctl%num_mplummer+1
		ctl%ini_model_list_in(i)=ctl%num_mplummer
		call readpar_dbl_sp(ctl%plummer(ctl%num_mplummer)%mtot,funit, "#", "=", ier)
		call readpar_dbl_sp(ctl%plummer(ctl%num_mplummer)%ra_crit,funit, "#", "=", ier)
		
		if(spp_new%mbh_dmless/ctl%plummer(ctl%num_mplummer)%mtot>0.01d0)then
			print*, "warnning!!, for plummer model with MBH mass/Mtot>0.01, self-consistent solution is not garentee!!"
		end if
		
		ctl%plummer(ctl%num_mplummer)%mtot=ctl%plummer(ctl%num_mplummer)%mtot/mconv
		ctl%plummer(ctl%num_mplummer)%ra_crit=ctl%plummer(ctl%num_mplummer)%ra_crit/rconv
 
	case default
		print*, "error! define ini den model", trim(adjustl(ctl%str_ini_den_model(i)))
		stop
	end select

	

	call read_check_point(funit,"end of this mass bin")
end subroutine
subroutine find_close_number(nin, nout)
	implicit none
	real(8) nin, nout
	integer dg
	if(nin<=0)then
		print*, "nin should be larger than zero! stoped!"
		stop
	end if
	if(log10(nin)>0)then
		dg=int(log10(nin))+1
	else
		dg=int(log10(nin))
	end if
	if(nin>5*10**(dg-1)) then
		nout=5*10d0**(dg-1)
	else
		nout=10d0**(dg-1)
	end if
	print*, "nin, nout, dg=",nin, nout, dg
end subroutine 

subroutine print_current_code_version()
	use model_basic
	use MPI_comu
	implicit none
	if(rid.eq.0)then
		print*, "current code version: 2.2-beta : updated in 2026-05-28"
	end if
	
end subroutine
   