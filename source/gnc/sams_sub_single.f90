!subroutine sams_get_weight_single(sps,e0)
!	use com_main_gw
!	implicit none
!	type(particle_samples_arr_type)::sps
!	real(8) e0
!	integer i
!	do i=1, sps%n
!		call particle_sample_get_weight(sps%sp(i), e0)
!	end do
!	!print*, "3"
!end subroutine
 
subroutine update_jm(dm,spp,ex,jph_dmless, jm, jc_xy)
    use com_main_gw
    implicit none
    type(diffuse_mspec)::dm
    type(star_pot_para)::spp
    real(8) ex,logex,jm
    real(8) r_c_iter,jc_dmless, rmax, jc_xy
    real(8) rc, jph_dmless
    integer ier

    !logex=log10(ex)
    !call get_rmax_accurate(spp ,  dm%fr_phi, logex,rmax)

    !if(ex>1)then
    !    call spp%fphi_star%print("fphi_star")
    !end if

    rc=r_c_iter(spp, ex,ier)

    jc_xy=jc_dmless(rc,spp)
    !===========================
    !if(ctl%debug.ge.1)then
    !    print*, "i,sp%jm,jcf,jci=", sp%jm, jc_xy, sp%jc/(ctl%v0*r0_cl)
    !end if
    
    jm=jph_dmless/jc_xy
    !if(ctl%debug.ge.1)then
    !    print*, "i1,sp%jm=", sp%jm
    !end if
    call set_jm_bound(jm)
    
end subroutine

subroutine get_sample_para(dm,bks,replace,spp)
    use com_main_gw
    implicit none
    type(particle_samples_arr_type)::bks
    type(diffuse_mspec)::dm
    type(star_pot_para)::spp
    integer i
    real(8) t1,t2,xmin,xmax
    logical::replace

    if(rid.eq.0)then
        call cpu_time(t1)
    end if
    !energy_range_covered=.true.
    !xmin=bks%sp(1)%x
    !xmax=xmin
    do i=1, bks%n
        !if(xmin>bks%sp(i)%x) xmin=bks%sp(i)%x
        !if(xmax<bks%sp(i)%x) xmax=bks%sp(i)%x
        if(bks%sp(i)%x>emax_factor )then
            if(spp%mbh_dmless.eq.0.and. replace)then
                bks%sp(i)%x=emax_factor*0.999999
                print*, "get_sample_para:bks%sp(i)%x,emax_dstr_factor=",i, bks%sp(i)%x,emax_factor
            else
                cycle
            end if
           ! energy_range_covered=.false.
        end if
        if(bks%sp(i)%x<emin_factor.and.spp%mbh_dmless.eq.0)then
            print*, "get_sample_para:bks%sp(i)%x,emin_dstr_factor=",i, bks%sp(i)%x,emin_factor
            cycle
        end if
        call get_sample_para_one(dm,bks%sp(i),spp)
    end do
    if(rid.eq.0)then
        call cpu_time(t2)
    end if
    if(rid.eq.0)then
        print*, "get_sample_para: used time=",t2-t1, " s"!, " xmin, xmax=", xmin, xmax
    end if
end subroutine

subroutine get_sample_para_no_pd(dm,bks,replace,spp)
    use com_main_gw
    implicit none
    type(particle_samples_arr_type)::bks
    type(diffuse_mspec)::dm
    type(star_pot_para)::spp
    integer i
    real(8) t1,t2,xmin,xmax
    logical::replace

    if(rid.eq.0)then
        call cpu_time(t1)
    end if
    !energy_range_covered=.true.
    !xmin=bks%sp(1)%x
    !xmax=xmin
    !print*, "get_sample_para_no_pd",rid
    self_correction_emax=0
    do i=1, bks%n
        !if(xmin>bks%sp(i)%x) xmin=bks%sp(i)%x
        !if(xmax<bks%sp(i)%x) xmax=bks%sp(i)%x
        if(bks%sp(i)%x>emax_factor)then
            if(replace.and.spp%mbh_dmless.eq.0)then
                bks%sp(i)%x=emax_factor*0.999999
            end if
            self_correction_emax=self_correction_emax+1
            print*, "get_sample_para_no_pd:bks%sp(i)%x,emax_dstr_factor=", i, bks%sp(i)%x,emax_factor
            cycle
           ! energy_range_covered=.false.
        end if
        if(bks%sp(i)%x<emin_factor.and.spp%mbh_dmless.eq.0)then
            print*, "get_sample_para_no_pd:bks%sp(i)%x,emin_dstr_factor=", i, bks%sp(i)%x,emax_factor
            cycle
        end if
        call get_sample_para_one_no_pd(dm,bks%sp(i),spp)
    end do
    call collection_int(self_correction_emax)
    if(rid.eq.0)then
        call cpu_time(t2)
    end if
    if(rid.eq.0)then
        print*, "get_sample_para_no_pd: used time=",t2-t1, " s"!, " xmin, xmax=", xmin, xmax
    end if
end subroutine
subroutine get_sample_para_one_xj_no_reset(ex,jm,spp,fr_phi,  &
     jc_xy, rp_xy,ra_xy,pd_xy)
! use com_main_gw
 use com_sts_type
 use md_star_pot
 implicit none
 type(star_pot_para)::spp
 real(8) ex, logex, jc, jc_dmless
 real(8) rmax, rc,jm,jc_xy,rp_xy,ra_xy
 real(8) pd_xy,p_EJ_dmless_fast,r_c_iter
 type(s1d_type)::fr_phi
 integer ier

 logex=log10(ex)
 
 call get_rmax_accurate(spp,  fr_phi, logex,rmax)
 
 rc=r_c_iter(spp,ex,ier)

 jc_xy=jc_dmless(rc,spp)
 !jm0=jph_xy/jc_xy
 call get_rpra_dmless(spp, ex, jm, jc_xy, &
             log10(rc), rmax, rp_xy,ra_xy)
 pd_xy=p_EJ_dmless_fast(spp, ex,jm,  jc_xy, rp_xy,ra_xy)   
 !jm_xy=jm0
end subroutine
 
 subroutine get_sample_para_one_xj_rpra(ex,jm,jph_xy,spp,fr_phi,  &
    jm_xy, jc_xy, rp_xy,ra_xy)
    ! use com_main_gw
    use com_sts_type
    use md_star_pot
    implicit none
    type(star_pot_para)::spp
    real(8) ex, logex, jc, jc_dmless,jph_xy,jm_xy
    real(8) rmax, rc,jm,jc_xy,rp_xy,ra_xy,jph
    real(8) pd_xy,p_EJ_dmless,jm0, r_c_iter
    type(s1d_type)::fphi_star,fr_phi,fma_star,frho_star
    integer ier

    logex=log10(ex)
    
    call get_rmax_accurate(spp,  fr_phi, logex,rmax)
 
    rc=r_c_iter(spp,ex,ier)
    
    jc_xy=jc_dmless(rc,spp)
    jm0=jph_xy/jc_xy
    call get_rpra_dmless(spp, ex, jm0, jc_xy, &
                log10(rc), rmax, rp_xy,ra_xy)
    jm_xy=jm0
end subroutine

subroutine get_kpl_ae(jm, jc, mbh, rp, ac,ec)
    implicit none
    real(8) jm, jc, mbh, rp, ac, ec
    ec=(jm*jc)**2/mbh/rp-1
    ac=rp/(1-ec)
end subroutine
!subroutine sams_get_weight_real_one(sp)
!    use com_main_gw
!    implicit none
!    type(particle_sample_type)::sp
!    call sams_get_weight_clone_single_one(sp)
!    call get_sample_weight_real(sp)
!end subroutine
subroutine sams_get_weight_clone_single_one(sp)
    use com_main_gw
	implicit none
	type(particle_sample_type)::sp
	integer i, amplifier, obidx
    real(8) en
	integer mass_idx,nlvl

    obidx=sp%obidx
    call get_mass_idx(sp%m, mass_idx)
    amplifier=ctl%clone_factor(mass_idx)
	  
    if(sp%exit_flag.eq.exit_boundary_max)then
        !print*, "sp%x,bx=",sp%x, sp%en/ctl%energy0, sp%byot_bf%e/ctl%energy0,sp%jm, sp%rp, sp%r_lc, &
        !    sp%byot%e_bin,sp%byot_bf%l
        !read(*,*)
        en=sp%byot_bf%e
        
    else
        en=sp%en
    end if
        
    call particle_sample_get_weight_clone(en, ctl%clone_scheme, &
        amplifier,ctl%clone_e0,sp%weight_clone,nlvl)
    sp%Lvl_clone=nlvl
    if(ieee_is_nan(sp%weight_clone))then
        print*, "sams_get_weight_clone_single:NaN", sp%weight_clone, &
        sp%source, sp%obtype, sp%obidx, sp%en, &
        sp%id
        stop
    end if

end subroutine
subroutine sams_get_weight_clone_single(sps)
	use com_main_gw
	implicit none
	type(particle_samples_arr_type)::sps
	integer i, amplifier, obidx
    real(8) en
	integer mass_idx
	!print*, "sps%n=",sps%n
	do i=1, sps%n
        call sams_get_weight_clone_single_one(sps%sp(i))
	end do
end subroutine
!subroutine get_wsi(spin, inc_in, wsi)
!	use model_basic
!	implicit none
!	real(8) spin, inc_in,inc, wsi
!	integer isel,i
!	if(inc_in>pi) then
!		inc=2*pi-inc_in
!	else
!		inc=inc_in
!	end if
!
!	if(inc<spin_data_0999(1)) then
!		isel=1
!		goto 100
!	end if
!	if(inc>spin_data_0999(11)) then
!		isel=11
!		goto 100
!	end if
!
!loop1: do i=1, 10
!		if(inc>=spin_data_0999(i).and.inc<=spin_data_0999(i+1))then
!			isel=i
!			exit loop1
!		end if
!	   end do loop1
!100	wsi=wsi_data_0999(isel)
!	!print*, "inc_in, inc, isel, wsi=",inc_in,inc,isel, wsi
!	!read(*,*)
!	return
!end subroutine
 
subroutine set_clone_weight(sms)
	use com_main_gw
	implicit none
	type(chain_type)::sms
	type(chain_pointer_type),pointer::pt
	!integer i, obidx, amplifier
    real(8) en
	!integer mass_idx

    pt=>sms%head
    do while(associated(pt).and.allocated(pt%ob))
        !obidx=pt%ob%obidx
        select type(ca=>pt%ob)
        type is(particle_sample_type)
            call sams_get_weight_clone_single_one(ca)
         
        end select
        pt=>pt%next
    end do		
    
end subroutine
 
  

subroutine sams_arr_select_type_single(sps, sps_out, obtype)
	use com_main_gw
	implicit none
	type(particle_samples_arr_type)::sps,sps_out
	integer nsel,i, obtype
	nsel=0
	do i=1, sps%n
		if(selection()) nsel=nsel+1
	end do
	!print*, "nsel=",nsel
	call sps_out%init(nsel)
	nsel=0
	do i=1, sps%n
		if(selection()) then
			nsel=nsel+1
			sps_out%sp(nsel)=sps%sp(i)
		end if		
	end do
contains
	function selection()
		implicit none
		logical selection
		selection=.false.
		
		if(sps%sp(i)%obtype.eq.obtype)then
				selection=.true.
		end if
	end function
end subroutine	 

subroutine get_obidx_from_type_sg(stellar_type,idx)
    use com_main_gw
    implicit none
    integer,intent(in):: stellar_type
    integer,intent(out)::idx
    select case(stellar_type)
    case(star_type_ms)
        idx=1
    case(star_type_bh)
        idx=2
    case(star_type_ns)
        idx=3
    case(star_type_wd)
        idx=4
    case(star_type_bd)
        idx=5
    case(star_type_rg)
        idx=6
    case(star_type_dark_matter)
        idx=7
    case(star_type_nakedHe)
        idx=8
    case default
        idx=-1
        !print*, "get_obidx_from_type_sg: error! define stellar type", &
        !    stellar_type
    end select
    !print*, "get_obidx_from_type_sg", stellar_type,idx
end subroutine 

subroutine get_type_from_ctl_obidx_sg(idx,stellar_type)
    use com_main_gw
    implicit none
    integer idx, stellar_type
    select case(idx)
    case(1)
        stellar_type=star_type_ms
    case(2)
        stellar_type=star_type_bh
    case(3)        
        stellar_type=star_type_ns
    case(4)
        stellar_type=star_type_wd
    case(5)
        stellar_type=star_type_bd
    case(6)
        stellar_type=star_type_rg
    case(7)
        stellar_type=star_type_dark_matter
    case(8)
        stellar_type=star_type_nakedHe
    !case(7)
     !   stellar_type=bytype_bhb
    !case(8)
    !    stellar_type=bytype_msb
    case default
        print*, "get_type_from_ctl_obidx_sg: error! define idx=", &
            idx
    end select
end subroutine
  
subroutine get_mass_idx(m,idx)
    use com_main_gw
    implicit none
    integer i, idx
    real(8) m
    idx=-1
    do i=1, ctl%m_bins
        if(m.ge.ctl%bin_mass_m1(i).and.m.lt.ctl%bin_mass_m2(i))then
            idx=i
			!print*, m, idx
            return
        end if
    end do
    ! if(m.ge.ctl%bin_mass_m1(ctl%m_bins).and.m.le.ctl%bin_mass_m2(ctl%m_bins))then
    !     idx=ctl%m_bins
    !     return
    ! end if
    print*, "get_mass_idx idx=-1: m=", m
    if(rid.eq.0)then
        do i=1, ctl%m_bins
            print*, "i,m1,m2=",i,ctl%bin_mass_m1(i),ctl%bin_mass_m2(i)
        end do
    end if
    !stop
    
end subroutine 
  

subroutine get_star_type(type_int, str)
    use com_main_gw
    implicit none
    integer type_int
    character*(*) str
    select case(type_int)
    case (star_type_BH)
        str="BH"
    case (star_type_MS)
        str="MS"
    case (star_type_NS)
        str="NS"
    case (star_type_WD)
        str="WD"
    case (star_type_BD)
        str="BD"		
    case (star_type_rg)
        str="RG"		
    case (star_type_nakedHe)
        str="nHe"		
    case default
        str='UNKNOWN'
      !  print*, "unknown type type_int=",type_int
     !   stop
     !   print*, "define star_type:",type_int
     !   stop        
    end select
end subroutine

subroutine deallocate_chains_arrs()
    use com_main_gw
    implicit none
    integer i,j
    !print*, "1"
    call bksams%destory()
    !print*, "2"
    call bysams%destory()
    !print*, "3"
    if(allocated(bksams_arr%sp))deallocate(bksams_arr%sp) 
    if(allocated(bksams_arr_norm%sp))deallocate(bksams_arr_norm%sp)
    if(allocated(bksams_pointer_arr%pt))deallocate(bksams_pointer_arr%pt)
    
    do i=1, dms%n
        do j=1, n_tot_comp
            if(allocated(dms%mb(i)%dsp(j)%p%nejw)) deallocate(dms%mb(i)%dsp(j)%p%nejw)
            call dms%mb(i)%dsp(j)%p%deallocation()
			Nullify(dms%mb(i)%dsp(j)%p)
		end do	
    end do
    deallocate(dms%mb)
   
    if(allocated(df)) deallocate(df)
    !print*, "6"  
end subroutine

subroutine set_star_radius(pr)
    use com_main_gw
     
    implicit none
    type(particle)::pr
    real(8) star_Radius,white_dwarf_radius
    real(8) rnd, xtmp
    select case(pr%obtype)
    case(star_type_MS)
        pr%radius=star_Radius(pr%m)
	case(star_type_BD)
		pr%radius=star_Radius(pr%m) !pr%m/my_unit_vel_c**2
    case(star_type_BH)
        pr%radius=pr%m/my_unit_vel_c**2
    case(star_type_WD)
        pr%radius=white_dwarf_radius(pr%m)
    case(star_type_NS)
        pr%radius=10d3/AU_SI 
    case(star_type_rg)
        print*, "error! red giant radius is not set! finished the code here" 
        stop
    case(star_type_nakedHe)
        print*, "error! naked He radius is not set! finished the code here"
        stop
    case(star_type_dark_matter)
        pr%radius=0d0
    case default
        !print*, "error, "
		print*, "set_star_radius*****************star type=",pr%obtype
        pr%radius=0d0
		stop
    end select
end subroutine

subroutine get_riso_compact_obj_given_spin_inc(spin,inc,riso)
    use com_main_gw
    implicit none
    real(8) inc, riso, wsi
    integer spin

    select case(spin)
    case(1)
        !ac=-mbh/2d0/sp%en
        call iso_kerr%rpnw_incgr%get_value_d(inc, wsi)
        riso=wsi*mbh_radius
        !print*, "ac, sp%byot%inc, wsi=",ac, sp%byot%inc, wsi
       ! print*, "sp%r_lc=",sp%r_lc
        !read(*,*)
    case(0)
        riso=8*mbh_radius
    case default
        print*, "ini:error! define mbh spin", spin
        stop
    end select
end subroutine
 