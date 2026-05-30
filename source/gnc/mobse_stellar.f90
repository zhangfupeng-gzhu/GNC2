module md_mobse_stellar_single
    use com_sts_type
    implicit none
    type stellar_object_arr
		integer,allocatable:: ktype(:)
		real(8),allocatable:: mass(:), age(:)
        real(8),allocatable:: w(:), radius(:)
		integer n
	contains
		procedure::init=>init_stellar_object
        procedure::write_info=>write_info_stellar_object_arr
        procedure::read_info=>read_info_stellar_object_arr
	end type
	
    type,extends(stellar_object_arr)::stellar_object_hs
    integer,allocatable:: kwtype(:)
    contains
        procedure::init=>init_stellar_object_hs
        procedure::print=>print_stellar_history
        procedure::write_info=>write_info_stellar_object_hs
        procedure::read_info=>read_info_stellar_object_hs
    end type

    type stellar_samples_sts
        type(sts_fc_type):: mstar, mbbh, msbh, mbstar,mbd, msrg, mfgb
		type(sts_fc_type):: mwd, mns, rsrg, rfgb, rstar, rbd
        type(sts_fc_type):: bbhac, bbhec, bbhq, bbhm1, bbhm2
    end type

	integer,parameter::mobse_ktype_fcms=0
	integer,parameter::mobse_ktype_ms=1
	integer,parameter::mobse_ktype_bh=14
	integer,parameter::mobse_ktype_ns=13
    integer,parameter::mobse_ktype_srg=5              ! supper red giant
    integer,parameter::mobse_ktype_hgap=2
    integer,parameter::mobse_ktype_hrbr=4
    integer,parameter::mobse_ktype_sabr=6
    integer,parameter::mobse_ktype_fgb=3              ! first gaint branch
	integer,parameter::mobse_ktype_hewd=10
	integer,parameter::mobse_ktype_cowd=11
	integer,parameter::mobse_ktype_onwd=12
    real(8)::mobse_orgsample_IMF_slope=0d0 
	integer,parameter::model_imf_pow=1, model_imf_can=2, model_imf_top_heavy=3
	integer,parameter::gsp_model_method_weight=1,gsp_model_method_org=2
	real(8),parameter::mass_brown_dwarf_max=0.1d0
!===========================================
    integer num_select_sg_for_sts, num_select_for_type_mbin
    integer num_select_by_for_sts
    integer num_sg_simu, num_by_simu, mobse_include_by
    real(8) simu_time_tot
    real(8) alpha_IMF, by_frac, metal_z, mobse_mass_min
	real(8) mobse_mass_max
	integer model_imf, mobse_idx_ref, gsp_model_method
!===========================================
	private::init_stellar_object_hs
	private::init_stellar_object
    private::write_info_stellar_object_arr
    private::read_info_stellar_object_hs
	private::write_info_stellar_object_hs
	private::print_stellar_history
contains

    subroutine init_stellar_object_hs(this, n)
        implicit none
        class(stellar_object_hs)::this
        integer n
        call this%stellar_object_arr%init(n)
        if(allocated(this%kwtype))then
            deallocate(this%kwtype)
        end if
        allocate(this%kwtype(n))
    end subroutine
    
    subroutine write_info_stellar_object_arr(this, funit)
        implicit none
        class(stellar_object_arr)::this
        integer funit
        write(unit=funit) this%n
		if(this%n<1) return
        write(unit=funit) this%ktype(1:this%n), &
            this%mass(1:this%n), this%age(1:this%n), &
            this%w(1:this%n), this%radius(1:this%n)
    end subroutine
    subroutine read_info_stellar_object_arr(this, funit)
        implicit none
        class(stellar_object_arr)::this
        integer funit
        read(unit=funit) this%n
        call this%init(this%n)
		if(this%n<1) return
        read(unit=funit) this%ktype(1:this%n), &
        this%mass(1:this%n), this%age(1:this%n) ,&
        this%w(1:this%n), this%radius(1:this%n)
    end subroutine
    subroutine write_info_stellar_object_hs(this, funit)
        implicit none
        class(stellar_object_hs)::this
        integer funit
        call write_info_stellar_object_arr(this,funit)
		if(this%n<1)return
        write(funit) this%kwtype(1:this%n)
    end subroutine
    subroutine read_info_stellar_object_hs(this, funit)
        implicit none
        class(stellar_object_hs)::this
        integer funit
        call read_info_stellar_object_arr(this, funit)
		if(this%n<1)return
        if(.not.allocated(this%kwtype))allocate(this%kwtype(1:this%n))
        read(unit=funit) this%kwtype(1:this%n)
    end subroutine
    
    subroutine init_stellar_object(this,n)
        implicit none
        class(stellar_object_arr)::this
        integer n,ktype
        if(allocated(this%mass))then
            deallocate(this%mass,this%age,this%ktype, this%w,this%radius)
        end if
        allocate(this%mass(n),this%age(n),this%ktype(n),this%w(n),this%radius(n))
        this%n=n
    end subroutine
    
    subroutine print_stellar_history(this)
        class(stellar_object_hs)::this
        integer i
        character*(9) kwtype,kstartype
        write(*,fmt="(8A15)") "time", "mass","radius", "type", "event"
        do i=1, this%n
            kstartype=get_kstar_type(this%ktype(i))
            kwtype=get_kw_type(this%kwtype(i) )
            write(*,fmt="(3F15.3, A15, A15)") this%age(i), this%mass(i), this%radius(i), &
                 kstartype,kwtype 
        end do
    end subroutine    

    subroutine get_snapshot_of_samples(sh,n, curtime, star, bh, ns, wd, bd, rg, fgb)
        implicit none
        integer nstar, nbh, i, j, n, nwd, nns,nbd, nrg, nfgb
        type(stellar_object_hs)::sh(n)
        type(stellar_object_arr)::star, bh, wd, ns, bd, rg, fgb
        real(8) curtime	
        real(8),parameter::tiny=1d-5
        integer,allocatable::idx(:)
        allocate(idx(n))
         
        nstar=0; nbh=0; nwd=0; nns=0;nbd=0; nrg=0; nfgb=0
        if(n.eq.0)then 
            print*, "warnning: stellar_history_object_n=0"
            return
        end if
        do i=1, n
            if(curtime<sh(i)%age(1).or.curtime>sh(i)%age(sh(i)%n)+tiny)then
                print*, "i=",i,"curtime not in rage", curtime, sh(i)%age(1),&
                 sh(i)%age(sh(i)%n)
                cycle
            end if
            
    loop1:		do j=1, sh(i)%n-1
                if((sh(i)%age(j).le.curtime).and.(sh(i)%age(j+1)>curtime))then
                    idx(i)=j
                    exit loop1
                end if
            end do loop1
            !call sh(i)%print()
            !print*, "ctime,idx=",curtime,idx(i)
            if(abs(sh(i)%age(sh(i)%n)-curtime)<tiny)then
                idx(i)=sh(i)%n
                !goto 100
            end if
            select case(sh(i)%ktype(idx(i)))
            case(mobse_ktype_bh)
                nbh=nbh+1
			case(mobse_ktype_fcms)
				if(sh(i)%mass(idx(i))<mass_brown_dwarf_max)then
					nbd=nbd+1
				else
					nstar=nstar+1
				end if
            case(mobse_ktype_ms)
                nstar=nstar+1
            case(mobse_ktype_cowd,mobse_ktype_hewd, mobse_ktype_onwd)
                nwd=nwd+1
            case(mobse_ktype_ns)
                nns=nns+1
            case(mobse_ktype_srg,mobse_ktype_hgap,mobse_ktype_sabr,mobse_ktype_hrbr)
                nrg=nrg+1
            case(mobse_ktype_fgb)
                nfgb=nfgb+1
            case default
                !print*, curtime, SH(i)%mass(idx(i)),&
                !	  get_kstar_type(sh(i)%ktype(idx(i)))
                !read(*,*)
            end select
        end do
        call bd%init(nbd)
        bd%ktype=mobse_ktype_fcms
        call star%init(nstar)
        star%ktype=mobse_ktype_ms
        call bh%init(nbh)
        bh%ktype=mobse_ktype_bh
        call wd%init(nwd)
        wd%ktype=mobse_ktype_cowd
        call ns%init(nns)
        ns%ktype=mobse_ktype_ns
        call rg%init(nrg)
        rg%ktype=mobse_ktype_srg
        call fgb%init(nfgb)
        fgb%ktype=mobse_ktype_fgb
        
        nstar=0; nbh=0;nwd=0; nns=0; nbd=0; nfgb=0
        do i=1, n
            select case (sh(i)%ktype(idx(i)))
            case (mobse_ktype_bh)
                nbh=nbh+1
                bh%mass(nbh)=sh(i)%mass(idx(i))
                bh%age(nbh)=sh(i)%age(idx(i))
			case(mobse_ktype_fcms)
				if(sh(i)%mass(idx(i))<mass_brown_dwarf_max)then
					nbd=nbd+1
					bd%mass(nbd)=sh(i)%mass(idx(i))
					bd%age(nbd)=sh(i)%age(idx(i))
				else
					nstar=nstar+1
					star%mass(nstar)=sh(i)%mass(idx(i))
					star%age(nstar)=sh(i)%age(idx(i))
				end if				
            case(mobse_ktype_ms)
                nstar=nstar+1
                star%mass(nstar)=sh(i)%mass(idx(i))
                star%age(nstar)=sh(i)%age(idx(i))
                star%radius(nstar)=sh(i)%radius(idx(i))
            case(mobse_ktype_cowd,mobse_ktype_hewd, mobse_ktype_onwd)
                nwd=nwd+1
                wd%mass(nwd)=sh(i)%mass(idx(i))
                wd%age(nwd)=sh(i)%age(idx(i))
            case(mobse_ktype_ns)
                nns=nns+1
                ns%mass(nns)=sh(i)%mass(idx(i))
                ns%age(nns)=sh(i)%age(idx(i))
            case(mobse_ktype_srg,mobse_ktype_hgap,mobse_ktype_sabr,mobse_ktype_hrbr)
                nrg=nrg+1
                rg%mass(nns)=sh(i)%mass(idx(i))
                rg%age(nns)=sh(i)%age(idx(i))
                rg%radius(nstar)=sh(i)%radius(idx(i))
            case(mobse_ktype_fgb)
                nfgb=nfgb+1
                fgb%mass(nns)=sh(i)%mass(idx(i))
                fgb%age(nns)=sh(i)%age(idx(i))
                fgb%radius(nstar)=sh(i)%radius(idx(i))
            end select			
        end do
    end subroutine
 
	character*(5) function get_kstar_type(ktype)
		implicit none
		integer ktype
		select case(ktype)
		case(0)
			get_kstar_type='FCMS'
		case(1)
			get_kstar_type='MS'
		case(2)
			get_kstar_type='HGap'
		case(3)
			get_kstar_type='FGBr'
		case(4)
			get_kstar_type='HrBr'
		case(5)
			get_kstar_type='RedG'
		case(6)
			get_kstar_type='SABr'
		case(7)
			get_kstar_type='MShe'
		case(8)
			get_kstar_type='HGhe'		
		case(9)
			get_kstar_type='GBhe'
		case(10)
			get_kstar_type='HeWD'
		case(11)
			get_kstar_type='COWD'							
		case(12)
			get_kstar_type='ONWD'
		case(13)
			get_kstar_type='NS'
		case(mobse_ktype_bh)
			get_kstar_type='BH'
		case(15)
			get_kstar_type='MSup'		
		end select
	end function

	character*(9) function get_kw_type(ktype)
		implicit none
		integer ktype
		select case(ktype)
        case(0)
            get_kw_type='NOTHING  '
		case(1)
			get_kw_type='INITIAL  '
		case(2)
           ! Change of stellar type
			get_kw_type='KW_CHANGE'
		case(3)
			get_kw_type='BEG_RCHE '
		case(4)
			get_kw_type='END_RCHE '
		case(5)
            ! Contacted system
			get_kw_type='CONTACT  '
		case(6)
			get_kw_type='COELESCE '
		case(7)
            !Common envelpe evolution
			get_kw_type='COMMONENV'
		case(8)
			get_kw_type='GNTAGE   '		
		case(9)
            !supernova with no remant left, e.g., accretion induced supernova
			get_kw_type='NO_REMNT '
		case(10)
			get_kw_type='MAX_TIME '
		case(11)
            !Binary dissolved by a supernova or tides
			get_kw_type='DISRUPT  '							
		case(12)
            !Begin of Symbiotic-type star
			get_kw_type='BEG_SYMB '
		case(13)
            !End of Symbiotic type star
			get_kw_type='END_SYMB '
		case(14)
            !Begin of Blue Struggler star
			get_kw_type='BEG_BSS  '
        case default
            print*, "unknown kw type:", ktype
		end select
	end function
	integer function get_kstar_integer(str_kstar)
		implicit none
		character*(*) str_kstar
		select case(trim(adjustl(str_kstar)))
		case('FCMS')
			get_kstar_integer=mobse_ktype_fcms   ! deeply or fully convective low mass MS star
		case('MS')
			get_kstar_integer=mobse_ktype_ms   ! Main sequence star
		case('HGap')
			get_kstar_integer=2   ! Hertzsprung Gap
		case('FGBr')
			get_kstar_integer=3   ! First Giant Branch
		case('HrBr')
			get_kstar_integer=4   ! Horizontal Branch / Core Helium Burning
		case('RedG')
			get_kstar_integer=5   ! First Asymptotic Giant Branch / Red Supergiant
		case('SABr')
			get_kstar_integer=6   ! Second Asymptotic Giant Branch 
		case('MShe')
			get_kstar_integer=7   ! Main sequence Naked Helium star
		case('HGhe')
			get_kstar_integer=8   ! Hertzsprung Gap Naked Helium star
		case('GBhe')
			get_kstar_integer=9   ! Giant Branch Naked Helium star
		case('HeWD')
			get_kstar_integer=10  ! Helium White Dwarf
		case('COWD')
			get_kstar_integer=11  ! Carbon / Oxygen White Dwarf
		case('ONWD')
			get_kstar_integer=12  ! Oxygen / Neon White Dwarf
		case('NS')
			get_kstar_integer=13  ! Neutron star
		case('BH')
			get_kstar_integer=mobse_ktype_bh  ! Black hole
		case('MSup')
			get_kstar_integer=15  ! Massless Supernova
		end select
	end function
end module