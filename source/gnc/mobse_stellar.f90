module md_mobse_stellar_single
    use com_sts_type
    implicit none
     

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
contains 
 
 
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