module md_iso
	use com_sts_type
	use constant
	type ISO_table_type
		type(s2d_type)::rpgr_ecgr_incgr, rpnw_ecgr_incgr,rpnw_acnw_incgr
		type(s2d_type)::incnw_ecgr_incgr, coitagr_ecgr_incgr, acnw_ecgr_incgr
 		type(s1d_type)::rpnw_incgr
		real(8) spin
		integer nx, ny 
	end type 
contains 
 	subroutine input_rp_iso(iso_table, fl)
		!use com_main_gw
		implicit none
		character*(*) fl
		type(ISO_table_type)::iso_table
		open(unit=199992,file=trim(adjustl(fl))//".bin",form="unformatted",status="old")
		read(unit=199992) iso_table%rpgr_ecgr_incgr,iso_table%rpnw_ecgr_incgr, iso_table%rpnw_acnw_incgr, &
			iso_table%incnw_ecgr_incgr,iso_table%coitagr_ecgr_incgr,iso_table%acnw_ecgr_incgr,iso_table%rpnw_incgr
		close(unit=199992)
	end subroutine 
end module