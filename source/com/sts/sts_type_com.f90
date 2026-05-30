!module md_sts
!	type(stats_1d_type)
!		real(8),allocatable::bin(:),fprb(:),fcmp(:),fncp(:)
!		integer,allocatable::ffra(:),fcfr(:),fncf(:)
!		real(8) mean,std
!	end type
!end module

module com_sts_type
	use md_s1d_type
	use md_s2d_type
	use md_s1d_ird_type
	use md_s1d_hst_type
	use md_s1d_hst_ird_type
	use md_s1d_adv_type
	use md_s2d_ird_type
	use md_s2d_hst_ird_type
	!use md_s1d_hst_type
	!use md_fc_type
end module
