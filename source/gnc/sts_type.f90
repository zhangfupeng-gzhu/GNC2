module md_sts
	use com_sts_type
	implicit none
    type sts_single_type
        type(sts_fc_type)::mtot
    end type
	type, extends(sts_single_type):: sts_basic_type
		type(sts_fc_type):: fabin, febin, fac, fec, frpout, fac_rsw,frpout_rsw, fr_rsw, fr
		type(sts_fc_type):: fe_acc, fvc, fdop, m2, m1
	contains
	 	procedure::write=>write_sts
		procedure::read=>read_sts
!		procedure::write_Hdf5
	end type
    type(sts_single_type)::sts_sbh_mgene
	type(sts_basic_type):: sts_merge
	private::write_sts, read_sts
contains

	subroutine write_sts(sts, fl)
		implicit none
		character*(*) fl
       class(sts_basic_type), intent(in) :: sts
	   open(unit=999,file=trim(adjustl(fl))//"_sts.bin",form='unformatted',access='stream')
      	 write(999) sts%fabin, sts%febin, sts%fac, sts%fec, sts%frpout, sts%fac_rsw,sts%frpout_rsw, &
				sts%fr, sts%fr_rsw, sts%fe_acc, sts%fvc, sts%fdop
		close(999)
	end subroutine
	subroutine read_sts(sts, fl)
		implicit none
		character*(*) fl
       class(sts_basic_type), intent(inout) :: sts
	   open(unit=999,file=trim(adjustl(fl))//"_sts.bin",form='unformatted',access='stream',status='old')
      	 read(999) sts%fabin, sts%febin, sts%fac, sts%fec, sts%frpout, sts%fac_rsw,sts%frpout_rsw, &
				sts%fr, sts%fr_rsw, sts%fe_acc, sts%fvc, sts%fdop
		close(999)
	end subroutine
end module
