
module md_dms_saving_data
	use com_sts_type
	implicit none
	type snap_event
		real(8)::n=0
		real(8)::nw=0
		real(8)::rate=0   ! in unit of Myr^{-1}
		type(s1d_hst_type)::fdstr_x
		type(s1d_hst_type)::fdstr_m  ! distribution of mass
		contains
		procedure::write_bin_se
		procedure::read_bin_se
		generic :: write(unformatted) => write_bin_se
		generic :: read(unformatted) => read_bin_se
		
	end type
	type obj_events_basic
		!integer type_idx
		type(snap_event)::se_emax
		type(snap_event)::se_emris
		type(snap_event)::se_lc  ! loss cone
		!type(snap_event)::se_stcold
	end type
	type,extends(obj_events_basic)::obj_events
		type(sts_fc_type)::fd_emris_ecc  ! distribution of eccentricity at 1mHz
		type(s2d_hst_ird_type)::fd_emris_nxj_ir
		type(snap_event)::se_td ! tidal disruption
		contains
		procedure::write_bin_dms_saving_data
		procedure::read_bin_dms_saving_data
		generic :: write(unformatted) => write_bin_dms_saving_data
		generic :: read(unformatted) => read_bin_dms_saving_data
		! contains
		! procedure::write_bin=>write_bin_dms_saving_data
		! procedure::read_bin=>read_bin_dms_saving_data
		
	end type
	! type,extends(obj_events_compact)::obj_events_star
		
	! end type
	type(obj_events)::oe_star
	type(obj_events)::oe_rg
	type(obj_events)::oe_sbh
	type(obj_events)::oe_wd
	type(obj_events)::oe_ns
	type(obj_events)::oe_bd	

	private::write_bin_dms_saving_data,read_bin_dms_saving_data
	private::read_bin_se,write_bin_se
contains 
	subroutine collection_snap_events(se,dt)
		implicit none
		class(snap_event)::se
		real(8) nr, nw,dt
		nr=se%n
		nw=se%nw
		call collection_event_numbers(nr,nw)
		se%n=nr
		se%nw=nw
		se%rate=nw/dt
	end subroutine
	subroutine write_bin_se(se,funit,iostat,iomsg)
		class(snap_event), intent(in) ::se
		integer, intent(in) :: funit
		integer, intent(out) :: iostat
		character(*), intent(inout) :: iomsg
		write(funit, iostat=iostat, iomsg=iomsg)se%n,se%nw,se%rate
		if(se%nw>0)then
			write(funit, iostat=iostat, iomsg=iomsg)&
			se%fdstr_m,se%fdstr_x
		end if
	end subroutine
	subroutine read_bin_se(se,funit,iostat,iomsg)
		class(snap_event), intent(inout) ::se
		integer, intent(in) :: funit
		integer, intent(out) :: iostat
		character(*), intent(inout) :: iomsg
		read(funit, iostat=iostat, iomsg=iomsg)se%n,se%nw,se%rate
		if(se%nw>0)then 
			read(funit, iostat=iostat, iomsg=iomsg) &
			se%fdstr_m,se%fdstr_x
		end if
	end subroutine
	subroutine write_bin_dms_saving_data(oe,funit,iostat,iomsg)
		implicit none
		class(obj_events), intent(in) ::oe
		integer, intent(in) :: funit
		integer, intent(out) :: iostat
		character(*), intent(inout) :: iomsg

		write(funit, iostat=iostat, iomsg=iomsg) oe%se_emax,oe%se_emris,oe%se_lc,oe%se_td
		!print*, "write dms saving:1"
		write(funit, iostat=iostat, iomsg=iomsg) oe%fd_emris_ecc,oe%fd_emris_nxj_ir
		!print*, "write dms saving:2"
	end subroutine
	subroutine read_bin_dms_saving_data(oe,funit,iostat,iomsg)
		implicit none
		class(obj_events), intent(inout) ::oe
		integer, intent(in) :: funit
		integer, intent(out) :: iostat
		character(*), intent(inout) :: iomsg
		! print*, "1"
		read(funit, iostat=iostat, iomsg=iomsg) oe%se_emax,oe%se_emris,oe%se_lc,oe%se_td
		! print*, "2"
		read(funit, iostat=iostat, iomsg=iomsg) oe%fd_emris_ecc,oe%fd_emris_nxj_ir
		! print*, "3"
	end subroutine
end module