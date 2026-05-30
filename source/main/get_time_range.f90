program get_trange
    use com_main_gw
    use md_dms_saving_data
    implicit none
    integer na
    integer i
    character*(200) fdir,fl
    character*(6) tmpi
    real(8)::rh_now
    real(8),allocatable::cur_time(:)
    fdir="./output/ecev/dms/"
    call get_num_of_runs(fdir,na)

    if (na>0)then
        allocate(cur_time(na))
        call head("model.in","mfrac.in")
        print*, "ctl%gw_radiation_otby=",ctl%gw_radiation_otby
        open(unit=1,file="./output/run_time")
        if(ctl%gw_radiation_otby.eq.0)then           
            
            write(unit=1,fmt="(20A20)") "idx", "time", "mbh", "Mcl", "rh(pc)"!, "reff(pc)"
            do i=1, na
                print*, "i=",i
                write(unit=tmpi,fmt="(I6)") i
                fl=trim(adjustl(fdir))//"dms_"//trim(adjustl(tmpi))
                call input_diffuse_mspec_bin(trim(adjustl(fl)))
                cur_time(i)=ctl%run_snap_time_f

                call get_rh_now(rh_now)
            ! call get_reff_now(nsc_radius_eff)
                write(unit=1,fmt="(I20, 10E20.10)") i, cur_time(i), spp_new%mbh, &
                spp_new%M_r_within_max*m0_cl, 10**rh_now*r0_cl/pc!, 10**nsc_radius_eff*r0_cl/pc
            end do
            
        else
            write(unit=1,fmt="(22A20)") "idx", "time(Myr)", "dt(Myr)","mbh", "Mcl", "rh(pc)", "SBH-emri", &
                "NS-emri", "WD-emri", "BD-emri", "MS-emri", "SBH-emax", "BD-emax", "MS-emax"
            do i=1, na
                print*, "i=",i
                write(unit=tmpi,fmt="(I6)") i
                fl=trim(adjustl(fdir))//"dms_"//trim(adjustl(tmpi))
                call input_diffuse_mspec_bin(trim(adjustl(fl)))
                cur_time(i)=ctl%run_snap_time_f

                call get_rh_now(rh_now)
            ! call get_reff_now(nsc_radius_eff)
                write(unit=1,fmt="(I20, 22E20.10)") i, cur_time(i), ctl%run_snap_time_f-ctl%run_snap_time_i, &
                    spp_new%mbh, &
                spp_new%M_r_within_max*m0_cl, 10**rh_now*r0_cl/pc,&
                oe_sbh%se_emris%rate, oe_ns%se_emris%rate, oe_wd%se_emris%rate, &
                oe_bd%se_emris%rate, oe_star%se_emris%rate, oe_sbh%se_emax%rate, &
                oe_bd%se_emax%rate, oe_star%se_emax%rate
            end do
        end if
        close(1)
    else
        print*, "na=0 do nothing"
    end if

end
