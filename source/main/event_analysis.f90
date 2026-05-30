program emris_analysis
    use com_main_gw
    use md_event_datas
    implicit none
    character*(200) fl
    character*(200) args(3)
    integer i, ibg, ied, istep
    character*(6) tmpi 

    call head("model.in","mfrac.in")
    !call init_mpi()
    
    call get_command_argument(1,args(1))  ! ibg
    if(args(1).eq.'')then
        open(unit=2,file="./output/run_summary.txt")
		read(unit=2,fmt=*) ied
		close(unit=2)
        ibg=1; istep=1
    else
        read(unit=args(1), fmt=*) ibg
        call get_command_argument(2,args(2))  ! ied
        call get_command_argument(3,args(3))  ! istep
        read(unit=args(2), fmt=*) ied
        read(unit=args(3), fmt=*) istep
    end if
    print*, "ibg,ied,istep=",ibg,ied,istep

    do i=ibg, ied, istep
        print*, "i=",i
        write(unit=tmpi,fmt="(I6)") i
        fl="./output/ecev/bin/event_data/"//trim(adjustl(tmpi))
        call input_all_emris_datas(fl)
        call input_diffuse_mspec_bin("./output/ecev/dms/dms_"//trim(adjustl(tmpi)))
        print*, "sbh"
        call output_event_datas_hdf5(data_sbh, &
        "./output/ecev/dms/event_data/"//trim(adjustl(tmpi))//"_sbh_event",star_type_bh)
        print*, "ns"
        call output_event_datas_hdf5(data_ns, &
        "./output/ecev/dms/event_data/"//trim(adjustl(tmpi))//"_ns_event",star_type_ns)
        print*, "wd"
        call output_event_datas_hdf5(data_wd, &
        "./output/ecev/dms/event_data/"//trim(adjustl(tmpi))//"_wd_event",star_type_wd)
        print*, "star"
        call output_event_datas_hdf5(data_star, &
        "./output/ecev/dms/event_data/"//trim(adjustl(tmpi))//"_star_event",star_type_ms)
        print*, "bd"
        call output_event_datas_hdf5(data_bd, &
        "./output/ecev/dms/event_data/"//trim(adjustl(tmpi))//"_bd_event",star_type_bd)
        print*, "rg"
        call output_event_datas_hdf5(data_rg, &
        "./output/ecev/dms/event_data/"//trim(adjustl(tmpi))//"_rg_event",star_type_rg)
    end do
    call end()
end