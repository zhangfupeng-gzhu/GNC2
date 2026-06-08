module md_mbh_evl_acc
    use md_mass_bins
    implicit none
    real(8),parameter::td_direct_frac(n_tot_comp_sg)=0d0 !(/0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0/)
    integer,parameter::mmg_num_of_arrays=12
    integer,parameter::mmg_real_length=n_tot_comp_sg*mmg_num_of_arrays+1!+5
    type mass_mbh_growth
        !integer lc_direct_nsample(1:n_tot_comp_sg)
        !integer 
        real(8) lc_direct(1:n_tot_comp_sg)
        real(8) td_direct(1:n_tot_comp_sg)
        real(8) emax_direct(1:n_tot_comp_sg)
        real(8) emax_direct_acum(1:n_tot_comp_sg)
        real(8) lc_direct_acum(1:n_tot_comp_sg)
        real(8) td_direct_acum(1:n_tot_comp_sg)

        real(8) td_disc(1:n_tot_comp_sg)
        real(8) td_disc_acum(1:n_tot_comp_sg) 
        real(8) emri(1:n_tot_comp_sg)
        real(8) emri_acum(1:n_tot_comp_sg)

        real(8) mass_loss_stellar_evolution
        real(8) mass_loss_stellar_evolution_acum

        real(8) gas_reservior_add_acum
        real(8) gas_reservior_left
        real(8) gas_reservior_add ! in this snap
        real(8) gas_reservior_consum ! in this snap
        real(8) mass_direct_swallow_acum
        real(8) mass_direct_swallow !in_this_snap
        real(8) mass_tot
        
    contains
        procedure::init=>init_mass_mbh_growth
        procedure::init_one_snap=>init_mass_mbh_growth_in_one_snap
        procedure::add_event_stype=>mmg_add_base_event_stype
        procedure::pack=>mmg_pack
        procedure::unpack=>mmg_unpack
        
    end type
    private::init_mass_mbh_growth,mmg_add_base_event_stype, mmg_pack, mmg_unpack
    private::init_mass_mbh_growth_in_one_snap
    type(mass_mbh_growth)::mbh_mmg
contains
    subroutine mmg_pack(mmg,fx)
        implicit none
        class(mass_mbh_growth)::mmg
        real(8) fx(mmg_real_length)
        integer i, idx_starter
        idx_starter=0
        fx(idx_starter+1:idx_starter+n_tot_comp_sg)=mmg%lc_direct
        idx_starter=idx_starter+n_tot_comp_sg
        fx(idx_starter+1:idx_starter+n_tot_comp_sg)=mmg%td_direct
        idx_starter=idx_starter+n_tot_comp_sg
        fx(idx_starter+1:idx_starter+n_tot_comp_sg)=mmg%emax_direct
        idx_starter=idx_starter+n_tot_comp_sg
        fx(idx_starter+1:idx_starter+n_tot_comp_sg)=mmg%emax_direct_acum
        idx_starter=idx_starter+n_tot_comp_sg
        fx(idx_starter+1:idx_starter+n_tot_comp_sg)=mmg%lc_direct_acum
        idx_starter=idx_starter+n_tot_comp_sg
        fx(idx_starter+1:idx_starter+n_tot_comp_sg)=mmg%td_direct_acum
        idx_starter=idx_starter+n_tot_comp_sg
        fx(idx_starter+1:idx_starter+n_tot_comp_sg)=mmg%td_disc
        idx_starter=idx_starter+n_tot_comp_sg
        fx(idx_starter+1:idx_starter+n_tot_comp_sg)=mmg%td_disc_acum 
        idx_starter=idx_starter+n_tot_comp_sg
        fx(idx_starter+1:idx_starter+n_tot_comp_sg)=mmg%emri
        idx_starter=idx_starter+n_tot_comp_sg
        fx(idx_starter+1:idx_starter+n_tot_comp_sg)=mmg%emri_acum

        idx_starter=idx_starter+n_tot_comp_sg
        fx(idx_starter+1)=mmg%mass_loss_stellar_evolution 
    end subroutine
    subroutine mmg_unpack(mmg,fx)
        implicit none
        class(mass_mbh_growth)::mmg
        real(8) fx(mmg_real_length)
        integer i, idx_starter
        idx_starter=0
        mmg%lc_direct=fx(idx_starter+1:idx_starter+n_tot_comp_sg)

        idx_starter=idx_starter+n_tot_comp_sg
        mmg%td_direct=fx(idx_starter+1:idx_starter+n_tot_comp_sg)

        idx_starter=idx_starter+n_tot_comp_sg
        mmg%emax_direct=fx(idx_starter+1:idx_starter+n_tot_comp_sg)

        idx_starter=idx_starter+n_tot_comp_sg
        mmg%emax_direct_acum=fx(idx_starter+1:idx_starter+n_tot_comp_sg)

        idx_starter=idx_starter+n_tot_comp_sg
        mmg%lc_direct_acum=fx(idx_starter+1:idx_starter+n_tot_comp_sg)

        idx_starter=idx_starter+n_tot_comp_sg
        mmg%td_direct_acum=fx(idx_starter+1:idx_starter+n_tot_comp_sg)

        idx_starter=idx_starter+n_tot_comp_sg
        mmg%td_disc=fx(idx_starter+1:idx_starter+n_tot_comp_sg)

        idx_starter=idx_starter+n_tot_comp_sg
        mmg%td_disc_acum=fx(idx_starter+1:idx_starter+n_tot_comp_sg) 

        idx_starter=idx_starter+n_tot_comp_sg
        mmg%emri=fx(idx_starter+1:idx_starter+n_tot_comp_sg)

        idx_starter=idx_starter+n_tot_comp_sg
        mmg%emri_acum=fx(idx_starter+1:idx_starter+n_tot_comp_sg)

        idx_starter=idx_starter+n_tot_comp_sg
        mmg%mass_loss_stellar_evolution=fx(idx_starter+1) 
        
    end subroutine
    
    subroutine init_mass_mbh_growth(mmg)
        implicit none
        class(mass_mbh_growth)::mmg
        call init_mass_mbh_growth_in_one_snap(mmg)
        mmg%lc_direct_acum=0; mmg%td_direct_acum=0; 
        mmg%td_disc_acum=0; 
        mmg%emri_acum=0; mmg%emax_direct_acum=0
        mmg%gas_reservior_add_acum=0
        mmg%mass_direct_swallow_acum=0
        mmg%gas_reservior_consum=0
        mmg%mass_loss_stellar_evolution_acum=0
    end subroutine
    subroutine init_mass_mbh_growth_in_one_snap(mmg)
        implicit none
        class(mass_mbh_growth)::mmg
        mmg%lc_direct=0; mmg%td_direct=0; mmg%td_direct=0;  
        mmg%emri=0; mmg%emax_direct=0; 
        mmg%gas_reservior_add=0
        mmg%mass_direct_swallow=0
        mmg%mass_tot=0
        mmg%mass_loss_stellar_evolution=0
    end subroutine
    subroutine update_acum_data(mmg)
        implicit none
        type(mass_mbh_growth)::mmg
        mmg%lc_direct_acum=mmg%lc_direct_acum+mmg%lc_direct
        mmg%td_direct_acum=mmg%td_direct_acum+mmg%td_direct
        mmg%td_disc_acum=mmg%td_disc_acum+mmg%td_disc
        mmg%td_direct_acum=mmg%td_direct_acum+mmg%td_direct 
        mmg%emax_direct_acum=mmg%emax_direct_acum+mmg%emax_direct
        mmg%emri_acum=mmg%emri_acum+mmg%emri
        mmg%mass_loss_stellar_evolution_acum=mmg%mass_loss_stellar_evolution_acum &
            +mmg%mass_loss_stellar_evolution
    end subroutine
    subroutine update_tot_data(mmg)
        implicit none
        type(mass_mbh_growth)::mmg

        mmg%gas_reservior_add=sum(mmg%td_disc)+ mmg%mass_loss_stellar_evolution
        mmg%gas_reservior_add_acum=mmg%gas_reservior_add_acum+mmg%gas_reservior_add
        mmg%mass_direct_swallow=sum(mmg%lc_direct)+sum(mmg%td_direct)+sum(mmg%emax_direct)+sum(mmg%emri)
        mmg%mass_direct_swallow_acum=mmg%mass_direct_swallow_acum+ mmg%mass_direct_swallow
    end subroutine
    subroutine write_hdf5_mass_mbh_growth(group_id,mmg)
        use md_hdf5
        implicit none
        type(mass_mbh_growth)::mmg
        integer(HID_T)::group_id,attr_id  
        call add_attr_dble(group_id,attr_id,"dmbh",mmg%mass_tot)
        call add_attr_dble(group_id,attr_id,"gas_reservior_consum",mmg%gas_reservior_consum)
        call add_attr_dble(group_id,attr_id,"gas_reservior_left",mmg%gas_reservior_left)
        call add_attr_dble(group_id,attr_id,"gas_reservior_add",mmg%gas_reservior_add)
        call add_attr_dble(group_id,attr_id,"mass_direct_swallow",mmg%mass_direct_swallow)
        call add_attr_dble(group_id,attr_id,"mass_direct_swallow_acum",mmg%mass_direct_swallow_acum)
        call add_attr_dble(group_id,attr_id,"gas_reservior_add_acum",mmg%gas_reservior_add_acum)
        call add_attr_dble(group_id,attr_id,"mass_loss_stellar_evolution",mmg%mass_loss_stellar_evolution)
        call add_attr_dble(group_id,attr_id,"mass_loss_stellar_evolution_acum",mmg%mass_loss_stellar_evolution_acum)

        call add_attr_dble_arr(group_id,attr_id,"emax_direct",mmg%emax_direct, n_tot_comp_sg)
        call add_attr_dble_arr(group_id,attr_id,"emax_direct_acum",mmg%emax_direct_acum, n_tot_comp_sg)
        call add_attr_dble_arr(group_id,attr_id,"emri",mmg%emri, n_tot_comp_sg)
        call add_attr_dble_arr(group_id,attr_id,"emri_acum",mmg%emri_acum, n_tot_comp_sg)

        call add_attr_dble_arr(group_id,attr_id,"td_direct", &
            mmg%td_direct, n_tot_comp_sg)
        call add_attr_dble_arr(group_id,attr_id,"td_direct_acum", &
            mmg%td_direct_acum, n_tot_comp_sg)
        call add_attr_dble_arr(group_id,attr_id,"lc_direct", &
            mmg%lc_direct, n_tot_comp_sg)
        call add_attr_dble_arr(group_id,attr_id,"lc_direct_acum", &
            mmg%lc_direct_acum, n_tot_comp_sg)
        call add_attr_dble_arr(group_id,attr_id,"td_disc", &
            mmg%td_disc, n_tot_comp_sg)
        call add_attr_dble_arr(group_id,attr_id,"td_disc_acum", &
            mmg%td_disc_acum, n_tot_comp_sg) 
        
    end subroutine
    subroutine mmg_add_base_event_stype(mmg,sp,mass)
        use md_particle_sample
        implicit none
        class(mass_mbh_growth)::mmg
        type(particle_sample_type)::sp
        real(8) mass
        integer idx
        call get_obidx_from_type_sg(sp%obtype,idx)
        select case(sp%exit_flag)
        case(exit_tidal_empty,exit_tidal_full)
            mmg%td_direct(idx)=mmg%td_direct(idx)+mass*td_direct_frac(idx)
            mmg%td_disc(idx)=mmg%td_disc(idx)+mass*(1-td_direct_frac(idx))
        case(exit_lc)
            mmg%lc_direct(idx)=mmg%lc_direct(idx)+mass
        case(exit_boundary_max)
            if(sp%state_emri_current.ge.1.and.sp%state_emri_last.ge.1)then
                mmg%emri(idx)=mmg%emri(idx)+mass
            else
                mmg%emax_direct(idx)=mmg%emax_direct(idx)+mass
            end if
        case(exit_emri_single)
            mmg%emri(idx)=mmg%emri(idx)+mass
        end select
    end subroutine
    subroutine write_mmg(funit)
        implicit none
        !type(mass_mbh_growth)::mmg
        integer funit
        write(funit) &
        mbh_mmg%lc_direct(1:n_tot_comp_sg),&
        mbh_mmg%td_direct(1:n_tot_comp_sg),&
        mbh_mmg%emax_direct(1:n_tot_comp_sg),&
        mbh_mmg%emax_direct_acum(1:n_tot_comp_sg),&
        mbh_mmg%lc_direct_acum(1:n_tot_comp_sg),&
        mbh_mmg%td_direct_acum(1:n_tot_comp_sg),&
        mbh_mmg%td_disc(1:n_tot_comp_sg),&
        mbh_mmg%td_disc_acum(1:n_tot_comp_sg),& 
        mbh_mmg%emri(1:n_tot_comp_sg),&
        mbh_mmg%emri_acum(1:n_tot_comp_sg),&
        mbh_mmg%mass_loss_stellar_evolution,&
        mbh_mmg%mass_loss_stellar_evolution_acum,&
        mbh_mmg%gas_reservior_add_acum,&
        mbh_mmg%gas_reservior_left,&
        mbh_mmg%gas_reservior_add,&
        mbh_mmg%gas_reservior_consum,&
        mbh_mmg%mass_direct_swallow_acum,&
        mbh_mmg%mass_direct_swallow,&
        mbh_mmg%mass_tot
    end subroutine
    subroutine read_mmg(funit)
        implicit none
        !type(mass_mbh_growth)::mmg
        integer funit
        read(funit) &
        mbh_mmg%lc_direct(1:n_tot_comp_sg),&
        mbh_mmg%td_direct(1:n_tot_comp_sg),&
        mbh_mmg%emax_direct(1:n_tot_comp_sg),&
        mbh_mmg%emax_direct_acum(1:n_tot_comp_sg),&
        mbh_mmg%lc_direct_acum(1:n_tot_comp_sg),&
        mbh_mmg%td_direct_acum(1:n_tot_comp_sg),&
        mbh_mmg%td_disc(1:n_tot_comp_sg),&
        mbh_mmg%td_disc_acum(1:n_tot_comp_sg),& 
        mbh_mmg%emri(1:n_tot_comp_sg),&
        mbh_mmg%emri_acum(1:n_tot_comp_sg),&
        mbh_mmg%mass_loss_stellar_evolution,&
        mbh_mmg%mass_loss_stellar_evolution_acum,&
        mbh_mmg%gas_reservior_add_acum,&
        mbh_mmg%gas_reservior_left,&
        mbh_mmg%gas_reservior_add,&
        mbh_mmg%gas_reservior_consum,&
        mbh_mmg%mass_direct_swallow_acum,&
        mbh_mmg%mass_direct_swallow,&
        mbh_mmg%mass_tot
    end subroutine
end module