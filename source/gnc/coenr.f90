

subroutine get_coenr(even, evjum, m, en, jc, coenr,idx,idy)
    use com_main_gw
    implicit none
    real(8) m
    type(coeff_type)::coenr
    real(8) even, evjum
    real(8) dj_111(20), de_110(20)
    integer idx, idy
    real(8) jc, de_0, dee,dj_rest,djj,dej, en, jmin, jmax!, rdx,rdy
    integer i, nbin
    
    
    select case(ctl%method_interpolate)
    case(method_int_nearst)
        
        do i=1, dms%n
            !evjum=dms%dc0%s2_de_110%ymin            
            de_110(i)=dms%mb(i)%dc%s2_de_110%fxy(idx,idy)
            dj_111(i)=dms%mb(i)%dc%s2_dj_111%fxy(idx,idy)
        end do
        associate(dc0=>dms%dc0)
        !    dc0=>dms%dc0
            de_0=dc0%s2_de_0%fxy(idx,idy)
            dee=dc0%s2_dee%fxy(idx,idy)
            dj_rest=dc0%s2_dj_rest%fxy(idx,idy)
            djj=dc0%s2_djj%fxy(idx,idy)
            dej=dc0%s2_dej%fxy(idx,idy)
            if(djj<0)then
                print*, "what happened???"
                print*, "djj=",djj
                call dc0%s2_djj%print("s2_djj")
                stop
            end if
        end associate
    case(method_int_linear)
        do i=1, dms%n
            !print*, "dms%df_coe_bins=",dms%df_coe_bins
            call linear_int_2d_xy(idx,idy,sample_table_rdx,sample_table_rdy,&
                dms%mb(i)%dc%s2_de_110%fxy,dms%df_coe_bins,dms%df_coe_bins,de_110(i))
            call linear_int_2d_xy(idx,idy,sample_table_rdx,sample_table_rdy,&
                dms%mb(i)%dc%s2_dj_111%fxy,dms%df_coe_bins,dms%df_coe_bins,dj_111(i))
        end do
        associate(dc0=>dms%dc0)
            call linear_int_2d_xy(idx,idy,sample_table_rdx,sample_table_rdy,&
                dc0%s2_de_0%fxy,dms%df_coe_bins,dms%df_coe_bins,de_0)

            !call linear_int_2d(dms%logemin,dms%jmin,dms%df_coe_bins,dms%df_coe_bins,dc_grid_xstep,dc_grid_ystep,&
            !dc0%s2_dee%fxy, even,evjum,dee)
            !print*, "dee=",dee
            !==test==
            !    block
            !            integer i, j
            !            real(8) xt,yt
            !            xt=dms%df_coe_bins+0.1; yt=dms%df_coe_bins+0.1
            !            do i=1, dms%df_coe_bins
            !                do j=1,dms%df_coe_bins
            !                common_dee_log%fxy(i,j)=i+j-1
            !                end do
            !            end do
            !            print*, "idx,idy,rdx,rdy=",int(xt)+1,int(yt)+1,xt,yt
            !            call linear_int_2d_xy(dms%df_coe_bins,dms%df_coe_bins,xt,yt,&
            !    common_dee_log%fxy,dms%df_coe_bins,dms%df_coe_bins,dee)
            !            print*, "dee=",dee, xt+yt
            !            stop
            !    end block
            !==
            call linear_int_2d_xy(idx,idy,sample_table_rdx,sample_table_rdy,&
                common_dee_log%fxy,dms%df_coe_bins,dms%df_coe_bins,dee)
            !print*, "idx,idy,rdx,rdy,dee=",idx,idy,rdx,rdy,dee

            dee=10**dee
            !print*, "dee=",dee
            call linear_int_2d_xy(idx,idy,sample_table_rdx,sample_table_rdy,&
                dc0%s2_dj_rest%fxy,dms%df_coe_bins,dms%df_coe_bins,dj_rest)
        ! call linear_int_2d(dms%logemin,dms%jmin,dms%df_coe_bins,dms%df_coe_bins,dc_grid_xstep,dc_grid_ystep,&
            !dc0%s2_djj%fxy, even,evjum,djj)
            !print*, "djj=",djj
            call linear_int_2d_xy(idx,idy,sample_table_rdx,sample_table_rdy,&
                common_djj_log%fxy,dms%df_coe_bins,dms%df_coe_bins,djj)
            djj=10**djj
            !print*, "djj=",djj
            !read(*,*)
            call linear_int_2d_xy(idx,idy,sample_table_rdx,sample_table_rdy,&
                dc0%s2_dej%fxy,dms%df_coe_bins,dms%df_coe_bins,dej)    
        end associate
    !end if
        !if(idx.eq.nbin.or.idx.eq.1)then
        !    call common_jc%print("jc")
        !    print*, "rdx,idx,even, common_jc%xb=", rdx,idx,even, common_jc%xb(idx)
        !    print*, "de_110,dj_111,de_0,dee,dj_rest, djj,dej=", de_110(1),dj_111(1),de_0,dee,dj_rest, djj,dej
        !    read(*,*)
        !end if
    
    !read(*,*)
    end select

    !sample%en=-mbh/(2*sample%byot%a_bin)
    select case (ctl%Dejmodel)
    case(dejmodel_EJ)
        coeNr%jj=djj*jc*jc; 
        coeNr%e=de_0; coeNr%j=dj_rest
        do i=1, dms%n
            coeNr%e=coeNr%e+m/dms%mb(i)%mc*de_110(i)
            coeNr%j=coeNr%j+dj_111(i)*(m+dms%mb(i)%mc)/dms%mb(i)%mc/2d0; 
        end do
        coeNr%ee=dee*en*en;
        !coeNr%ee=dee*(10**even*ctl%energy0)**2

        coeNr%e=coeNr%e*en; 
        !coeNr%e=coeNr%e*(10**even*ctl%energy0)

        coeNR%j=  coeNR%j*jc
        coeNr%ej=  dej*en*jc
        if(ctl%gw_radiation_otby.ge.1)then
            sample_rlx_e_time=1d0/dee
            sample_rlx_j_time=10**(evjum*2)/djj
        end if
    case(dejmodel_xj)
        coeNr%jj=djj; 
        coeNr%e=de_0; coeNr%j=dj_rest
        do i=1, dms%n
            coeNr%e=coeNr%e+m/dms%mb(i)%mc*de_110(i)
            coeNr%j=coeNr%j+dj_111(i)*(m+dms%mb(i)%mc)/dms%mb(i)%mc/2d0; 
        end do
        coeNr%ee=dee*even*even;
        !coeNr%ee=dee*(10**even*ctl%energy0)**2

        coeNr%e=coeNr%e*even; 
        coeNr%ej=  dej*even
    case default
        print*, "error! define dejmodel", ctl%dejmodel
        stop
    end select
    
    if(ctl%chattery.ge.4)then
        print*, "=========get cej NR========================="
        print*, "enev, evjum, en, jc=", even, evjum, en, jc
        print*, "idx,idy=",idx,idy
        print*, "de_0, dee=", de_0, dee
        print*, "djj=", djj
        print*, "coenr%e, ee, ej=", coenr%e, coenr%ee,coenr%ej
        print*, "coenr%j, jj=", coenr%j, coenr%jj
        if(ctl%method_interpolate.eq.method_int_linear)then
            print*, "rdx,rdy,xstep,ystep,xmin,ymin=",sample_table_rdx,sample_table_rdy,&
                dc_grid_xstep,dc_grid_ystep,&
                log10emin_factor,log10jmin_value
        end if
        !print*, "=========end of get cej NR========"
    end if

end subroutine
 
subroutine get_steps_nr_EJ(en, jm,coenr, jc,time_dt_e, time_dt_j)
    use com_main_gw
    implicit none
    real(8) steps,period, jc, jm,en, enev
    type(coeff_type)::coenr
    real(8) time_dt_e, time_dt_j!, time_dt_nr
    if(coenr%ee.ne.0)then
        time_dt_e=min((en*0.15)**2/coenr%ee, abs(en*0.15)/abs(coenr%e))
        !if(mbh_dmless.eq.0)then
            !if(en/ctl%energy_max>0.95)then
                !print*, "time_dt_e:i=",time_dt_e,(en*0.15)**2/coenr%ee,abs(en*0.15)/abs(coenr%e)
                !time_dt_e=min(time_dt_e,((1.01-en/ctl%energy_max)*en)**2/coenr%ee)
                !print*, "time_dt_e:f=",time_dt_e, ((1.01-en/ctl%energy_max)*en)**2/coenr%ee
                !read(*,*)
            !end if
        !end if
        
    else
        time_dt_e=1d6!*period
    end if
    if(coenr%jj.ne.0)then
        time_dt_j=min((min(jc*jm*2d0,jc*0.1))**2/coenr%jj, &
        (0.4d0*(1.0075-jm)*jc)**2/coenr%jj)
    else
        time_dt_j=1d6!*period
    end if
    !time_dt_nr=min(time_dt_e,time_dt_j)
    if(ctl%chattery.ge.4)then
        print*, "====get steps NR==========================="
        print*, "en, jm=", en, jm
        print*, "time_dt_e, time_dt_j=", time_dt_e, time_dt_j
        !print*, "=========end of get steps NR========"
    end if
end subroutine

subroutine get_steps_nr_xj(en, jm,coenr, time_dt_nr)
    use com_main_gw
    implicit none
    real(8) steps,period,  jm,en, enev
    type(coeff_type)::coenr
    real(8) time_dt_e, time_dt_j, time_dt_nr

    if(coenr%ee.ne.0)then
        enev=en/ctl%energy0
        time_dt_e=min((enev*0.15)**2/coenr%ee, abs(enev*0.15)/abs(coenr%e))
    else
        time_dt_e=1d6!*period
    end if
    if(coenr%jj.ne.0)then
        time_dt_j=min(0.1d0**2/coenr%jj, &
        (0.4d0*(1.0075-jm))**2/coenr%jj, &
        (0.25*abs(jm))**2/coenr%jj)
    else
        time_dt_j=1d6!*period
    end if
    
    time_dt_nr=min(time_dt_e,time_dt_j)

    if(ctl%chattery.ge.4)then
        print*, "=========get steps NR==============="
        print*, "en, jm=", en, jm
        print*, "time_dt_nr=", time_dt_nr
        print*, "=========end of get steps NR========"
    end if
end subroutine

subroutine get_de_dj_nr(coenr, dt, steps, jum, den, djp,djp0)
    use com_main_gw
    implicit none
    real(8),intent(out):: den, djp, djp0
    type(coeff_type)::coenr
    real(8) jum, rho, y1, y2, dt,n2j,steps
    real(8) gen_gaussian, y3,y4

    !jum=sample%jm*sqrt(mbh*sample%byot%a_bin)
    rho=coeNR%ej/sqrt(abs(coeNR%ee*coeNR%jj))  
    !print*, "jm,rho=", sample%jm, rho
    !if(abs(rho)>=1)then
    !    print*, "rho,ej,ee,jj=", rho, coeNR%ej, coeNR%ee, coeNR%jj
    !    print*, "sample%idx,idy=",sample_table_idx,sample_table_idy
    !    block 
    !        integer idx, idy
    !        idx=sample_table_idx
    !        idy=sample_table_idy
    !        print*, "ej, ee, jj=", dms%dc0%s2_dej%fxy(idx,idy), &
    !        dms%dc0%s2_dee%fxy(idx,idy), dms%dc0%s2_djj%fxy(idx,idy)
    !        print*, "ej/(ee*jj)**0.5=", dms%dc0%s2_dej%fxy(idx,idy)/ &
    !            sqrt(dms%dc0%s2_dee%fxy(idx,idy)*dms%dc0%s2_djj%fxy(idx,idy))
    !        stop
    !    end block
    !endif
    call gen_gaussian_correlate(y1,y2,rho)
    y1=max(min(y1,6d0),-6d0)
    y2=max(min(y2,6d0),-6d0)
    den=coeNR%e*dt+y1*sqrt(coeNR%ee*dt)
    !sample_den=coeNR%e*dt+y1*sqrt((coeNR%ee+coeNr%e**2)*dt)
    n2j=sqrt(coeNR%jj*dt)
    if(n2j<jum/4d0)then
        djp=coeNR%j*dt+y2*n2j
    !	print*, "jum=",jum
    !    print*, "1:djp=",djp, n2j, dt, y2,coeNR%jj
    else
        y3=gen_gaussian(1d0)
        y4=gen_gaussian(1d0)
        djp=sqrt((jum+n2j*y3)**2+(n2j*y4)**2)-jum
        !print*, "1:dsp=",coeNR%j*dt+y2*n2j
        !print*, "2:djp=",djp, y3,y4, jum, n2j
        !read(*,*)
    !    !print*, sample_djp, jum, sample%jm
    !    !read(*,*)
    end if
    
    if(ctl%chattery.ge.4)then
        print*, "======start dedj==========================="
        print*, "dedrift, scatter=", coenr%e*dt, y1*sqrt(coeNR%ee*dt)
        print*, "djdrift, scatter=", coenr%j*dt, y2*n2j
        print*, "den,djp=", den, djp
    end if
    djp0=coeNR%j*dt/steps+y2*sqrt(coeNR%jj*dt/steps)
end subroutine


