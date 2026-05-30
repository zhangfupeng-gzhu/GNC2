subroutine get_fden_ir_2d(dm,spp,gxj_ir,fden,rho_rmin)
    use constant
	use com_sts_type
	use my_intgl
    use md_dms
    use model_basic,only:ctl,sample_logrmin
    use md_star_pot
    use md_coeff
    use MPI_comu,only:rid,mpi_comm_world
	use, intrinsic :: ieee_arithmetic
    implicit none
    type(diffuse_mspec)::dm
	type(s1d_type)::fden
    type(s2d_type)::gxj_ir
	real(8)  rho_rmin
	integer i
    integer::debug=0
    type(star_pot_para)::spp

	do i=1, fden%nbin
		call get_fden_ird_2d_one(fden%xb(i),fden%fx(i),gxj_ir,dm,spp)
	end do
    !call fden%print("fden")
    
    call get_fden_ird_2d_one(sample_logrmin,rho_rmin,gxj_ir,dm,spp)
    !print*, "rho_rmin=",rho_rmin
    !read(*,*)
end subroutine
subroutine get_fden_ird_2d_one(logr,fden,gxj_ir,dm,spp)
    use com_main_gw
    implicit none
    real(8) logr, fden
    integer flag_out
    type(s2d_type)::gxj_ir
    type(diffuse_mspec)::dm
    type(star_pot_para)::spp
    type(s1d_ird_type)::common_gx_ir    
    interface 
        subroutine get_fden_ird_one(logr, fx,common_gx_ir,dm,spp,rho_max,flag_out)
            use com_main_gw
            implicit none
            real(8) logr, fx,rho_max
            type(s1d_ird_type),target::common_gx_ir
            type(diffuse_mspec)::dm
            type(star_pot_para)::spp
            integer flag_out
        end subroutine
    end interface

    call common_gx_ir%init(sample_logemin,sample_logemax,ctl%dstr_bins_e,sts_type_dstr)
    common_gx_ir%xb=gxj_ir%xcenter
    common_gx_ir%nbin=gxj_ir%nx

    !call get_gxr_bak(common_gx_ir,gxj_ir,logr,spp,dm,.false.)
    !call common_gx_ir%print('common_gx_ir_bak')
    !ctl%debug=1
    call get_gxr(common_gx_ir,gxj_ir,logr,spp,dm,.false.)
    !call common_gx_ir%print('common_gx_ir')
    !read(*,*)
    call get_fden_ird_one(logr,fden,common_gx_ir,dm,spp,0d0,flag_out)
    !print*, "logr,fden=",logr,fden
    !read(*,*)
end subroutine
subroutine get_gxr_bak(gxr,gxj_ir,logr,spp,dm,jc_pre)
    use com_main_gw
    implicit none
    type(s1d_ird_type)::gxr
    type(s2d_type)::gxj_ir
    type(diffuse_mspec)::dm
    type(star_pot_para)::spp
    logical::jc_pre
    real(8) t0,r,ex,logphi,jx,theta
    real(8) logex, logr,jc_xy
    real(8) rc, r_c_iter,jc_dmless
    real(8) theta_max,theta_min
    integer ier, i,j,idx,n
    integer,parameter::nbin_size=300
    !real(8),parameter:: yp1=1d30, ypn=1d30
    !real(8) y2(gxj_ir%ny),yout,y(gxj_ir%ny), jx1, jmid,jx0,theta1,theta0,logt0
    real(8) jx1, jmid,jx0,theta1,theta0,logt0
    real(8) sintheta1,sintheta0,costheta1,costheta0


    gxr%fx=0
    r=10**logr
    do i=1, gxr%nbin
        ex=10**gxj_ir%xcenter(i)
        call get_phi_star_full_range(spp,logr,logphi)
        if(jc_pre)then
            jc_xy=dms%jc_sample_erange%fx(i)
            
            !print*, "1, xb, jc_xy=",10**dms%jc_sample_erange%xb(i),jc_xy
            !rc=r_c_iter(spp,ex,ier)
            !jc_xy=jc_dmless(rc,spp)
           ! print*, "2, ex, jc_xy=",ex, jc_xy
           !read(*,*)
        else
            rc=r_c_iter(spp,ex,ier)
            jc_xy=jc_dmless(rc,spp)
        end if
        !end block

        t0=r/jc_xy*(2*(10**logphi+spp%mbh_dmless/r-ex))**0.5
        !dlny=gxj_ir%ystep
        !print*, "t0,r,phi,dlny=",t0,r,10**logphi!,dlny
        !print*, "gj=", gxj_ir%ycenter
        !if(jmax_value/t0<1d0)then
        !    theta_max=asin(jmax_value/t0)
        !else
        !    theta_max=pi/2d0
        !end if
        !if(jmin_value/t0<1d0)then
        !    theta_min=asin(jmin_value/t0)
        !else
        !    theta_min=pi/2d0
        !end if
        !n=gxj_ir%ny
        !y(1:n)=gxj_ir%fxy(i,1:n)
        !call spline_mylib(gxj_ir%ycenter(1:n),y(1:n),n,yp1,ypn, y2(1:n))
        logt0=log10(t0)
        if(logt0>log10jmin_value)then
            if(logt0<log10jmax_value)then
                call return_idx(logt0, log10jmin_value,log10jmax_value,gxj_ir%ny,idx,0)
            else
                idx=gxj_ir%ny
            end if
            !print*, "gxj_ir%yb=",gxj_ir%ycenter
            !print*, "t0,log10(t0),idx=",t0,log10(t0),idx
            
            do j=1, idx-1
                jx0=10**(gxj_ir%ycenter(j)-gxj_ir%ystep/2d0)
                !jmid=10**gxj_ir%ycenter(j)
                jx1=10**(gxj_ir%ycenter(j)+gxj_ir%ystep/2d0)
                sintheta1=jx1/t0
                sintheta0=jx0/t0
                !theta1=asin(jx1/t0)
                !theta0=asin(jx0/t0)
                costheta1=(1-sintheta1*sintheta1)**0.5
                costheta0=(1-sintheta0*sintheta0)**0.5
                gxr%fx(i)=gxr%fx(i)+&
                gxj_ir%fxy(i,j)*(costheta0-costheta1)
                !gxj_ir%fxy(i,j)*(cos(theta0)-cos(theta1))

            end do
            jx0=10**(gxj_ir%ycenter(idx)-gxj_ir%ystep/2d0)
            !jmid=10**((gxj_ir%ycenter(idx)-gxj_ir%ystep/2d0+log10(t0))/2d0)
            !print*, "jx0,jx1, logt0=",jx0, jx1, log10(t0)
            !print*, "t0, log(t0), jx0,jx1,ystep=",t0, log10(t0),jx0,jx1, gxj_ir%ystep
            !theta1=pi/2d0
            !theta0=asin(jx0/t0)
            sintheta0=jx0/t0
            costheta0=(1-sintheta0*sintheta0)**0.5
            !print*, "theta0,theta1=",theta0,theta1
            gxr%fx(i)=gxr%fx(i)+gxj_ir%fxy(i,idx)*costheta0!jmid/t0*(theta1-theta0)
            !print*, "gxr=", gxj_ir%fxy(i,idx)*jx1*(theta1-theta0)
            !read(*,*)
            !do j=1, idx-1
            !    theta=(theta_max-theta_min)/dble(nbin_size)*(j-0.5)+theta_min
            !    jx=sin(theta)
            !    call return_idx(log10(jx*t0), log10jmin_value,log10jmax_value,gxj_ir%ny,idx,0)
            !    gxr%fx(i)=gxr%fx(i)+gxj_ir%fxy(i,idx)*jx*pi/2d0/dble(nbin_size)
                !if(idx<0)then
                !    print*, "jmin,jmax=",jmin_value,jmax_value
                !    print*, "theta_min,max=",theta_min,theta_max
                !    print*, "jx,t0=",jx,t0
                !    print*, "log10(jx*t0),idx=",log10(jx*t0), idx
                !    stop
                !end if
                !if(gxj_ir%ycenter(j)-gxj_ir%ystep/2d0<t0)then
                !    gxr%fx(i)=gxr%fx(i)+&
                !    gxj_ir%fxy(i,j)*jx*pi/2d0/dble(nbin_size)
                !end if
                !call splint_mylib(gxj_ir%ycenter(1:n),y(1:n),y2(1:n),n,log10(jx*t0), yout,ier)
                !call linear_int(gxj_ir%ycenter(1:n),y(1:n),n, log10(jx*t0), yout)
                !gxr%fx(i)=gxr%fx(i)+max(yout,0d0)*jx*pi/2d0/dble(nbin_size)
            !end do
            !print*, "ex,gx0, gx1=", ex,dm%mb(1)%star%barge_ir%fx(i), common_gx_ir%fx(i)
            !read(*,*)
        else
            if(logt0<=log10jmin_value)then
                gxr%fx(i)=0
            !else
            !    print*, "logt0??=",logt0, log10jmax_value
            !    stop
            end if
        end if
    end do
end subroutine
subroutine get_gxr(gxr,gxj_ir,logr,spp,dm,jc_pre)
    use com_main_gw
    implicit none
    type(s1d_ird_type)::gxr
    type(s2d_type)::gxj_ir
    type(diffuse_mspec)::dm
    type(star_pot_para)::spp
    logical::jc_pre
    real(8) r,ex,logphi,jx,theta
    real(8) logex, logr,jc_xy
    real(8) rc, r_c_iter,jc_dmless
    real(8) theta_max,theta_min
    integer ier, i,j,idx,n
    integer,parameter::nbin_size=300
    !real(8),parameter:: yp1=1d30, ypn=1d30
    !real(8) y2(gxj_ir%ny),yout,y(gxj_ir%ny), jx1, jmid,jx0,theta1,theta0,logt0
    real(8) jx1, jmid,jx0,theta1,theta0,logt0
    real(8) sintheta1,sintheta0,costheta1,costheta0
	external::func    
	real(8) exmax,exmin,fout,ratio,res, f0,f1
	real(8) tmin,tmax,fx_tmp,xmin,xmax,j0
	real(8) logx0, logx1,logxmid, logxminmid,logxmaxmid,t0,t1,phi,logxmin,logxmax
	real(8) line_slope, line_c,yp0,yp1,xp0,xp1,  res_tmp,dint, term_tmp
	integer debug, nbin,flag_end
	integer idx0,idx1,idxmin,idxmax
    
	nbin=ctl%dstr_bins_j
	tmin=0d0
	tmax=pi/2d0
    r=10**logr
    gxr%fx=0
    do i=1, ctl%dstr_bins_e
        ex=10**gxj_ir%xcenter(i)
        !print*, "i=",i, ex
        call get_phi_star_full_range(spp,logr,logphi)
        if(jc_pre)then
            !print*, size(dms%jc_sample_erange%fx)
            jc_xy=dms%jc_sample_erange%fx(i)
            
           ! print*, "1, xb, jc_xy=",10**dms%jc_sample_erange%xb(i),jc_xy
            !rc=r_c_iter(spp,ex,ier)
            !jc_xy=jc_dmless(rc,spp)
           ! print*, "2, ex, jc_xy=",ex, jc_xy
           !read(*,*)
        else
            rc=r_c_iter(spp,ex,ier)
            jc_xy=jc_dmless(rc,spp)
        end if
        !end block
        term_tmp=10**logphi+spp%mbh_dmless/r-ex
        if(term_tmp>0)then
            j0=r/jc_xy*(2*term_tmp)**0.5
            if(j0<=0)then
                print*, "r, jc_xy, term_tmp=",j0, r, jc_xy, term_tmp 
                stop
            end if
        else
            gxr%fx(i)=0
            cycle
        end if
        
        logxmin=log10jmin_value
        logxmax=min(log10jmax_value,log10(j0))
        idxmin=1
        if(logxmax<logxmin)then
            gxr%fx(i)=0
            cycle
        end if
        if(logxmax<log10jmax_value)then
            call return_idx(logxmax, log10jmin_value,log10jmax_value,gxj_ir%ny,idxmax,0)
        else
            idxmax=nbin
        end if

        if(idxmin<0.or.idxmax<0)then
            print*, "error!"
            print*, "logxmin,logxmax, logjmin,logjmax, nbin"
            print*, logxmin,logxmax, log10jmin_value,log10jmax_value, nbin
            print*, "idxmin,idxmax=",idxmin,idxmax
            stop
        end if
        logxminmid=gxj_ir%ycenter(idxmin)
        logxmaxmid=gxj_ir%ycenter(idxmax)
        !print*, "logr=",logr
        !print*, "logjmax, idxmax, logjmaxmid,j0=",logxmax, idxmax, logxmaxmid,j0

        
        logx0=logxmin
        res_tmp=0
        flag_end=0
100	    if(logx0<logxminmid)then 
            !if(ctl%debug.eq.1)then
            !	print*, "=logx0<logxminmid", logx0, logxminmid
            !end if
            if(idxmin.ne.1)then
                idx1=idxmin
                idx0=idxmin-1
                logxmid=gxj_ir%ycenter(idx1)
            else
                idx1=2
                idx0=1
                logxmid=gxj_ir%ycenter(1)
            end if
            if(logxmid<logxmax)then
                logx1=logxmid
            else
                logx1=logxmax
                flag_end=1
            end if
        else
            !if(ctl%debug.eq.1)then
            !	print*, "=logx0<logxminmid", logx0, logxminmid
            !end if
            if(idxmin.ne.nbin)then
                idx1=idxmin+1
                idx0=idxmin
                logxmid=gxj_ir%ycenter(idx1)
                if(logxmid<logxmax)then
                    logx1=logxmid
                else
                    logx1=logxmax
                    flag_end=1
                end if
            else
                idx1=nbin
                idx0=nbin-1
                logxmid=gxj_ir%ycenter(nbin)
                logx1=logxmax
                flag_end=1
            end if
            
        end if
        !if(ctl%debug.eq.1)then
        !	print*, "xmid,logxmaxmid",logxmid,logxmaxmid,logxmax
        !	print*, "x0,x1,flag_end=",logx0,logx1,flag_end
        !end if
        yp1=gxj_ir%fxy(i,idx1)
        yp0=gxj_ir%fxy(i,idx0)
        xp1=10**gxj_ir%ycenter(idx1)
        xp0=10**gxj_ir%ycenter(idx0)
        if(ctl%debug.eq.1)then
        	print*, "xp0,xp1,yp0,yp1=", xp0,xp1,yp0,yp1,idx0,idx1
        end if
        line_slope=(yp1-yp0)/(xp1-xp0)
        line_c=yp0-line_slope*xp0
        !print*, "l_slope,c=",line_slope,line_c
        !if(line_slope.eq.0)then
        !end if
        if(flag_end.eq.1)then
            if(10**logx1>-line_c/line_slope.and.line_slope.lt.0d0)then
                if(-line_c/line_slope>0)then
                    logx1=log10(-line_c/line_slope)
                else
                    print*, "error! -line_c/line_slope<0??"
                    print*, "2:10**logx1, jmin_value, jmax_value,line_slope,line_c=", &
                    10**logx1, jmin_value, jmax_value,line_slope,line_c, 10**logx1*line_slope+line_c
                    stop
                end if
            end if
        end if
        if(10**logx0<-line_c/line_slope.and.line_slope.gt.0d0)then
            if(-line_c/line_slope>0)then
                logx0=log10(-line_c/line_slope)
            else
                print*, "error! -line_c/line_slope>0??"
                print*, "2:10**logx1, jmin_value, jmax_value,line_slope,line_c=", &
                10**logx0, jmin_value, jmax_value,line_slope,line_c, 10**logx0*line_slope+line_c
                stop
            end if
        end if
        if(10**logx1<j0)then
            t1=asin(10**logx1/j0)
        else
            t1=pi/2d0    
        end if
        if(10**logx0<j0)then
            t0=asin(10**logx0/j0)    
        else
            t0=pi/2d0
        end if       
        
        
       ! print*, "t0,t1=",t0,t1
        f0=(t0/2d0-sin(2*t0)/4d0)*line_slope*j0-line_c*cos(t0)
        f1=(t1/2d0-sin(2*t1)/4d0)*line_slope*j0-line_c*cos(t1)
        dint=f1-f0
        res_tmp=res_tmp+dint
        !print*, "f0, f1=",f0, f1
        !read(*,*)
        !if(ctl%debug.eq.1)then
        !	print*, "res_tmp=",res_tmp
        !end if
        if(flag_end.ne.1)then
            logx0=logx1
            if(idxmin<nbin)then
                idxmin=idxmin+1
                logxminmid=gxj_ir%ycenter(idxmin)
            end if
            goto 100
        end if
        if(res_tmp.lt.0)	then
            print*, "res_tmp=",res_tmp, logx0, logx1, log10jmin_value, log10jmax_value
            print*, "f1,f0=",f1,f0, t0,t1,line_slope,j0,line_c
            ! read(*,*)
            gxr%fx(i)=1d-99
        else
            gxr%fx(i)=res_tmp
        end if
        !read(*,*)
    end do
	

end subroutine

subroutine get_barge_stellar_2d(so,jbtype)
    use md_dms
    implicit none
    type(dms_stellar_object)::so
    !real(8) sums
    integer i, j, idid, jbtype
    real(8) x,int_out
    if(so%n.eq.0)return
    do i=1, so%barge%nbin
        int_out=0
        !print*,"xmax=",mb%barge%xmax
        if(so%barge%nbin.ne.so%gxj%nx)then
            print*, "error! barge%nbin should = gxj%nx"
            stop
        end if
        so%barge%xb(i)=so%gxj%xcenter(i)
        select case(jbtype)
        case (Jbin_type_lin)
            do j=1, so%gxj%ny
                !==original==
                int_out=int_out+so%gxj%fxy(i,j)*so%gxj%ycenter(j)*so%gxj%ystep*2d0
                !===test=====
                !int_out=int_out+so%gxj%fxy(i,j)*so%gxj%ycenter(j)*so%gxj%ystep&
                !    /(1-so%gxj%ycenter(j)**2)**0.5d0
            end do
            so%barge%fx(i)=int_out
            !so%barge%nsam=so%n
            if(ieee_is_nan(so%barge%fx(i)))then
                print*, "get_barge_stellar:fx is NaN:", int_out
                call so%gxj%print()
                read(*,*)
            end if
        case(Jbin_type_log)
            do j=1, so%gxj%ny
                !==original==
                int_out=int_out+so%gxj%fxy(i,j)*(10**so%gxj%ycenter(j))**2*so%gxj%ystep*2d0*log(10d0)
                !===test=====
                !int_out=int_out+so%gxj%fxy(i,j)*so%gxj%ycenter(j)*so%gxj%ystep&
                !    /(1-so%gxj%ycenter(j)**2)**0.5d0
            end do
            so%barge%fx(i)=int_out
            !so%barge%nsam=so%n
            if(ieee_is_nan(so%barge%fx(i)))then
                print*, "get_barge_stellar:fx is NaN:", int_out
                call so%gxj%print()
                read(*,*)
            end if
        case(Jbin_type_sqr)
            do j=1, so%gxj%ny
                !==original==
                int_out=int_out+so%gxj%fxy(i,j)*so%gxj%ycenter(j)*so%gxj%ystep*2d0
                !===test=====
                !int_out=int_out+so%gxj%fxy(i,j)*so%gxj%ycenter(j)*so%gxj%ystep&
                !    /(1-so%gxj%ycenter(j)**2)**0.5d0
            end do
            so%barge%fx(i)=int_out
            !so%barge%nsam=so%n
            if(ieee_is_nan(so%barge%fx(i)))then
                print*, "get_barge_stellar:fx is NaN:", int_out
                call so%gxj%print()
                read(*,*)
            end if
        end select	
    end do
    call so%barge%prepare_spline()
end subroutine

subroutine get_barge_stellar_2d_ir(so,jbtype)
    use md_dms
    use MPI_comu,only:rid
    implicit none
    type(dms_stellar_object)::so
    !real(8) sums
    integer i, j, idid, jbtype
    real(8) x,int_out
    if(so%n.eq.0)return
    do i=1, so%barge_ir%nbin
        int_out=0
        !print*,"xmax=",mb%barge%xmax
        if(so%barge_ir%nbin.ne.so%gxj_ir%nx)then
            print*, "get_barge_stellar_2d_ir: error! barge%nbin should = gxj%nx"
            stop
        end if
        so%barge_ir%xb(i)=so%gxj_ir%xcenter(i)
        select case(jbtype)
        case (Jbin_type_lin)
            do j=1, so%gxj_ir%ny
                !==original==
                int_out=int_out+so%gxj_ir%fxy(i,j)*so%gxj_ir%ycenter(j)*so%gxj_ir%ystep*2d0
                !===test=====
                !int_out=int_out+so%gxj%fxy(i,j)*so%gxj%ycenter(j)*so%gxj%ystep&
                !    /(1-so%gxj%ycenter(j)**2)**0.5d0
            end do
            so%barge_ir%fx(i)=int_out
            !so%barge%nsam=so%n
            if(ieee_is_nan(so%barge_ir%fx(i)))then
                print*, "get_barge_stellar:Jbin_type_lin:fx is NaN:", int_out
                call so%gxj_ir%print()
                read(*,*)
            end if
        case(Jbin_type_log)
            !print*, "i=",i
            do j=1, so%gxj_ir%ny
                !==original==
                int_out=int_out+so%gxj_ir%fxy(i,j)*(10**so%gxj_ir%ycenter(j))**2*&
                    so%nxj_ir%ysteps(j)*2d0*log(10d0)
               ! print*, "so%gxj_ir%fxy(i,j)=",so%gxj_ir%fxy(i,j)
                !===test=====
                !int_out=int_out+so%gxj%fxy(i,j)*so%gxj%ycenter(j)*so%gxj%ystep&
                !    /(1-so%gxj%ycenter(j)**2)**0.5d0
            end do
            !print*, "int_out=",int_out 
            so%barge_ir%fx(i)=int_out
            !read(*,*)
            !so%barge%nsam=so%n
            if(ieee_is_nan(so%barge_ir%fx(i)))then
                print*, "get_barge_stellar:Jbin_type_log:fx is NaN:", int_out
                call so%gxj_ir%print()
                read(*,*)
            end if
        case(Jbin_type_sqr)
            do j=1, so%gxj%ny
                !==original==
                int_out=int_out+so%gxj%fxy(i,j)*so%gxj%ycenter(j)*so%gxj%ystep*2d0
                !===test=====
                !int_out=int_out+so%gxj%fxy(i,j)*so%gxj%ycenter(j)*so%gxj%ystep&
                !    /(1-so%gxj%ycenter(j)**2)**0.5d0
            end do
            so%barge%fx(i)=int_out
            !so%barge%nsam=so%n
            if(ieee_is_nan(so%barge%fx(i)))then
                print*, "get_barge_stellar:fx is NaN:", int_out
                call so%gxj%print()
                read(*,*)
            end if
        end select	
    end do
    !call so%barge_ir%print("barge_ir")
    !read(*,*)
    !call so%barge%prepare_spline()
end subroutine
subroutine dms_so_get_fxj_spt_ir(so, n0, pd,jc, r0,jbtype)
    use md_star_pot
    use md_stellar_object
    use model_basic,only:dms
    implicit none
    type(dms_stellar_object)::so
    integer i, j,jbtype
    real(8) jm ,x, n0,  r0
    type(s2d_type)::pd
    type(s1d_type) jc
    real(8) jc_xy,rp_xy,ra_xy,pd_xy

    if(so%n.eq.0) return
    select case(jbtype)
    case(Jbin_type_lin)

    case(Jbin_type_log)
        do i=1, so%nxj_ir%nx
            x=10**so%nxj_ir%xcenter(i)
            do j=1, so%nxj_ir%ny
                jm=10**so%nxj_ir%ycenter(j)
                !call jc%get_value_s(so%nxj%xcenter(i),jc_xy)
                !call pd%get_value_l(so%nxj%xcenter(i),so%nxj%ycenter(j),pd_xy)
                if(so%nxj_ir%nxyw(i,j).gt.0d0)then
                    call get_sample_para_one_xj_no_reset(x,jm,spp_new, dms%fr_phi,jc_xy,&
                        rp_xy,ra_xy,pd_xy)
                    !so%gxj%fxy(i,j)=so%nxj%nxyw(i,j)/(x*log(10d0))&
                    !/so%nxj%xstep/so%nxj%ystep &
                    !*pi**(-0.5d0)*2**(-1.5d0)/r0**3/(jm**2*log(10d0))/n0/jc%fx(i)**2/pd%fxy(i,j)
                    so%gxj_ir%fxy(i,j)=so%nxj_ir%nxyw(i,j)/(x*log(10d0))&
                    /so%nxj_ir%xsteps(i)/so%nxj_ir%ysteps(j) &
                    *pi**(-0.5d0)*2**(-1.5d0)/r0**3/(jm**2*log(10d0))/n0/jc_xy**2/pd_xy
                    !if(i.eq.so%nxj_ir%nx)then
                    !    print*, "gxj%fxy,nxyw,xstep,ystep,jm, jc_xy, pd_xy=", so%nxj_ir%nxyw(i,j), &
                    !        so%nxj_ir%xsteps(i), so%nxj_ir%ysteps(i), &
                    !        jm, jc_xy, pd_xy
                    !end if
                else
                    so%gxj_ir%fxy(i,j)=0
                end if
                !print*, "x,jm,nxyw=",x,jm,so%nxj_ir%nxyw(i,j)
                !print*, "xstep, ystep=",so%nxj_ir%xsteps(i), so%nxj_ir%ysteps(j)
                !print*, "r0,n0, pd_xy, jc_xy=", r0, n0, pd_xy, jc_xy
                !read(*, *)
                
            end do 
        end do
        !call so%gxj_ir%print("gxj_ir")
        !read(*,*)
    case(Jbin_type_sqr)

    case default
        print*, "fxj error!"
        stop
    end select
end subroutine