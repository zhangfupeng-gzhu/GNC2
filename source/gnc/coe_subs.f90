subroutine get_coegw(rp, ec, m, period, coegw)
	use com_main_gw
	implicit none
	real(8) rp, ec, m,period,delta_e,delta_l
	type(coeff_type)::coegw

	delta_e=-64d0/5d0*pi*(spp_new%mbh**2*m*(spp_new%mbh+m)**0.5)/(rp**3.5d0*my_unit_vel_c5)&
			/(1+ec)**(3.5d0)*(1+73d0/24d0*ec**2+37d0/96d0*ec**4)
	! delta_e=-64d0/5d0*pi*spp_new%mbh**2.5*m/(rp**3.5d0*my_unit_vel_c5)&
	! 		/(1+ec)**(3.5d0)*(1+73d0/24d0*ec**2+37d0/96d0*ec**4)
	coegw%e=delta_e/period
	
	delta_l=-64d0/5d0*pi*spp_new%mbh**2*m/my_unit_vel_c5&
			/rp**2d0/(1+ec)**2*(1+7d0/8d0*ec*ec)
	coegw%j=delta_l/period
end subroutine

subroutine get_coegw_ec(en0,jp0,jm0,m1,m2,dt, de,coegwdjexp,coegwdj)
	use com_main_gw
	implicit none
	type(binary)::by
	real(8) ac,ec,dt,de,coegwdj,af
	real(8) c0,m1,m2,eout, en0, enf, jp0,jm0,jf,ef,coegwdjexp

	! ac=by%a_bin;ec=by%e_bin
	! m1=by%ms%m; m2=by%mm%m
	! en0=-m2/2d0/ac
	! j0=(m2*ac*(1-ec**2))**0.5
	enf=en0+de
	af=-m2/2d0/enf
	ec=(1-jm0**2)**0.5
	ac=-m2/2d0/en0
	call get_c0(ac,ec,m1,m2,c0)
	call get_e_given_ca(c0,af,m1,m2,ef)
	jf=(m2*af*(1-ef**2))**0.5
	coegwdj=jf-jp0
	! if(ac>100)then
		! print*, "ac,ec,exp,now=",-m2/2/en0,ec,coegwdjexp,coegwdj
		! read(*,*)
	! end if

end subroutine

subroutine get_coegw_unbound(rp, ec, m,Period, coegw)
	use com_main_gw
	implicit none
	real(8) rp, ec, m,theta0,Period,delta_E, delta_L
	type(coeff_type)::coegw
	
	theta0=acos(1/ec)
	
	delta_E=-8d0/15d0*(spp_new%mbh**2*m*(spp_new%mbh+m)**0.5d0)/(rp**3.5*my_unit_vel_c5)&
			/(1+ec)**(3.5d0)*(24*(pi-theta0)*(1+73d0/24d0*ec**2+37d0/96d0*ec**4)+&
			(ec**2-1)**0.5*(301d0/6d0+673d0/12d0*ec**2))
	coegw%e=delta_E/period
	delta_l=-8d0/5d0*spp_new%mbh**2*m/my_unit_vel_c5&
			/rp**(2d0)/(1+ec)**2*((pi-theta0)*(8+7d0*ec*ec)+ec*sin(theta0)*(13+ec**2))
	coegw%j=delta_L/period
end subroutine
    
subroutine get_step_gw(coegw, en, jph, jm, period)
	use com_main_gw
	implicit none
	type(coeff_type)::coegw
	real(8) en, jph, period, jm,nmingwe, nmingwj,jgwr

	nmingwe=0.01d0*abs(en/(coegw%e*period))!*abs(1.1-jm)
	jgwr=abs(coegw%j)*period
	nmingwj=0.005d0*jph/jgwr!*abs(1.1-jm)
	if(jm>0.1d0)then
		! steps=min(steps, nmingwe, nmingwj)
		sample_step_gw=min(nmingwe, nmingwj)
	else
		sample_step_gw=max(min(nmingwe, nmingwj),0.2d0)
	end if
	
	if(ctl%chattery.ge.4)then
		if(ctl%ntasks.gt.1)then 
		else
			print*, "coegw%e, coegw%j=",coegw%e, coegw%j                    
			print*, "sample_step_gw, nmingwe, nmingwj=", sample_step_gw, nmingwe, nmingwj	
		end if
		!read(*,*)	
	end if
end subroutine


subroutine get_de_dj(sample,coeNr, coeRR, coeGW,  time, dt,steps, period)        !use com_main_gw
        use com_main_gw
        implicit none
        type(particle_sample_type)::sample
        real(8),intent(in) ::dt, steps,time
        real(8) den, gen_gaussian
        real(8) y1,y2,y3,y4,rho, deb2, jc, jum, n2j
        real(8) dpe,ai, ei, Eni,Enf, Ji,Jf, af, ef
        integer jb, ju
        real(8) ipdi, ipdf
        real(8) period, npi, npf, coegwde,coegwdj
        !real(8),save::t0=0
        type(coeff_type)::coeNr, coeRR, coeGW
        
        if(ieee_is_nan(dt).or.(.not.ieee_is_finite(dt)).or.(.not.ieee_is_finite(steps)))then
            print*, "dt=!, dt, steps=", dt, steps
            stop
        end if

        if(ieee_is_nan(sample%jm).or.ieee_is_nan(sample%en))then
            print*
            print*, "ac,ec=",sample%byot%a_bin,sample%byot%e_bin
            print*, "en,jm=",sample%en, sample%jm,dt
            read(*,*)
        end if 

		call get_de_dj_nr(coenr, dt, steps, sample%jph, sample_den, sample_djp, &
			sample_djp0)
        
        if(ieee_is_nan(sample_djp))then
            print*, "djp=",sample_djp, coeNR%jj, dt, coeNR%j
        end if
        !print*, "den, djp, djp0=", sample_den, sample_djp, sample_djp0
        !stop
        if(ctl%gw_radiation_otby.ge.1 )then
            if(steps.ge.1.or.sample%jm>0.4d0)then
                sample_den=sample_den+coeGW%e*dt
                if(sample%state_emri_current.ge.1.and.sample%byot%a_bin>0d0.and.sample_alpha<ctl%md_alpha_cri)then
                    call get_coegw_ec(sample%en, sample%jph,sample%jm,sample%m, spp_new%mbh,&
                        dt,coeGW%e*dt, coegw%j*dt,coegwdj)
                    sample_djp=sample_djp+coegwdj
                    sample_djp0=sample_djp0+coegwdj/steps
                else
                    sample_djp=sample_djp+coeGW%j*dt
                    sample_djp0=sample_djp0+coeGW%j*dt/steps
                end if


                if(ctl%chattery.ge.4)then
                    print*, "========get_dedj:coegw====================="
                    print*, "steps,ec=",steps,sample%byot%e_bin
                    print*, "delta_egw,delta_jgw=",coeGW%e*dt,coeGW%j*dt
                end if
            else
                !ipdx=(time-t0)/period-steps
                !ipd=(time-t0)/period
				npi=sample%byot%me
				npf=npi+steps*pi*2
				ipdi=npi/2d0/pi-int(npi/2d0/pi)
				ipdf=npf/2d0/pi-int(npf/2d0/pi)
                !print*, "GW:ipdi,ipdf=", ipdi, ipdf
				!print*, sample%byot%me, sample%byot%e_bin
				if(ipdi<0.5.and.ipdf>0.5)then
                    sample_den=sample_den+coeGW%e*period
                    sample_djp=sample_djp+coeGW%j*period
					!if(sample%byot%a_bin<10)then
					!	print*, "coegw%e, period="
					!end if
					!print*, "GW decayed"
                    !if(sample%byot%a_bin<1d3)then
                    !    print*, "en, den, gwe,gwedt=", sample%en, sample_den, coeGW%e, coeGW%e*dt
                    !    print*, "jm, djp, gwj,gwjdt=", sample%jm, sample_djp, coeGW%j, coeGW%j*dt
                    !    read(*,*)
                    !end if
                    sample_djp0=sample_djp0+coeGW%j
             !       t0=time
                end if
                if(ctl%chattery.ge.4)then
                    print*, "========get_dedj:coegw====================="
                    print*, "steps,ipdi,ipdf=",steps,ipdi,ipdf
                    print*, "delta_egw,delta_jgw=",coeGW%e*period,coeGW%j*period
                end if
            end if

        end if		

        !print*, "sample%jm=", sample%jm

    
    end subroutine


subroutine move_de_dj_one(spp,sample, eni,enf, jf, mef,af)
	use com_main_gw
	implicit none
	type(particle_sample_type)::sample
	real(8) den, djp,steps
	real(8) enf, jf, mef,af
	real(8) eni, ai, ei,jmfi, even
	integer idx
	type(star_pot_para)::spp 
	eni=sample%en
	sample%byot_bf%e=sample%en
	sample%byot_bf%l=sample%jm
	sample%byot_bf%Jc=sample%jc
	sample%en=enf
	sample%jph=jf
	!sample%x=enf/ctl%energy0
	
	!sample%jm=jmf
	!sample%byot%e_bin=sqrt(1-sample%jm**2)
	if(spp%mbh_dmless.eq.0)then
		call set_ebound_for_samples(sample)
	end if
	sample%byot_bf%me=sample%byot%me
	if(mef>2*pi)then
		sample%byot%me=mef-2*pi
	else
		sample%byot%me=mef
	end if
	
	if(ieee_is_nan(sample%jm).and.sample%en<0d0)then
		print*, "error! sample%jm=NaN", sample%jm
		print*, "sample%m,id=",sample%m,sample%id

		print*, "enmin,enmax=",ctl%energy_min,ctl%energy_max
		print*, "jm, en=", sample%jm, sample%en, sample_den
		print*, "sample%ngene=",sample%n_gene
		print*, "sample%byot%a_bin,djp=",sample%byot%a_bin,sample_djp
		stop
	end if

	if(ctl%chattery.ge.4)then
		print*, "======start de_dj==========================="
		print*, "befor: sample%en,jm=",sample%byot_bf%e/ctl%energy0,sample%byot_bf%l
		print*, "after: sample%en,jm(estimate)=",sample%en/ctl%energy0, sample%jph/sample%jc
	end if
	end subroutine
	 