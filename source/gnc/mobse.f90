
  	subroutine mbos_set_model_init()
		implicit none
		include 'const_mobse.h'
		neta=0.5
        !neta is the Reimers mass-loss coefficent (neta*4x10^-13: 0.5 normally). 
        bwind=0
        !bwind is the binary enhanced mass loss parameter (inactive for single).
        hewind=1.0
        !hewind is a helium star mass loss factor (normally inactive).
        alpha1=3.0
        !alpha1 is the common-envelope efficiency parameter (1.0).  
        lambda=0.1
        !lambda is the binding energy factor for common envelope evolution (0.1).
		ceflag=0
        !ceflag > 0 activates spin-energy correction in common-envelope (0). #defunct#
        tflag=1
        !tflag > 0 activates tidal circularisation (1).
        ifflag=0
        !ifflag > 0 uses WD IFMR of HPE, 1995, MNRAS, 272, 800 (0). 
        wdflag=1
        !wdflag > 0 uses modified-Mestel cooling for WDs (1). 
        bhflag=3
        !bhflag > 0 allows velocity kick at BH formation (3). 
        !     1 --> old case
        !     2 --> fallback case: Fryer et al. 2012, ApJ, 749, 91
        !     3 --> mass ejected case: Giacobbo & Mapelli 2020, ApJ, 891
        !     4 --> full kick case

        nsflag=2
        !     nsflag > 0 takes NS/BH mass from (default 3): 
        !     1 --> Belczynski et al. 2002, ApJ, 572, 407
        !     2 --> Rapid supernova model Fryer et al. 2012, ApJ, 749, 91
        !     3 --> Delayed supernova model Fryer et al. 2012, ApJ, 749, 91
        piflag=1
        !piflag > 0 activates the PPISNs and PISNe (1) Spera et al. 2017, MNRAS, 470, 4739.  

		mxns=3.0; 
        !mxns is the maximum NS mass (3.0).
		pts1=0.05; pts2=0.01; pts3=0.02
        ! Next come the parameters that determine the timesteps chosen in each
        ! evolution phase:
        !                 pts1 - MS                  (0.05) 
        !                 pts2 - GB, CHeB, AGB, HeGB (0.01)
        !                 pts3 - HG, HeMS            (0.02)

		sigma1=265d0; sigma2=15d0
        ! sigma1 is the dispersion in the Maxwellian for the SN kick speed (265. km/s)
        !           from Hobbs et al. 2005, ApJ, 591, 288.
        ! sigma2 is the dispersion in the Maxwellian for the SN kick speed (15. km/s)
        !           to considere the different mechanims involve in ECS. 
		beta=0.125
        !beta is wind velocity factor: proportional to vwind**2 (1/8). 
        xi=1.0
        !xi is the wind accretion efficiency factor (1.0). 
		acc2=1.5
        !acc2 is the Bondi-Hoyle wind accretion factor (3/2). 
		epsnov=0.001
        !epsnov is the fraction of accreted matter retained in nova eruption (0.001). 
		eddfac=1.0
        !eddfac is Eddington limit factor for mass transfer (1.0).
		gamma=-1.0
        !gamma is the angular momentum factor for mass lost during Roche (-1.0). 
	end subroutine

	

	subroutine get_evl_single_final_mobse(mass0_in,  tmax_evl_myr, kstar_in, z, &
		random_seed,kstar_f, mass_f, kw_f, radius_f, nrecord_max, output_flag)
		implicit none
		integer nrecord_max
		real(8) mass0_in, tmax_evl_myr,z
		integer random_seed, output_flag, kstar_in, kstar_f
		real(8) time_seq(nrecord_max)
		real(8) mass_seq(nrecord_max),mass_f, radius_seq(nrecord_max), radius_f
		integer kstar_seq(nrecord_max), nseq,  kw_seq(nrecord_max), kw_f
		
		call  get_evl_single_mobse(mass0_in, tmax_evl_myr, kstar_in, z, &
		random_seed,time_seq, kstar_seq, mass_seq, kw_seq, radius_seq, nseq, nrecord_max, .false., output_flag)
		kstar_f=kstar_seq(nseq)
		mass_f=mass_seq(nseq)
		kw_f=kw_seq(nseq)
		radius_f=radius_seq(nseq)
	end subroutine

	subroutine get_evl_single_mobse(mass0_in, tmax_evl_myr, kstar_in, z, &
				random_seed,time_seq, kstar_seq, mass_seq, kw_seq, radius_seq, nseq,nrecord_max, save_details, output_flag)
		use md_mobse_stellar_single,only:get_kstar_type,get_kw_type
		implicit none
		include 'const_mobse.h'
		integer nrecord_max
		real(8) mass0(2), mass0_in, tmax_evl_myr, t_orb_day, z, ecc, tphys,dtp, tphysf
		integer kstar(2), random_seed, j, kw,kw2, output_flag, kstar_in
		real(8) rad(2),radc(2),menv(2),renv(2),massc(2), lum(2), ospin(2)
		real(8) epoch(2), mass(2), tms(2),zpars(20)
		real(8) time_seq(nrecord_max)
		real(8) mass_seq(nrecord_max),radius_seq(nrecord_max)
		integer kstar_seq(nrecord_max), nseq,  kw_seq(nrecord_max)	
		logical::save_details
		character*(8) label(14)
		!character*(5) get_kstar_type
		!character*(9) get_kw_type
		character*(9) kwtype,kstartype
		integer,parameter::KSTAR_ZAMS=1, KSTAR_EVLOVE=-1

		kstar=(/kstar_in,1/)
		mass0=(/mass0_in,0d0/)
		t_orb_day=0d0; ecc=0d0
		!call get_evl_binary(mass0_input, tmax_evl_myr, t_orb_day, kstar_input, &
		!z, ecc, random_seed, output_flag)

		if(z>0.03.or.z<0.0001) then 
			print*, "z should be between 0.0001<z<0.03", z
			stop
		end if
		!print*, "start",mass0, tmax_evl_myr, t_orb_day, kstar, &
		!z, ecc, random_seed, output_flag
		call mbos_set_model_init()
		idum=random_seed
        !idum is the random number seed used by the kick routine. 
		if(kstar(1).eq.KSTAR_EVLOVE)then
			print*, "not yet adopted! finished the code when needed!"
			stop
			!READ(22,*)tphys
			!READ(22,*)aj,mass(1),ospin(1)
			!epoch(1) = tphys - aj
			!kstar(1) = ABS(kstar(1))
			!READ(22,*)aj,mass(2),ospin(2)
			!epoch(2) = tphys - aj
			!kstar(2) = ABS(kstar(2))
		else
			tphys=0d0
			mass(1)=mass0(1)
			mass(2)=0d0
			epoch(1)=0d0
			epoch(2)=0d0
			ospin(1)=0d0
			ospin(2)=0d0
		end if
		if(idum.gt.0) idum=-idum
		tphysf=tmax_evl_myr
		!Set parameters which depend on the metallicity 
		!print*, "1"
		call zcnsts(z, zpars)
		!Set the collision matrix.
		!print*, "2"
		call instar
		!print*, "2.5"
		! Set the array with the labels.
		!


	! Set the data-save parameter. If dtp is zero then the parameters of the 
	! star will be stored in the bcm array at each timestep otherwise they 
	! will be stored at intervals of dtp. Setting dtp equal to tphysf will 
	! store data only at the start and end while a value of dtp greater than 
	! tphysf will mean that no data is stored.	
		dtp = 0.d0
	! Evolve the binary.
	
		CALL evolve(kstar,mass0,mass,rad,lum,massc,radc, &
					menv,renv,ospin,epoch,tms, &
					tphys,tphysf,dtp,z,zpars,t_orb_day,ecc)
		nseq=0
		if(save_details)then
			
			do while(bcm(nseq+1,1).ge.0)
				nseq=nseq+1
				if(nseq>nrecord_max) then
					print*, "nseq>nrecord_max expand the arr!"
					print*, "mass0_in=",mass0_in
					stop
				end if
				kstar_seq(nseq) = INT(bcm(nseq,2))
				j=1
				do while (bpp(j,1).ge.0)
					if(bpp(j,1).eq. bcm(nseq,1))then
						exit 
					end if
					j=j+1
				end do
				!print*,"j=",j
				kw_seq(nseq) = INT(bpp(j,33))
				time_seq(nseq)=bcm(nseq,1)
				mass_seq(nseq)=bcm(nseq,4)
				radius_seq(nseq)=bcm(nseq,6)
				!print*, "radius_seq=",radius_seq(nseq), bpp(nseq, 1:7)
			end do	
			nseq=nseq-1
		else
			do while(bpp(nseq+1,1).ge.0)
				nseq=nseq+1
				if(nseq>100) then
					print*, "nseq>100 expand the arr!"
					stop
				end if
				kstar_seq(nseq) = INT(bpp(nseq,2))
				kw_seq(nseq) = INT(bpp(nseq,33))
				time_seq(nseq)=bpp(nseq,1)
				mass_seq(nseq)=bpp(nseq,4)
				radius_seq(nseq)=bpp(nseq,6)
				!print*, "radius_seq=",radius_seq(nseq), bpp(nseq, 1:7)
			end do	
		end if
	!***********************************************************************
	! Output:
	! First check that bcm is not empty.
	!
		!print*, "3"
		if(output_flag.ge.1)then
			if(bcm(1,1).lt.0.0) goto 50
	! The bcm array stores the stellar and orbital parameters at the 
	! specified output times. The parameters are (in order of storage in bcm):
	!
	!    1-Time, 
	!    2-15[stellar type, initial mass, current mass, log10(L), log10(r),
	!    log10(T), core mass, core radius, mass of any convective 
	!    envelope, radius of the envelope, epoch, spin, mass loss rate and 
	!    ratio of radius to roche lobe radius],
	!    16-29(repeated for secondary)
	!    30-32[period[year], separation(in unit of solar radius), eccentricity].
	!
	!	OPEN(24,file='../mobse.out',status='unknown')
			j = 0
			WRITE(*,98)"Time", "KW", "mass", "mcore", "log10(L)", "log10(r)", &
				"log10(T)", "TYPE"
	30	   j = j + 1
			kw = INT(bcm(j,2))
			kw2 = INT(bcm(j,16))
			kstartype=adjustr(get_kstar_type(int(bcm(j,2))))
			WRITE(*,99)bcm(j,1),kw,bcm(j,4), &
				bcm(j,8),bcm(j,5),bcm(j,6),bcm(j,7), kstartype
			if(bcm(j+1,1).ge.0.0) goto 30
	!	CLOSE(24)
	98	   FORMAT(A10, A3, 5A10, A8)	
	99	   FORMAT(f10.4,i3,5f10.4,A8)
	!
	! The bpp array acts as a log, storing parameters at each change
	! of evolution stage (it has the same information of bcm plus the 
	! labels at the end).
	!
	50   j = 0
			WRITE(*, 101)'TIME',"M","R","K" , 'TYPE'
	52   j = j + 1
			if(bpp(j,1).lt.0.0) goto 60
				kstar(1) = INT(bpp(j,2))
				kstar(2) = INT(bpp(j,16))
				kw = INT(bpp(j,33))
				kwtype=adjustr(get_kw_type(kw))
				kstartype=adjustr(get_kstar_type(kstar(1)))
				WRITE(*,100)bpp(j,1),bpp(j,4), bpp(j,6), kstartype,kwtype 
			goto 52
	60    continue
		read(*,*)
	101    FORMAT(a11,a11,a11,a5,a10)
	100   FORMAT(f11.4,2f11.5,a11,a5,a10)
		end if

	end subroutine
	subroutine get_evl_binary_mobse(mass0, tmax_evl_myr, t_orb_day, kstar, &
			z, ecc, random_seed,time_seq, kstar_seq, mass_seq, pd_seq, r_seq, &
			ec_seq, kw_seq, nseq, nrecord_max, output_flag)
		use md_mobse_stellar_single,only:get_kstar_type
		implicit none
		include 'const_mobse.h'
		integer nrecord_max
		real(8) mass0(2), tmax_evl_myr, t_orb_day, z, ecc, tphys,dtp, tphysf
		integer kstar(2), random_seed, j, kw,kw2, output_flag
		real(8) rad(2),radc(2),menv(2),renv(2),massc(2), lum(2), ospin(2)
		real(8) epoch(2), mass(2), tms(2),zpars(20)
		real(8) time_seq(nrecord_max), r_seq(nrecord_max)
		real(8) mass_seq(2, nrecord_max), pd_seq(nrecord_max), ec_seq(nrecord_max)
		integer kstar_seq(2, nrecord_max), nseq,  kw_seq(nrecord_max)
		character*(8) label(14)
		!character*(5) get_kstar_type
		integer,parameter::KSTAR_ZAMS=1, KSTAR_EVLOVE=-1
		! Note that this routine can be used to evolve a single star if you 
		! simply set mass0(2) = 0.0 or t_orb_day = 0.0 (setting both is advised as  
		! well as some dummy value for ecc). 

		if(z>0.03.or.z<0.0001) then 
			print*, "z should be between 0.0001<z<0.03", z
			stop
		end if
		!print*, "start",mass0, tmax_evl_myr, t_orb_day, kstar, &
		!z, ecc, random_seed, output_flag
		call mbos_set_model_init()
		idum=random_seed
		if(kstar(1).eq.KSTAR_EVLOVE.or.kstar(2).eq.KSTAR_EVLOVE)then
			print*, "not yet adopted! finished the code when needed!"
			stop
			!READ(22,*)tphys
			!READ(22,*)aj,mass(1),ospin(1)
			!epoch(1) = tphys - aj
			!kstar(1) = ABS(kstar(1))
			!READ(22,*)aj,mass(2),ospin(2)
			!epoch(2) = tphys - aj
			!kstar(2) = ABS(kstar(2))
		else
			tphys=0d0
			mass(1)=mass0(1)
			mass(2)=mass0(2)
			epoch(1)=0d0
			epoch(2)=0d0
			ospin(1)=0d0
			ospin(2)=0d0
		end if
		if(idum.gt.0) idum=-idum
		tphysf=tmax_evl_myr
		!Set parameters which depend on the metallicity 
		!print*, "1"
		call zcnsts(z, zpars)
		!Set the collision matrix.
		!print*, "2"
		call instar
		!print*, "2.5"
		! Set the array with the labels.
		!
		label(1) = 'INITIAL '
		label(2) = 'KW_CHNGE'
		label(3) = 'BEG_RCHE'
		label(4) = 'END_RCHE'
		label(5) = 'CONTACT '
		label(6) = 'COELESCE'
		label(7) = 'COMENV  '
		label(8) = 'GNTAGE  '
		label(9) = 'NO_REMNT'
		label(10) = 'MAX_TIME'
		label(11) = 'DISRUPT '
		label(12) = 'BEG_SYMB'
		label(13) = 'END_SYMB'
		label(14) = 'BEG_BSS'

	! Set the data-save parameter. If dtp is zero then the parameters of the 
	! star will be stored in the bcm array at each timestep otherwise they 
	! will be stored at intervals of dtp. Setting dtp equal to tphysf will 
	! store data only at the start and end while a value of dtp greater than 
	! tphysf will mean that no data is stored.	
		dtp = 0.d0
	! Evolve the binary.
	
		CALL evolve(kstar,mass0,mass,rad,lum,massc,radc, &
					menv,renv,ospin,epoch,tms, &
					tphys,tphysf,dtp,z,zpars,t_orb_day,ecc)

	!***********************************************************************
	! Output:
	! First check that bcm is not empty.
	!
		nseq=0
		do while(bpp(nseq+1,1).ge.0)
			nseq=nseq+1
			if(nseq>nrecord_max) then
				print*, "nseq>nrecord_max expand the arr!"
				stop
			end if
			kstar_seq(1, nseq) = INT(bpp(nseq,2))
			kstar_seq(2, nseq) = INT(bpp(nseq,16))
			kw_seq(nseq) = INT(bpp(nseq,33))
			time_seq(nseq)=bpp(nseq,1)
			mass_seq(1,nseq)=bpp(nseq,4)
			mass_seq(2,nseq)=bpp(nseq,18)
			
			r_seq(nseq)=bpp(nseq,31)
			ec_seq(nseq)=bpp(nseq,32)
			pd_seq(nseq)=((r_seq(nseq)*0.00465)**3/(mass_seq(1,nseq)+mass_seq(2,nseq)))**0.5&
				*365.2425
		end do	
		!if(time_seq(nseq).ne.tmax_evl_myr)then
		!	output_flag=1
		!end if
		!print*, "3"
		if(output_flag.ge.1)then
			if(bcm(1,1).lt.0.0) goto 50
	! The bcm array stores the stellar and orbital parameters at the 
	! specified output times. The parameters are (in order of storage in bcm):
	!
	!    1-Time, 
	!    2-15[stellar type, initial mass, current mass, log10(L), log10(r),
	!    log10(T), core mass, core radius, mass of any convective 
	!    envelope, radius of the envelope, epoch, spin, mass loss rate and 
	!    ratio of radius to roche lobe radius],
	!    16-29(repeated for secondary)
	!    30-32[period[year], separation(in unit of solar radius), eccentricity].
	!
	!	OPEN(24,file='../mobse.out',status='unknown')
			j = 0
			WRITE(*,98)"Time", "KW1", "mass1", "mcore1", "log10(L)", "log10(r)", &
				"log10(T)","KW2", "mass2", "mcore2", "log10(L)", "log10(r)", &
				"log10(T)", "period[year]", "r", "ecc"
	30	   j = j + 1
			kw = INT(bcm(j,2))
			kw2 = INT(bcm(j,16))
			
			WRITE(*,99)bcm(j,1),kw,bcm(j,4), &
				bcm(j,8),bcm(j,5),bcm(j,6),bcm(j,7), &
				kw2,bcm(j,18),bcm(j,22), &
				bcm(j,19),bcm(j,20),bcm(j,21), &
				bcm(j,30),bcm(j,31),bcm(j,32)
			if(bcm(j+1,1).ge.0.0) goto 30
	!	CLOSE(24)
	98	   FORMAT(A10, A3, 5A10, A3, 5A10, 2A16, A7, A8)	
	99	   FORMAT(f10.4,i3,5f10.4,i3,5f10.4,2f16.4,f7.3,a8)
	!
	! The bpp array acts as a log, storing parameters at each change
	! of evolution stage (it has the same information of bcm plus the 
	! labels at the end).
	!
	50   j = 0
			WRITE(*, 101)'TIME',"M1","M2","K1", "K2", "SEP", 'ECC', &
					'R1/ROL1', 'R2/ROL2',  'TYPE'
	52   j = j + 1
			if(bpp(j,1).lt.0.0) goto 60
				kstar(1) = INT(bpp(j,2))
				kstar(2) = INT(bpp(j,16))
				kw = INT(bpp(j,33))
				WRITE(*,100)bpp(j,1),bpp(j,4),bpp(j,18),get_kstar_type(kstar(1)),&
				get_kstar_type(kstar(2)), &
				bpp(j,31), bpp(j,32),bpp(j,15),bpp(j,29),label(kw)

			goto 52
	60    continue
	101    FORMAT(a11,2a11,2a5,a15,a6,2a8,2x,a8)
	100   FORMAT(f11.4,2f11.5,2a5,f15.5,f6.2,2f8.3,2x,a8)
		end if
	!
	!***********************************************************************

	end subroutine
	
