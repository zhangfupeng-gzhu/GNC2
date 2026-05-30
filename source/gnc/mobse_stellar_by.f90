
! module md_mobse_stellar_by
! 	use md_mobse_stellar_single
! 	implicit none
! 	type, extends(stellar_object_arr):: stellar_object_arr_by
! 		real(8),allocatable:: mass2(:), pd(:),ec(:),r(:),ac(:)
! 		real(8),allocatable::mass_tot(:)
! 		integer,allocatable::ktype2(:)
! 	contains
! 		procedure::init=>init_stellar_object2
! 		procedure::output_txt=>output_txt_stellar_object
! 		procedure::write_info=>write_info_stellar_object_by
! 		procedure::read_info=>read_info_stellar_object_by
! 	end type
! 	type, extends(stellar_object_arr_by)::stellar_object_hs_by
!         integer,allocatable:: kwtype(:)
!         contains
!         procedure::init=>init_stellar_object_hs2
!         procedure::write_info=>write_info_stellar_object_hs_by
!         procedure::print=>print_stellar_history2
! 		procedure::read_info=>read_info_stellar_object_hs_by
!     end type
! 	!(stellar_object) :: mobse_star, mobse_bh
! 	private::init_stellar_object_hs2
! 	private::init_stellar_object2
! 	private::output_txt_stellar_object
! 	private::read_info_stellar_object_by
! 	private::write_info_stellar_object_by
!     private::print_stellar_history2
! 	private::write_info_stellar_object_hs_by, read_info_stellar_object_hs_by
! contains
! 	subroutine write_info_stellar_object_by(this, funit)
!         implicit none
!         class(stellar_object_arr_by)::this
!         integer funit
!         call this%stellar_object_arr%write_info(funit)
! 		if(this%n<1) return
!         write(unit=funit) this%ktype2(1:this%n), this%mass2(1:this%n), &
!         this%ac(1:this%n), &
!         this%ec(1:this%n), this%pd(1:this%n),this%mass_tot(1:this%n), this%r(1:this%n)
!     end subroutine
! 	subroutine read_info_stellar_object_by(this, funit)
!         implicit none
!         class(stellar_object_arr_by)::this
!         integer funit
!         call this%stellar_object_arr%read_info(funit)
!         if(this%n<1) return
!         allocate(this%ktype2(this%n), this%mass2(this%n), &
!           this%ac(this%n), this%ec(this%n), this%pd(this%n), &
!           this%mass_tot(this%n),this%r(this%n))

!         read(unit=funit) this%ktype2(1:this%n), this%mass2(1:this%n), this%ac(1:this%n), &
!         this%ec(1:this%n), this%pd(1:this%n),this%mass_tot(1:this%n), this%r(1:this%n)
!     end subroutine
! 	subroutine write_info_stellar_object_hs_by(this, funit)
!         implicit none
!         class(stellar_object_hs_by)::this
!         integer funit
!         call write_info_stellar_object_by(this, funit)
! 		if(this%n<1)return
!         write(unit=funit) this%kwtype(1:this%n)
!     end subroutine
!     subroutine read_info_stellar_object_hs_by(this, funit)
!         implicit none
!         class(stellar_object_hs_by)::this
!         integer funit
!         call read_info_stellar_object_by(this, funit)
! 		if(this%n<1)return
!         if(.not.allocated(this%kwtype))allocate(this%kwtype(1:this%n))
!         read(unit=funit) this%kwtype(1:this%n)
!     end subroutine
    
!     subroutine init_stellar_object_hs2(this, n)
!         implicit none
!         class(stellar_object_hs_by)::this
!         integer n
!         call init_stellar_object2(this, n)
!         if(allocated(this%kwtype)) then
!             deallocate(this%kwtype)
!         end if
!         allocate(this%kwtype(n))
!     end subroutine
! 	subroutine init_stellar_object2(this,n)
!         implicit none
!         class(stellar_object_arr_by)::this
!         integer n,ktype, ktype2
!         call this%stellar_object_arr%init(n)
!         if(allocated(this%mass2))then
!             deallocate(this%mass2, this%mass_tot, this%pd,this%ec,this%ac, &
!                 this%ktype2, this%r)
!         end if
!         allocate(this%mass2(n), this%mass_tot(n), this%pd(n),this%ec(n),&
!             this%ac(n),this%ktype2(n), this%r(n))
!     end subroutine
! 	subroutine print_stellar_history2(this)
!         class(stellar_object_hs_by)::this
!         integer i
!         write(*,fmt="(10A15)") "time", "m1", "m2", "pd", "e", "type1", "type2", "event"
!         do i=1, this%n
!             write(*,fmt="(3F15.3, 2F15.6, 2A15, A15)") this%age(i), this%mass(i),this%mass2(i), &
!             this%pd(i),this%ec(i), get_kstar_type(this%ktype(i)), get_kstar_type(this%ktype2(i)), &
!                 get_kw_type(this%kwtype(i) )
!         end do
!     end subroutine
 
!     subroutine output_txt_stellar_object(this,fl)
!         implicit none
!         class(stellar_object_arr_by)::this
!         character*(*) fl
!         integer i
!         open(unit=19999,file=trim(adjustl(fl))//".txt")
!         write(unit=19999,fmt="(11A20)") "mtot", "m1", "m2", "ac(AU)", "ec", "pd(days)", "w"
!         do i=1, this%n
!             write(unit=19999,fmt="(11E20.10)") this%mass_tot(i), this%mass(i), this%mass2(i), &
!                 this%ac(i), this%ec(i), this%pd(i), this%w(i)
!         end do
!         close(unit=19999)
!     end subroutine
! 	subroutine get_snapshot_of_samples_by(sh,n, curtime, bstar, bbh,sbh)
!         implicit none
!         integer i, j, n
!         type(stellar_object_hs_by)::sh(n)
!         type(stellar_object_arr_by)::bstar, bbh
!         type(stellar_object_arr)::sbh
!         real(8) curtime		
!         integer nbstar, nbbh, nsbh
!         integer,allocatable::idx(:)

!         allocate(idx(n))

!         !star%ktype=mobse_ktype_ms
!         !bh%ktype=mobse_ktype_bh
!         nbstar=0; nbbh=0; nsbh=0
!         if(n.eq.0)then 
!             print*, "warnning: stellar_history_object_n=0"
!             return
!         end if
!         do i=1, n
!             if(curtime<sh(i)%age(1).or.curtime>sh(i)%age(sh(i)%n))then
!                 !print*, "i=",i,"curtime not in rage", curtime, sh(i)%age(1), sh(i)%age(sh(i)%n)
!                 cycle
!             end if
            
! loop1:		do j=1, sh(i)%n-1
!                 if((sh(i)%age(j).le.curtime).and.(sh(i)%age(j+1)>curtime))then
!                     idx(i)=j
!                     exit loop1
!                 end if
!             end do loop1
!             !call sh(i)%print()
!             !print*, "ctime,idx=",curtime,idx(i)
!             if(sh(i)%age(sh(i)%n).eq.curtime)then
!                 idx(i)=sh(i)%n
!                 !goto 100
!             end if
!             associate(ktype1=>sh(i)%ktype(idx(i)), ktype2=>sh(i)%ktype2(idx(i)))
!                 select case(ktype1)
!                 case(mobse_ktype_ms, mobse_ktype_fcms)
!                     select case(ktype2)
!                     case(mobse_ktype_ms, mobse_ktype_fcms)
!                         if(sh(i)%pd(idx(i)).ne.0d0.or.abs(sh(i)%ec(idx(i)))<1d0)then
!                             nbstar=nbstar+1
!                         end if
!                         !print*, "1:nbstar=", nbstar
!                     end select
!                 case(mobse_ktype_bh)
!                     select case(ktype2)
!                     case(mobse_ktype_bh)
!                         if(sh(i)%pd(idx(i)).ne.0d0.or.abs(sh(i)%ec(idx(i)))<1d0)then
!                             nbbh=nbbh+1
!                         end if
!                     case default
!                         nsbh=nsbh+2
!                     end select
!                 end select
!             end associate
!         end do
!         !print*, "nbstar=", nbstar
!         call bstar%init(nbstar)
!         bstar%ktype=mobse_ktype_ms
!         bstar%ktype2=mobse_ktype_ms
!         call bbh%init(nbbh)
!         bbh%ktype=mobse_ktype_bh
!         bbh%ktype2=mobse_ktype_bh
!         call sbh%init(nsbh)
!         sbh%ktype=mobse_ktype_bh

!         nbstar=0; nbbh=0;nsbh=0
!         do i=1, n
!             if(curtime<sh(i)%age(1).or.curtime>sh(i)%age(sh(i)%n))then
!                 !print*, "i=",i,"curtime not in rage", curtime, sh(i)%age(1), sh(i)%age(sh(i)%n)
!                 cycle
!             end if
!             associate(ktype1=>sh(i)%ktype(idx(i)), ktype2=>sh(i)%ktype2(idx(i)))
!                 select case(ktype1)
!                 case(mobse_ktype_ms, mobse_ktype_fcms)
!                     select case(ktype2)
!                     case(mobse_ktype_ms, mobse_ktype_fcms)
!                         if(sh(i)%pd(idx(i)).ne.0d0.or.abs(sh(i)%ec(idx(i)))<1d0)then
!                             nbstar=nbstar+1
!                             !print*, "2:nbstar=", nbstar
!                             bstar%mass(nbstar)=sh(i)%mass(idx(i))
!                             bstar%mass2(nbstar)=sh(i)%mass2(idx(i))
!                             bstar%age(nbstar)=sh(i)%age(idx(i))
!                             bstar%pd(nbstar)=sh(i)%pd(idx(i))
!                             bstar%ec(nbstar)=sh(i)%ec(idx(i))
!                             bstar%ktype(nbstar)=mobse_ktype_ms
!                             bstar%ktype2(nbstar)=mobse_ktype_ms
!                         end if
!                     end select
!                 case(mobse_ktype_bh)
!                     select case(ktype2)
!                     case(mobse_ktype_bh)
!                         if(sh(i)%pd(idx(i)).ne.0d0.or.abs(sh(i)%ec(idx(i)))<1d0)then
!                             nbbh=nbbh+1
!                             bbh%mass(nbbh)=sh(i)%mass(idx(i))
!                             bbh%mass2(nbbh)=sh(i)%mass2(idx(i))
!                             bbh%age(nbbh)=sh(i)%age(idx(i))
!                             bbh%pd(nbbh)=sh(i)%pd(idx(i))
!                             bbh%ec(nbbh)=sh(i)%ec(idx(i))
!                             bbh%ktype(nbbh)=mobse_ktype_bh
!                             bbh%ktype2(nbbh)=mobse_ktype_bh
!                             !bbh%kwtype(nbbh)=sh(i)%kwtype(idx(i))
!                             !call sh(i)%print()
!                             !read(*,*)
!                         end if
!                         if(sh(i)%pd(idx(i)).eq.0d0.or.abs(sh(i)%ec(idx(i))).ge.1d0)then
!                             !print*, "idx=",idx(i)
!                             !print*, sh(i)%pd(idx(i)), abs(sh(i)%ec(idx(i))), sh(i)%kwtype(idx(i))
!                             !call sh(i)%print()
!                             !read(*,*)
!                             nsbh=nsbh+1
!                             sbh%mass(nsbh)=sh(i)%mass(idx(i))
!                             sbh%age(nsbh)=sh(i)%age(idx(i))
!                             !sbh%kwtype(nsbh)=sh(i)%kwtype(idx(i))
!                             nsbh=nsbh+1
!                             sbh%mass(nsbh)=sh(i)%mass(idx(i))
!                             sbh%age(nsbh)=sh(i)%age(idx(i))
!                             !sbh%kwtype(nsbh)=sh(i)%kwtype(idx(i))
!                         end if
!                     end select
!                 end select
!             end associate
!         end do
!     end subroutine

	
! end module

! module md_mobse_stellar
! 	use md_mobse_stellar_single
! 	! use md_mobse_stellar_by
! 	implicit none
!     type stellar_hs_samples
!         !type(stellar_object_arr)::cstar, csbh, cwd, cns
!         !type(stellar_object_arr_by)::cbstar, cbbh, cbbh_u
!         !type(stellar_object_hs)::hstar, hsbh, hwd, hns
!         !type(stellar_object_hs_by)::hbstar, hbbh
!         type(stellar_object_hs), allocatable::shsg(:)
!         ! type(stellar_object_hs_by), allocatable::shby(:)
!         type(stellar_object_hs)::star0
!         ! type(stellar_object_hs_by):: bstar0
!         integer n_sg
! 		! integer n_by
!         contains
!         procedure::output_bin=>output_bin_stellar_hs_samples
!         procedure::input_bin=>input_bin_stellar_hs_samples
!     endtype

!     type stellar_samples
!         type(stellar_object_arr)::sbh, star, ns, wd, bd, srg, fgb
!         ! type(stellar_object_arr_by)::bstar,bbh, bbh_u
!         contains
!         procedure::output_bin=>output_bin_stellar_samples
!         procedure::input_bin=>input_bin_stellar_samples
        
!     end type
! 	type obj_massbin
! 		type(stellar_object_arr_by),allocatable::ba(:)
! 		type(stellar_object_arr),allocatable::a(:)
! 		real(8),allocatable:: mc(:), m1(:),m2(:)
! 		!integer,allocatable:: nstar(:), nsbh(:), nbstar(:), nbbh(:), nbbh_u(:)
! 		integer n
! 		real(8),allocatable:: nstar(:), nsbh(:), nbstar(:), nbbh(:), nbbh_u(:)
! 		real(8),allocatable:: nbd(:), nns(:), nwd(:), nrg(:), nfgb(:)
! 		real(8) nstar_tot, nbbh_tot, nsbh_tot, nbstar_tot
! 		real(8) nbd_tot, nns_tot, nwd_tot, nrg_tot, nfgb_tot
! 		type(stellar_samples_sts),allocatable::ss_sts(:)
! 		type(stellar_samples_sts)::ss_sts0
! 		real(8),allocatable::asymptot_ini(:), asymp(:,:)
! 		contains
! 			procedure:: init=>init_obj_massbin
! 			procedure::get_splited_mbins
! 			procedure::write_bin=>write_bin_obj_massbin
! 			procedure::read_bin=>read_bin_obj_massbin
! 			procedure::print=>print_obj_massbin
! 			procedure::write_asym=>write_asym_obj_massbin
! 	end type
! 	private::init_obj_massbin
! 	private::write_bin_obj_massbin
! 	private::read_bin_obj_massbin
!     private::print_obj_massbin
!     private::output_bin_stellar_hs_samples
! 	private::input_bin_stellar_hs_samples
!     private::write_asym_obj_massbin
!     private::output_bin_stellar_samples
!     private::input_bin_stellar_samples
!     real(8)::IMF_t(4),IMF_c(3), IMF_q
!     integer,parameter::nrecord_max=2000
!     real(8),parameter::alpha_IMF_triple_break(3)=(/-0.3d0,-1.3d0,-1.6d0/)
!     real(8),parameter::xb_IMF_triple_break(4)=(/0.01d0,0.08d0,0.5d0,150d0/)
! contains
! 	!=======obj_massbin===========================
! 	subroutine init_obj_massbin(this,n)
! 		implicit none
! 		class(obj_massbin)::this
! 		integer n
! 		if(allocated(this%ba))deallocate(this%ba,this%a, &
! 			this%m1, this%m2, this%mc, &
! 			this%nsbh, this%nbstar, this%nstar, this%nbbh,this%asymptot_ini, &
! 			this%asymp, this%nbbh_u, this%nbd, this%nwd,this%nns)
! 		allocate(this%ba(n),this%a(n), this%m1(n), this%m2(n), this%mc(n), &
! 			this%nsbh(n), this%nbstar(n), this%nstar(n), this%nbbh(n), &
! 			this%ss_sts(n),this%asymp(10,n),this%asymptot_ini(n),this%nbbh_u(n),&
! 			this%nbd(n),this%nns(n),this%nwd(n),this%nrg(n), this%nfgb(n))
! 		this%n=n
! 		this%nsbh=0; this%nbstar=0; this%nstar=0; this%nbbh=0; this%nbbh_u=0
! 		this%nbd=0; this%nwd=0; this%nns=0; this%nrg=0; this%nfgb=0
! 	end subroutine
! 	subroutine write_bin_obj_massbin(this, fl)
! 		implicit none
! 		class(obj_massbin)::this
! 		character*(*) fl
! 		integer::funit=19999, i

! 		open(unit=funit,file=trim(adjustl(fl))//".bin", form='unformatted', access='stream')
! 		write(funit) this%n, this%nstar_tot, this%nbstar_tot, this%nsbh_tot, this%nbbh_tot, &
! 			this%nns_tot, this%nwd_tot, this%nbd_tot
! 		write(funit) this%mc(1:this%n), this%m1(1:this%n), this%m2(1:this%n), &
! 			this%nstar(1:this%n), this%nbstar(1:this%n), this%nsbh(1:this%n), &
! 			this%nbbh(1:this%n), this%nns(1:this%n), this%nwd(1:this%n)&
!             ,this%nbd(1:this%n), this%nrg(1:this%n), this%nfgb(1:this%n)
! 		!print*, "writebinobjmassbin", this%n, this%nstar_tot, this%nbstar_tot, this%nsbh_tot, &
! 		!	 this%nbbh_tot
! 		do i=1, this%n			
! 			call this%ba(i)%write_info(funit)
! 			!print*, i, this%ba(i)%n
! 			call this%a(i)%write_info(funit)
! 			!print*, i, this%a(i)%n
! 		end do
! 		close(funit)
! 	end subroutine
! 	subroutine read_bin_obj_massbin(this, fl)
! 		implicit none
! 		class(obj_massbin)::this
! 		character*(*) fl
! 		integer::funit=19999, i

! 		open(unit=funit,file=trim(adjustl(fl))//".bin", form='unformatted', access='stream', &
! 			status='old')
! 		read(funit) this%n, this%nstar_tot, this%nbstar_tot, this%nsbh_tot, this%nbbh_tot, &
! 		this%nns_tot, this%nwd_tot, this%nbd_tot
! 		!print*, "readbinobjmassbin", this%n, this%nstar_tot, this%nbstar_tot, this%nsbh_tot, &
! 		!	 this%nbbh_tot
! 		call this%init(this%n)
! 		read(funit) this%mc(1:this%n), this%m1(1:this%n), this%m2(1:this%n), &
! 			this%nstar(1:this%n), this%nbstar(1:this%n), this%nsbh(1:this%n), &
! 			this%nbbh(1:this%n), this%nns(1:this%n), this%nwd(1:this%n),&
!             this%nbd(1:this%n), this%nrg(1:this%n), this%nfgb(1:this%n)
! 			!print*, "this%n=",this%n
! 		do i=1, this%n			
! 			call this%ba(i)%read_info(funit)
! 			!print*, i, this%ba(i)%n
! 			!print*, this%ba(i)%n
! 			call this%a(i)%read_info(funit)
! 			!print*, i, this%a(i)%n
! 		end do
! 		close(funit)
! 	end subroutine
! 	!=============================================================
! 	subroutine output_bin_stellar_hs_samples(this, fl)
!         implicit none
!         class(stellar_hs_samples)::this
!         character*(*) fl
!         integer::funit=19999, i

!         open(unit=funit,file=trim(adjustl(fl))//".bin", form='unformatted', access='stream')
!        ! print*, "1"
!         !print*, this%star0%n
!         !print*, allocated(this%star0%mass), allocated(this%star0%age), &
!         !allocated(this%star0%w)
!         write(funit) mobse_include_by, this%n_sg, this%n_by
!         call this%star0%write_info(funit)
!        ! print*, "2"
! 		if(mobse_include_by.ge.1)then
!         	call this%bstar0%write_info(funit)
! 		end if
!         !call this%hstar%write_info(funit)
!         !call this%hbstar%write_info(funit)
!         !call this%hsbh%write_info(funit)
!         !call this%hbbh%write_info(funit)
!        ! print*, "3"
!         do i=1, this%n_sg
!             call this%shsg(i)%write_info(funit)
!         end do
!         do i=1, this%n_by
!             call this%shby(i)%write_info(funit)
!             !print*, this%shby(i)%mass(:)
!         end do
!         !call this%cstar%write_info(funit)
!         !call this%cbstar%write_info(funit)
!         !call this%csbh%write_info(funit)
!         !call this%cbbh%write_info(funit)        
!         close(funit)
!     end subroutine
!     subroutine input_bin_stellar_hs_samples(this, fl)
!         implicit none
!         class(stellar_hs_samples)::this
!         character*(*) fl
!         integer::funit=19999, i

!         open(unit=funit,file=trim(adjustl(fl))//".bin", form='unformatted', access='stream', &
!             status='old')
        
!         !print*, this%z, this%by_frac, this%alpha_IMF, this%mmin, this%mmax, &
!         !        this%n_sg_simu, this%n_by_simu
!         !read(*,*)
!         read(funit) mobse_include_by, this%n_sg, this%n_by
!         allocate(this%shby(this%n_by))
!         allocate(this%shsg(this%n_sg))
!         call this%star0%read_info(funit)
!         !print*, "this%star0=",this%star0%mass(:)
! 		if(mobse_include_by.ge.1)then
!         	call this%bstar0%read_info(funit)
! 		end if
!         !print*, "this%bstar0=",this%bstar0%mass(:)
!         !read(*,*)
!         !call this%hstar%read_info(funit)
!         !call this%hbstar%read_info(funit)
!         !call this%hsbh%read_info(funit)
!         !call this%hbbh%read_info(funit)
!         do i=1, this%n_sg
!             call this%shsg(i)%read_info(funit)
!         !    print*, this%shsg(i)%mass(:)
!         end do
!         !read(*,*)
!         do i=1, this%n_by
!             call this%shby(i)%read_info(funit)
!         !    print*, this%shby(i)%mass(:)
!         end do
!         !read(*,*)
!         !call this%cstar%read_info(funit)
!         !call this%cbstar%read_info(funit)
!         !call this%csbh%read_info(funit)
!         !call this%cbbh%read_info(funit)    
!         close(funit)
!     end subroutine
! 	subroutine output_bin_stellar_samples(this, fl)
!         implicit none
!         class(stellar_samples)::this
!         character*(*) fl
!         integer::funit=19999, i

!         open(unit=funit,file=trim(adjustl(fl))//".bin", form='unformatted', access='stream')
!         call this%star%write_info(funit)
!         call this%bstar%write_info(funit)
!         call this%sbh%write_info(funit)
! 		call this%wd%write_info(funit)
! 		call this%ns%write_info(funit)
! 		call this%bd%write_info(funit)
!         call this%bbh%write_info(funit)
!         call this%bbh_u%write_info(funit)
!         call this%srg%write_info(funit)
!         call this%fgb%write_info(funit)
!         close(funit)
!     end subroutine
!     subroutine input_bin_stellar_samples(this, fl)
!         implicit none
!         class(stellar_samples)::this
!         character*(*) fl
!         integer::funit=19999, i

!         open(unit=funit,file=trim(adjustl(fl))//".bin", form='unformatted',&
!              access='stream',status='old')
!         call this%star%read_info(funit)
!         call this%bstar%read_info(funit)
!         call this%sbh%read_info(funit)
! 		call this%wd%read_info(funit)
! 		call this%ns%read_info(funit)
! 		call this%bd%read_info(funit)
!         call this%bbh%read_info(funit)
!         call this%bbh_u%read_info(funit)
!         call this%srg%read_info(funit)
!         call this%fgb%read_info(funit)
!         close(funit)
!     end subroutine
! 	subroutine write_asym_obj_massbin(om, fl)
!         implicit none
!         class(obj_massbin)::om
!         character*(*) fl
!         integer::i, funit=19989
!         open(file=trim(adjustl(fl))//"_asym.txt",unit=funit)
!         write(unit=funit, fmt="(20A15)") "m1", "mc", "m2", "asymptot_ini", "asym_star", &
!            "asym_sbh","asym_ns","asym_wd","asym_bd", "asym_rg", "asym_fgb", &
!            "asym_bstar", "asym_bbh"
!         do i=1, om%n
!             write(unit=funit, fmt="(20F15.5)") om%m1(i), om%mc(i), om%m2(i),&
!                 om%asymptot_ini(i), om%asymp(:,i)
!         end do
!     end subroutine
! 	subroutine get_splited_mbins(mb, star, bstar,sbh, ns, wd, bd,rg, fgb, bbh, bbh_u)
! 		implicit none
! 		class(obj_massbin)::mb
! 		type(stellar_object_arr)::star,sbh, ns, wd, bd, rg, fgb
! 		type(stellar_object_arr_by)::bstar, bbh, bbh_u
! 		integer i, j,  idxa,idxba
! 		!idxa=0;
! 		do i=1, mb%n
! 			mb%ba(i)%n=0
!             mb%a(i)%n=0
! 			do j=1, star%n
! 				if(star%mass(j).ge.mb%m1(i).and.star%mass(j).le.mb%m2(i))then
! 					mb%a(i)%n=mb%a(i)%n+1			
!                     mb%nstar(i)=mb%nstar(i)+star%w(j)		
! 				end if
! 			end do
! 			do j=1, bstar%n
! 				if(bstar%mass_tot(j).ge.mb%m1(i).and.bstar%mass_tot(j).le.mb%m2(i))then
! 					!mb%ba(i)%n=mb%ba(i)%n+1	
!                     mb%nbstar(i)=mb%nbstar(i)+bstar%w(j)				
! 				end if
! 			end do
! 			do j=1, sbh%n
! 				if(sbh%mass(j).ge.mb%m1(i).and.sbh%mass(j).le.mb%m2(i))then
! 					mb%a(i)%n=mb%a(i)%n+1					
!                     mb%nsbh(i)=mb%nsbh(i)+sbh%w(j)
! 				end if
! 			end do
! 			do j=1, ns%n
! 				if(ns%mass(j).ge.mb%m1(i).and.ns%mass(j).le.mb%m2(i))then
!                     mb%a(i)%n=mb%a(i)%n+1	
!                     mb%nns(i)=mb%nns(i)+ns%w(j)
! 				end if
! 			end do
! 			do j=1, wd%n
! 				if(wd%mass(j).ge.mb%m1(i).and.wd%mass(j).le.mb%m2(i))then
!                     mb%a(i)%n=mb%a(i)%n+1	
!                     mb%nwd(i)=mb%nwd(i)+wd%w(j)
! 				end if
! 			end do
! 			do j=1, bd%n
! 				!print*, "bd%mass(j)=", bd%mass(j)
! 				if(bd%mass(j).ge.mb%m1(i).and.bd%mass(j).le.mb%m2(i))then
                    
!                     mb%nbd(i)=mb%nbd(i)+bd%w(j)
! 				end if
! 			end do
!             do j=1, rg%n
! 				!print*, "bd%mass(j)=", bd%mass(j)
! 				if(rg%mass(j).ge.mb%m1(i).and.rg%mass(j).le.mb%m2(i))then
!                     mb%a(i)%n=mb%a(i)%n+1	
!                     mb%nrg(i)=mb%nrg(i)+rg%w(j)
! 				end if
! 			end do
!             do j=1, fgb%n
! 				!print*, "bd%mass(j)=", bd%mass(j)
! 				if(fgb%mass(j).ge.mb%m1(i).and.fgb%mass(j).le.mb%m2(i))then
!                     mb%a(i)%n=mb%a(i)%n+1	
!                     mb%nfgb(i)=mb%nfgb(i)+fgb%w(j)
! 				end if
! 			end do

! 			do j=1, bbh%n
! 				if(bbh%mass_tot(j).ge.mb%m1(i).and.bbh%mass_tot(j).le.mb%m2(i))then
! 					!mb%ba(i)%n=mb%ba(i)%n+1		
!                     mb%nbbh(i)=mb%nbbh(i)+bbh%w(j)		
! 				end if
! 			end do
!             !print*, "bbh_u%n=", bbh_u%n
!             do j=1, bbh_u%n
!                ! print*, "bbh_u%mass_tot(j)=",bbh_u%mass_tot(j), mb%m1(i), mb%m2(i)
! 				if(bbh_u%mass_tot(j).ge.mb%m1(i).and.bbh_u%mass_tot(j).le.mb%m2(i))then
! 					mb%ba(i)%n=mb%ba(i)%n+1		
!                     mb%nbbh_u(i)=mb%nbbh_u(i)+bbh_u%w(j)		
!                     !print*, "+"
! 				end if
! 			end do
!             !print*, "mb%nbbh_u(i)=",mb%nbbh_u(i)
! 			call mb%ba(i)%init(mb%ba(i)%n)
!             call mb%a(i)%init(mb%a(i)%n)
! 		end do
!         mb%nstar_tot=sum(mb%nstar)
!         mb%nsbh_tot=sum(mb%nsbh)
!         mb%nbstar_tot=sum(mb%nbstar)
!         mb%nbbh_tot=sum(mb%nbbh)
! 		mb%nns_tot=sum(mb%nns)
! 		mb%nwd_tot=sum(mb%nwd)
! 		mb%nbd_tot=sum(mb%nbd)
!         mb%nrg_tot=sum(mb%nrg)
!         mb%nfgb_tot=sum(mb%nfgb)
! 		do i=1, mb%n
! 			idxa=0;idxba=0
! 			do j=1, star%n
! 				if(star%mass(j).ge.mb%m1(i).and.star%mass(j).le.mb%m2(i))then
! 					idxa=idxa+1					
! 					mb%a(i)%ktype(idxa)=star%ktype(j)
! 					mb%a(i)%mass(idxa)=star%mass(j)
! 					!mb%ba(i)%ktype2(idx)=-1		
!                     mb%a(i)%w(idxa)=star%w(j)	
!                     mb%a(i)%radius(idxa)=star%radius(j)
! 				end if
! 			end do
!             do j=1, rg%n
! 				if(rg%mass(j).ge.mb%m1(i).and.rg%mass(j).le.mb%m2(i))then
! 					idxa=idxa+1					
! 					mb%a(i)%ktype(idxa)=rg%ktype(j)
! 					mb%a(i)%mass(idxa)=rg%mass(j)
! 					!mb%ba(i)%ktype2(idx)=-1		
!                     mb%a(i)%w(idxa)=rg%w(j)	
!                     mb%a(i)%radius(idxa)=rg%radius(j)
! 				end if
! 			end do
!             do j=1, fgb%n
! 				if(fgb%mass(j).ge.mb%m1(i).and.fgb%mass(j).le.mb%m2(i))then
! 					idxa=idxa+1					
! 					mb%a(i)%ktype(idxa)=fgb%ktype(j)
! 					mb%a(i)%mass(idxa)=fgb%mass(j)
! 					!mb%ba(i)%ktype2(idx)=-1		
!                     mb%a(i)%w(idxa)=fgb%w(j)	
!                     mb%a(i)%radius(idxa)=fgb%radius(j)
! 				end if
! 			end do
!             do j=1, ns%n
! 				if(ns%mass(j).ge.mb%m1(i).and.ns%mass(j).le.mb%m2(i))then
! 					idxa=idxa+1					
! 					mb%a(i)%ktype(idxa)=ns%ktype(j)
! 					mb%a(i)%mass(idxa)=ns%mass(j)
! 					!mb%ba(i)%ktype2(idx)=-1		
!                     mb%a(i)%w(idxa)=ns%w(j)	
!                     mb%a(i)%radius(idxa)=ns%radius(j)
! 				end if
! 			end do
!             do j=1, wd%n
! 				if(wd%mass(j).ge.mb%m1(i).and.wd%mass(j).le.mb%m2(i))then
! 					idxa=idxa+1					
! 					mb%a(i)%ktype(idxa)=wd%ktype(j)
! 					mb%a(i)%mass(idxa)=wd%mass(j)
! 					!mb%ba(i)%ktype2(idx)=-1		
!                     mb%a(i)%w(idxa)=wd%w(j)	
!                     mb%a(i)%radius(idxa)=wd%radius(j)
! 				end if
! 			end do
! 			do j=1, bstar%n
! 				!if(bstar%mass_tot(j).ge.mb%m1(i).and.bstar%mass_tot(j).le.mb%m2(i))then
! 				!	idx=idx+1					
! 				!	mb%ba(i)%ktype(idx)=bstar%ktype(j)
! 				!	mb%ba(i)%mass(idx)=bstar%mass(j)
! 				!	mb%ba(i)%ktype2(idx)=bstar%ktype2(j)
! 				!	mb%ba(i)%mass2(idx)=bstar%mass2(j)
! 				!	mb%ba(i)%ac(idx)=bstar%ac(j)
! 				!	mb%ba(i)%ec(idx)=bstar%ec(j)
!                 !    mb%ba(i)%w(idx)=bstar%w(j)
! 				!end if
! 			end do
! 			do j=1, sbh%n
! 				if(sbh%mass(j).ge.mb%m1(i).and.sbh%mass(j).le.mb%m2(i))then
! 					idxa=idxa+1					
! 					mb%a(i)%ktype(idxa)=sbh%ktype(j)
! 					mb%a(i)%mass(idxa)=sbh%mass(j)
! 					!mb%ba(i)%ktype2(idx)=-1		
!                     mb%a(i)%w(idxa)=sbh%w(j)	
!                     mb%a(i)%radius(idxa)=sbh%radius(j)
! 				end if
! 			end do
! 			do j=1, bbh_u%n
! 				if(bbh_u%mass_tot(j).ge.mb%m1(i).and.bbh_u%mass_tot(j).le.mb%m2(i))then
! 					idxba=idxba+1					
! 					mb%ba(i)%ktype(idxba)=bbh_u%ktype(j)
! 					mb%ba(i)%mass(idxba)=bbh_u%mass(j)
! 					mb%ba(i)%ktype2(idxba)=bbh_u%ktype2(j)
! 					mb%ba(i)%mass2(idxba)=bbh_u%mass2(j)
!                     mb%ba(i)%mass_tot(idxba)=bbh_u%mass_tot(j)
!                    ! print*, "mb%mass_tot=", bbh_u%mass_tot(j), mb%ba(i)%mass_tot(idxba)
! 					mb%ba(i)%ac(idxba)=bbh_u%ac(j)
! 					mb%ba(i)%ec(idxba)=bbh_u%ec(j)		
!                     mb%ba(i)%w(idxba)=bbh_u%w(j)	
! 				end if
! 			end do
! 		end do
! 	end subroutine
!     subroutine print_obj_massbin(mb)
!         implicit none
!         class(obj_massbin)::mb
!         integer i
!         print*, "weighted"
!         write(*, fmt="(13A15)") "m1","mc","m2", "nstar","nsbh","nns","nwd","nbd", &
!             "nrg", "nfgb", "nbstar", "nbbh", "nbbh_u"
!         do i=1, mb%n
!             write(*, fmt="(3F15.3, 10F15.3)") mb%m1(i),mb%mc(i), mb%m2(i), &
!                 mb%nstar(i), mb%nsbh(i), mb%nns(i), mb%nwd(i), mb%nbd(i), &
!                 mb%nrg(i), mb%nfgb(i), mb%nbstar(i), mb%nbbh(i),mb%nbbh_u(i)
!         end do
!         write(*, fmt="(10A15)")  "tot","nstar", "nsbh", "nns","nwd","nbd",&
!         "nrg", "nfgb", "nbstar", "nbbh"
!         write(*, fmt="(A15,10F15.3)")  "", mb%nstar_tot, mb%nsbh_tot, mb%nns_tot,&
!          mb%nwd_tot,mb%nbd_tot, mb%nrg_tot, mb%nfgb_tot, &
! 					 mb%nbstar_tot, mb%nbbh_tot
!         print*, "unweighted"
!         write(*, fmt="(3A15,5A9)") "m1","mc","m2","nsgtot", "nbbh"
!         do i=1, mb%n
!             write(*, fmt="(3F15.4, 5I9)") mb%m1(i),mb%mc(i), mb%m2(i), mb%a(i)%n,&
!              mb%ba(i)%n
!         end do
!     end subroutine
! end module