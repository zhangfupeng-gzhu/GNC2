module md_particle
!	use init_constant
!	use commonfunc
	use ieee_arithmetic
	use constant
	type::particle
        real(8) :: x(3)
        real(8) :: vx(3)
        real(8) spin(3) ! spin vector in x, y, z, 
                        ! for stars, WD and NS, the amplitude is the rotation angular velocity
                        ! for black holes, the amplitude is the spin in unit of M^2 (maximum=1)
        real(8) :: M
        real(8) radius
 		integer id, obtype, obidx !, source  ! where it comes from
		integer N_gene!, I_FOR_ALLIGMENT
	end type
	integer::PARTICLE_MPI_TYPE_
contains
	real(8) function get_distance(P1,P2)
		type(particle)::P1,P2
		get_distance=sqrt((P1%x(1)-P2%x(1))**2+(P1%x(2)-P2%x(2))**2+(P1%x(3)-P2%x(3))**2)
	end function
	real(8) function get_VelMag(P)
		type(particle)::P
		get_VelMag=sqrt(P%vx(1)**2+P%vx(2)**2+P%vx(3)**2)
	end function	
	real(8) function get_RMag(P)	
		type(particle)::P
		get_RMag=sqrt(P%x(1)**2+P%x(2)**2+P%x(3)**2)
	end function
end module

module md_pts
	use md_particle
	type psys_type
		type(particle),allocatable:: pts(:)
		integer N
	end type
contains
	subroutine init_psys(psys,N)
		implicit none
		type(psys_type)::psys
		integer n
		if(allocated(psys%pts))deallocate(psys%pts)
		psys%N=N
		allocate(psys%pts(N))
	end subroutine
end module

module md_binary
	use md_particle
	type::binary
		type(particle)::Ms
		type(particle)::Mm
		type(particle)::rd           !the relative position and velocity, from massive to light one
		real(8) E,l,k,miu,Mtot, Jc!, tl  !tl is the true longitude (tl=pe+f0)
		real(8) a_bin,e_bin, lum(3), f0
		real(8) Inc, Om, pe, t0, me ! Inc, Om, pe in unit of degree
		integer an_in_mode
		character*(100) bname
	contains
	!	procedure::print=>print_binary 
	!	procedure::get_period=>by_get_period
	!---
	!---be careful when e=0
	!---
	end type
	integer,parameter::an_in_mode_f0=1, an_in_mode_t0=2, an_in_mode_mean=3 
	integer::BINARY_MPI_TYPE_
contains  
 

	subroutine by_em2st(by)
		use constant
		implicit NONE
		type(binary)::by
		integer way_an_in
		real(8) gm, q,e,inc,p,n,l,x,y,z,u,v,w
		real(8) ecc_ano, r,theta,phi
		real(8) pos(3),vel(3), outvel(3)
		!if by%e_bin=0, we should avoid an_in_mode to be an_in_mode_f0

!  gm = grav const * (central + secondary mass)
!  q = perihelion distance
!  e = eccentricity
!  i = inclination                 )
!  p = longitude of perihelion !!! )   in
!  n = longitude of ascending node ) radians
!  l = mean anomaly       
!  x axis points to n=0
!  velocity in units of sqrt(gm/a) 

		gm=by%ms%m+by%mm%m;
		e=by%e_bin
		q=by%a_bin*(1-e)
		inc=by%inc
		n=by%om
		p=by%pe+n
		
		select case(by%an_in_mode)
			case (an_in_mode_t0)
				l=-sqrt(by%Mtot/by%a_bin**3)*by%t0
			case (an_in_mode_mean)
				l=by%me
				by%t0=-l/sqrt(by%Mtot/(by%a_bin**3))
			case (an_in_mode_f0)
				if(by%e_bin<1d0)then
					ecc_ano=atan(tan(by%f0/2d0)*(1-by%e_bin)**0.5/(1+by%e_bin)**0.5)*2
					if (by%f0.gt.pi-30d0/180d0*pi.and.by%f0.le.pi)then
						ecc_ano=acos((cos(by%f0)+by%e_bin)/(1+by%e_bin*cos(by%f0)))
					else if (by%f0.lt.pi+30d0/180d0*pi.and.by%f0.gt.pi)then
						ecc_ano=-acos((cos(by%f0)+by%e_bin)/(1+by%e_bin*cos(by%f0)))
					end if
					if(by%f0.eq.pi) ecc_ano= pi
					l=ecc_ano-by%e_bin*sin(ecc_ano)	
				else if(by%e_bin>1d0)then

					ecc_ano=atanh(tan(by%f0/2d0)*(by%e_bin-1)**0.5/(1+by%e_bin)**0.5)*2
					if (by%f0.gt.pi-30d0/180d0*pi.and.by%f0.lt.pi)then
						ecc_ano=acosh((cos(by%f0)+by%e_bin)/(1+by%e_bin*cos(by%f0)))
					else if (by%f0.lt.pi+30d0/180d0*pi.and.by%f0.gt.pi)then
						ecc_ano=-acosh((cos(by%f0)+by%e_bin)/(1+by%e_bin*cos(by%f0)))
					end if
					if(by%f0.eq.pi) ecc_ano=pi
					l=by%e_bin*sinh(ecc_ano)-ecc_ano	
				else
					ecc_ano=tan(by%f0/2d0)
					l=ecc_ano+ecc_ano**3/3d0
				end if
				by%t0=-l/sqrt(by%mtot/by%a_bin**3)
			case default
				print*, "ERROR: an_in_mode not defined!"
				stop
		end select
		by%me=l
		!print*, "l=",l
		call mco_el2x(gm,q,e,inc,p,n,l,x,y,z,u,v,w)

		if(by%an_in_mode.eq.an_in_mode_t0) then
			if(by%e_bin.gt.1d-6)then
				call get_true_anomaly(gm, x,y,z,u,v,w,by%f0)
			else
				by%f0=l
			end if
		end if
		!print*, x,y,z,u,v,w,by%f0,by%an_in_mode.eq.an_in_mode_t0
!		stop
		by%rd%x=(/x,y,z/);
		by%rd%vx=(/u,v,w/);
		!print*,"I, P, N, L=", inc, p ,n ,l 
		!print*, "x=",by%rd%x
		!stop
		by%rd%m=by%ms%m*by%mm%m/gm
!		by%mm%x=0;	by%mm%vx=0
		!by%rd%x=by%ms%x; by%rd%vx=by%ms%vx; 
		!print*, by%ms%x,by%ms%vx, by%mm%x,by%mm%vx
!		call by_move_to_mass_center(by)
		!print*, by%ms%x,by%ms%vx, by%mm%x,by%mm%vx

		if(ieee_is_nan(by%f0).and.by%an_in_mode.eq.an_in_mode_f0)then
			print*, "em2st:nan!"
			call print_binary(by)
!			stop
		end if

	end subroutine

	subroutine by_split_from_rd(by)
		implicit none
		type(binary)::by
		
 		by%ms%x=by%rd%x*by%mm%m/(by%ms%m+by%mm%m)
		by%mm%x=-by%rd%x*by%ms%m/(by%ms%m+by%mm%m)
		by%ms%vx=by%rd%vx*by%mm%m/(by%ms%m+by%mm%m)
		by%mm%vx=-by%rd%vx*by%ms%m/(by%ms%m+by%mm%m)
		
	end subroutine
	  
	subroutine print_binary(by,funit)
		implicit none
		class(binary)::by
        integer,optional::funit
        character*(10) str_ms_obtype, str_mm_obTYPe

        if(present(funit))then
            write(unit=funit,fmt="(2A20)") "NAME=", trim(adjustl(by%bname))
            write(unit=funit,fmt="(A20,I4)") "ANMODE=",by%an_in_mode
            write(unit=funit,fmt="(A20,1P20E20.8)") "MS,MM=",by%ms%m,by%mm%m
            call get_star_type(by%ms%obtype, str_ms_obtype)
            call get_star_type(by%mm%obtype, str_mm_obtype)
            write(unit=funit,fmt="(3A20)") "TYPE_MS,TYPE_MM=",str_ms_obtype,str_mm_obtype
            
            write(unit=funit,fmt="(A20,1P20E20.8)") "Sma,Ecc, Inc=",by%a_bin,by%e_bin,by%Inc
            write(unit=funit,fmt="(A20,1P20E20.8)") "E,L=",by%e,by%l
            if(by%a_bin<0d0)then
                write(*,fmt="(A20,1P20E20.8)") "Vinf(kms)=",sqrt(-(by%ms%m+by%mm%m)/by%a_bin)*29.79
            end if
            write(unit=funit,fmt="(A20,20F20.10)") "Ome,Pe,Me=",by%Om,by%pe,by%me
            write(unit=funit,fmt="(A20,20F20.10)") "f=",by%f0

            write(unit=funit,fmt="(A20,20F20.10)") "RD X=",by%rd%x
            write(unit=funit,fmt="(A20,20F20.10)") "RD MAG(X)=",sqrt(by%rd%x(1)**2+by%rd%x(2)**2+by%rd%x(3)**2)
            write(unit=funit,fmt="(A20,20F20.10)") "RD VX=",by%rd%vx

            write(unit=funit,fmt="(A20,20F20.10)") "MS X=",by%ms%x
            write(unit=funit,fmt="(A20,20F20.10)") "MS VX=",by%ms%vx

            write(unit=funit,fmt="(A20,20F20.10)") "MM X=",by%mm%x
            write(unit=funit,fmt="(A20,20F20.10)") "MM VX=",by%mm%vx
        else
            write(*,fmt="(2A20)") "NAME=", trim(adjustl(by%bname))
            write(*,fmt="(A20,I4)") "ANMODE=",by%an_in_mode
            write(*,fmt="(A20,1P20E20.8)") "MS,MM=",by%ms%m,by%mm%m
            call get_star_type(by%ms%obtype, str_ms_obtype)
            call get_star_type(by%mm%obtype, str_mm_obtype)
            write(*,fmt="(3A20)") "TYPE_MS,TYPE_MM=",str_ms_obtype,str_mm_obtype
            
            write(*,fmt="(A20,1P20E20.8)") "Sma,Ecc, Inc=",by%a_bin,by%e_bin,by%Inc
            write(*,fmt="(A20,1P20E20.8)") "E,L=",by%e,by%l
            if(by%a_bin<0d0)then
                write(*,fmt="(A20,1P20E20.8)") "Vinf(kms)=",sqrt(-(by%ms%m+by%mm%m)/by%a_bin)*29.79
            end if
            write(*,fmt="(A20,20F20.10)") "Ome,Pe,Me=",by%Om,by%pe,by%me
            write(*,fmt="(A20,20F20.10)") "f=",by%f0

            write(*,fmt="(A20,20F20.10)") "RD X=",by%rd%x
            write(*,fmt="(A20,20F20.10)") "RD MAG(X)=",sqrt(by%rd%x(1)**2+by%rd%x(2)**2+by%rd%x(3)**2)
            write(*,fmt="(A20,20F20.10)") "RD VX=",by%rd%vx

            write(*,fmt="(A20,20F20.10)") "MS X=",by%ms%x
            write(*,fmt="(A20,20F20.10)") "MS VX=",by%ms%vx

            write(*,fmt="(A20,20F20.10)") "MM X=",by%mm%x
            write(*,fmt="(A20,20F20.10)") "MM VX=",by%mm%vx
        end if
	end subroutine
end module
 
 