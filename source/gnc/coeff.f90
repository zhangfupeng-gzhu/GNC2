module md_coeff
	!use my_intgl
	!use model_basic
	use com_sts_type
	implicit none
	integer,parameter::Invns=101
	real(8),parameter::inf=-50.1d0, tiny_=1d-9
    real(8):: emin_factor, emax_factor
	real(8)::emin_dstr_factor,emax_dstr_factor
	real(8):: true_log10emax_factor
    real(8):: log10emin_factor, log10emax_factor
	type coeff_type
		real(8) e,j,ee,jj,ej, e_110,e_0,j_111,j_rest, m_avg
	end type 

    integer,parameter::Jbin_type_log=2 
	integer,parameter::ebin_type_log=2
	type diffuse_coeffient_type
		type(s2d_type)::s2_de, s2_dee,s2_dj, s2_djj,s2_dej
 
		type(s2d_type)::s2_de_110, s2_de_0, s2_dj_111, s2_dj_rest
        type(s1d_type)::s1_de, s1_dee
        real(8) emin, emax, jmin, jmax
		integer nbin, jbin_type		
		contains
		procedure::init=>init_diffuse_coeffient_grid
		procedure::write_grid=>write_diffuse_coeffient_grid
		procedure::read_grid=>read_diffuse_coeffient_grid
		!procedure::get_diffuse_coeff=>get_ba16_diffuse_coeff
	end type
	type(diffuse_coeffient_type),allocatable::df(:)
	type(diffuse_coeffient_type)::df_tot
	 
	integer coeff_chattery
	private::init_diffuse_coeffient_grid
	integer,parameter::coeff_sts_type_dc=sts_type_dstr
contains
	 
	subroutine init_diffuse_coeffient_grid(this, nbin, emin, emax, jmin,jmax, ebin_type,jbin_type)
		implicit none
		class(diffuse_coeffient_type)::this
		!integer,parameter::sts_type_dc=sts_type_grid
		integer ::nbin
		real(8) emin, emax
		real(8) jmin, jmax
        real(8) tmin, tmax, smin, smax
        integer jbin_type,ebin_type

		this%nbin=nbin
		this%emin=emin; this%emax=emax; this%jmin=jmin; this%jmax=jmax
        this%jbin_type=jbin_type
        select case(jbin_type) 
        case(jbin_type_log)
            tmin=log10(jmin)
            tmax=log10(jmax) 
        end select
		select case(ebin_type)
		case(ebin_type_log)
			smin=log10(emin)
			smax=log10(emax) 
		case default
			print*, "error! ini_dffuse_coeffient"
			stop
		end select

		call this%s2_de_110%init(nbin, nbin, smin, smax, tmin,tmax, coeff_sts_type_dc)
		call this%s2_de_110%set_range()

		call this%s2_de_0%init( nbin, nbin, smin, smax, tmin,tmax, coeff_sts_type_dc)
		call this%s2_de_0%set_range()
        !print*, tmin,tmax
        !print*, this%s2_de_0%ycenter
        !stop
		call this%s2_dee%init(nbin, nbin, smin, smax, tmin,tmax, coeff_sts_type_dc)
		call this%s2_dee%set_range()

		call this%s2_dj_111%init(nbin, nbin, smin, smax, tmin,tmax, coeff_sts_type_dc)
		call this%s2_dj_111%set_range()

		call this%s2_dj_rest%init(nbin, nbin, smin, smax, tmin,tmax, coeff_sts_type_dc)
		call this%s2_dj_rest%set_range()

		call this%s2_djj%init(nbin, nbin, smin, smax, tmin,tmax, coeff_sts_type_dc)
		call this%s2_djj%set_range()

		call this%s2_dej%init(nbin, nbin, smin, smax, tmin,tmax, coeff_sts_type_dc)
		call this%s2_dej%set_range()
  
		
        call  this%s1_de%init(smin,smax, nbin, coeff_sts_type_dc)
        this%s1_de%xb=this%s2_de_0%xcenter
!
        call this%s1_dee%init(smin,smax, nbin, coeff_sts_type_dc)
        this%s1_dee%xb=this%s2_dee%xcenter
	end subroutine
	subroutine write_diffuse_coeffient_grid(this, file_unit)
		implicit none
		class(diffuse_coeffient_type)::this
		integer file_unit
		write(unit=file_unit) this%nbin
		write(unit=file_unit) this%s2_de_110,this%s2_de_0,this%s2_dee
		write(unit=file_unit) this%s2_dj_111, this%s2_dj_rest, this%s2_djj
		write(unit=file_unit) this%s2_dej 
	end subroutine
	subroutine read_diffuse_coeffient_grid(this, file_unit)
		implicit none
		class(diffuse_coeffient_type)::this
		integer file_unit
		read(unit=file_unit) this%nbin
		read(unit=file_unit) this%s2_de_110,this%s2_de_0,this%s2_dee
		read(unit=file_unit) this%s2_dj_111, this%s2_dj_rest, this%s2_djj
		read(unit=file_unit) this%s2_dej 
	end subroutine
	      
end module
           
