module Astron_constant
	real(8),parameter::pc=206264.98d0      !AU
	real(8),parameter::rd_sun=4.65d-3      !AU
	real(8),parameter::AU_GS=1.4959787d13  !cm
	real(8),parameter::AU_SI=1.4959787d11  !m  
	real(8),parameter::m_sun_GS=1.98855d33    !g
	real(8),parameter::m_sun_SI=1.98855d30    !kg
	real(8),parameter::one_year=3.1556d7    !s
	real(8),parameter::pi4degree2=41252.96125 ! deg**2
end module
module unitless_value
	real(8),parameter::PI=3.141592653589793d0   
	real(8),parameter::TWO_PI=3.141592653589793d0*2d0   
	real(8),parameter::e_nature=2.718281828459045d0
	real(8),parameter::moer_const=6.02214129e23 !avogadro constant
!	real(8),parameter::PI=acos(-1d0)
	real(8),parameter::alpha_FC=7.297352570d-3     !Fine-structure constant
end module
module cosmology_constant
	use unitless_value
	use Astron_constant
	!planck + WP+highL+BAO
	real(8),parameter::cos_Hub_0=67.77		!kms-1/Mpc
	real(8),parameter::cos_h_0=cos_Hub_0/100d0
	real(8),parameter::Omega_b=0.022161/cos_h_0**2  !baryonic matter
	real(8),parameter::Omega_c=0.11889/cos_h_0**2	!cold dark matter
	real(8),parameter::Omega_m=Omega_b+Omega_c
	real(8),parameter::Omega_Lambda=0.6914			!dark energy
	real(8),parameter::Omega_0=Omega_b+Omega_c+Omega_Lambda
	real(8),parameter::Sigma_8=0.8288		!
	real(8),parameter::cos_n_s=0.9611		!initial power index
	real(8),parameter::rho_crit0=3*(100/29.79)**2/8/ &
					pi*(pc*1d6)*cos_h_0**2 ! Msun/Mpc^3
end module
module SI_unit_value
	real(8),parameter::G_SI=6.67384d-11 !m^3 kg^-1 s^-2
	real(8),parameter::lumi_sun_SI=3.83d26 !W
	real(8),parameter::vel_c_SI=2.99792458d8     !m S-1
	real(8),parameter::plank_SI=6.6260696d-34  !J S
	real(8),parameter::Boltzmann_SI=1.38d-23      !J K-1       
	real(8),parameter::c_electron_SI=1.6021766d-19 !C
	real(8),parameter::m_electron_SI=9.109383d-31    !Kg
	real(8),parameter::m_proton_SI=1.672621898d-27   !Kg
	real(8),parameter::Coulomb_constant_SI=8.987552d9   !N m^2 c^-2
	real(8),parameter::electric_constant=8.854188d-12  !F m^-1
	real(8),parameter::sigma_SB_SI=5.67d-8 !W m-2 k-4
	real(8),parameter::Pulsar_time_T_odot=4.925490947e-6 ! in unit of s 
	!! T_odot GM_dot/c^3 = 4.925490947 us see 2002ASPC..278..251S
end module
module GS_unit_value
	real(8),parameter::plank_GS=6.62606885d-27     !erg*s
	real(8),parameter::plank_bar_GS=1.0545716d-27  !erg*s
	real(8),parameter::vel_c_GS=2.99792458d10      !cm/s
	real(8),parameter::Boltzmann_GS=1.3806504d-16  !erg/K
	real(8),parameter::G_GS=6.67428d-8             !cm^3/(g*s^2)
	real(8),parameter::Coulomb_constant_GS=1       !for ESU system
	real(8),parameter::c_electron_GS=4.80320427d-10 !for ESU system
	real(8),parameter::m_electron_GS=9.109383d-28   !g
	real(8),parameter::m_proton_GS=1.672621898d-24   !g
	real(8),parameter::Thomson_elec_scatter=6.6524587158d-25 ! cm^2 Thomson cross-section for election
	real(8),parameter::sigma_SB_GS=5.670367d-5   !erg cm-2 k-4
	real(8),parameter::lumi_sun_GS=3.83d33   !erg s^-1
	real(8),parameter::pc_GS=3.08568016635926d18 !cm
end module
module my_unit
	use unitless_value
	!G=1,M=1msun, Length=1AU
	!then the velocity_unit=29.79kms-1
	real(8),parameter::my_unit_vel_c=3d5/29.784d0
    real(8),parameter::my_unit_of_time=58.123 !days
	real(8),parameter::my_unit_of_energy=6.138d5 !solar_lumi also =3.53d32 W
contains
	real(8) function myu_conv_t2second(t)
		implicit none
		real(8) t
		myu_conv_t2second=t*365.2425d0/2/pi*86400      !(second)
	end function 
	 
end module
module constant 
use unitless_value     !Unitless values
use SI_unit_value      !International System of Units
use GS_unit_value      !gaussian unit constant
use Astron_constant        
use my_unit
use cosmology_constant
end module

