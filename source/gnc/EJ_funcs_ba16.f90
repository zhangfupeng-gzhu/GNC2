
subroutine fgx_mb_star(x, fx)
	use model_basic
	implicit none
	real(8) x, fx
	!real(8) yout(6)
	!print*, fc_share%xb(fc_share%nbin)
	!read(*,*)
	if(x>=10**fc_share%xb(1).and.x<=10**fc_share%xb(fc_share%nbin))then
		!call get_value_at_x_fc(fc_share, log10(x), yout, 1)
		call fc_share%get_value_l(log10(x),fx)
		!fx=yout(1)
	else if(x<=0) then
		fx=exp(x)*fgx_g0
	else if(x>10**fc_share%xb(fc_share%nbin).or.(x.ge.0.and.x<10**fc_share%xb(1)))then
		fx=0d0		
	end if
	!call fc_share%print()
	!print*, "x, fx=", x, fx
	!read(*,*)
end subroutine

subroutine fgx_mb_ir(x, fx,xmin,xmax)
	use model_basic
	implicit none
	real(8) x, fx
	!real(8) yout(6)
    real(8) yout, xmin,xmax


	if(x>=xmin.and.x<=xmax)then 
		select case(ctl%ebin_type)
		case(ebin_type_log)
        	call linear_int_fast(fc_ir_share%xb(1:fc_ir_share%nbin),fc_ir_share%fx(1:fc_ir_share%nbin),&
			fc_ir_share%nbin, log10(x),yout) 
		end select
 
		fx=yout
		if(yout<0)fx=0 
	else if(x<xmin) then 
		fx=0d0
	else if(x>xmax)then
		call linear_int_fast(fc_ir_share%xb(1:fc_ir_share%nbin),fc_ir_share%fx(1:fc_ir_share%nbin),&
			fc_ir_share%nbin, fc_ir_share%xmax,yout) 
		fx=max(abs(yout)*(1-(x-xmax)*1d7),0d0) 
	end if 
end subroutine
subroutine fgx_mb(x, fx)
	use model_basic
	implicit none
	real(8) x, fx 
    real(8) yout 
	if(x>=sample_emin.and.x<=sample_emax)then 
		select case(ctl%ebin_type)
		case(ebin_type_log)
        	call linear_int_fast(fc_share%xb(1:fc_share%nbin),fc_share%fx(1:fc_share%nbin),&
	            fc_share%nbin, log10(x),yout) 
		end select
		 
		fx=yout
		if(yout<0)fx=0 
	else if(x<sample_emin) then
		fx=0d0!exp(x)*fgx_g0
 
	else if(x>sample_emax)then
		fx=0d0
	end if
 
end subroutine
  