 
     
subroutine get_sts_single(sps, sts)
    use md_sts
	use com_main_gw
	implicit none
	type(particle_samples_arr_type)::sps
	type(sts_single_type)::sts
    integer,parameter::rxn=30
	integer i,n
	real(8) xmin,xmax,rsw,mmin, mmax
    real(8),allocatable::w(:),x(:), gw_distance(:)

   	if(sps%n.eq.0)then
		print*, "no samples for sts, skiped"
		return
	end if
	allocate(w(sps%n),x(sps%n), gw_distance(sps%n))
	w=sps%sp(1:sps%n)%weight_real
    
	do i=1, sps%n
		if(ieee_is_nan(w(i))) then
			print*, "ieee_is_nan:i=",i
			stop
		end if
	end do
	n=sps%n

	x(1:n)=sps%sp(1:n)%byot%ms%m
    mmin=minval(x(1:n))
    mmax=maxval(x(1:n))
	call sts%mtot%init(mmin,mmax, rxn, fc_spacing_log,use_weight=.true.)
    call get_fc_weight(x(1:n), w(1:n) , n, sts%mtot)

end subroutine
