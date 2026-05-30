module MPI_comu
	use model_basic
	use md_sts
	use mpi
	implicit none
	!include 'mpif.h'
	integer rid, proc_id
	integer:: mpi_master_id=0

!$omp threadprivate(rid)
contains 
	subroutine init_mpi()
		implicit none
		integer ierr, getpid
		!print*, "start"
		call mpi_init(ierr)
		!print*, "init"
		call mpi_comm_rank(mpi_comm_world, rid,ierr)
		!print*, "rid=", rid
		!print*, "size",ctl%ntasks
		call mpi_comm_size(mpi_comm_world, ctl%ntasks, ierr)
		sample_threads_number=ctl%ntasks
		proc_id=getpid()
	! print*, "size", rid
		if(ctl%chattery>4)then 
			print*, "mpi initializtion finished"
		end if
        !print*, "declare particle"
		call declare_particle_mpi_type()
        !print*, "declare binary"
		call declare_binary_mpi_type()
		
	end subroutine
	subroutine stop_mpi()
		implicit none
		integer ierr
		call MPI_Type_free(PARTICLE_MPI_TYPE_, ierr)
		call MPI_Type_free(BINARY_MPI_TYPE_, ierr)
		call mpi_finalize(ierr)
	end subroutine
	subroutine collect_data_mpi(fx,nbin,ibg,ied,ntasks)
		implicit none
		integer nbin
		real(8) fx(nbin)
		real(8),allocatable::reivbuffer(:)
		real(8),allocatable::sentbuffer(:)
		integer i, j, k, ibg,ied, sent_count, ierr
		integer idbg, ided,nblocks, ntasks

		allocate(reivbuffer(nbin))
		allocate(sentbuffer(nbin))
		sentbuffer=0
		reivbuffer=0
		if(mod(nbin, ntasks).ne.0)then
			print*, "error in collect data mpi!: ntasks is not times of nbin"
			stop
		end if
		sent_count=nbin/ntasks
		do i=1, ntasks
			idbg=(i-1)*sent_count+1
			ided=i*sent_count
			sentbuffer(idbg:ided)=fx(ibg:ied)
		end do
		
		call mpi_alltoall(sentbuffer,sent_count, MPI_DOUBLE, reivbuffer, sent_count, &
			MPI_DOUBLE, mpi_comm_world, ierr)

		!do i=1, ntasks
		!	idbg=(i-1)*ntasks+1
		!	ided=i*ntasks
		!	fx(idbg:ided)=reivbuffer(idbg:ided)
		!end do
		fx=reivbuffer
		
	end subroutine
	subroutine collect_data_mpi_y(fxy, nbin, nbg, ned, nblocks,ntasks)
		!use com_main_gw
		implicit none
		!type(diffuse_coeffient_type)::dc
		integer nbin
		real(8) fxy(nbin, nbin)
		real(8),allocatable::reivbuffer(:)
		real(8),allocatable::sentbuffer(:)
		integer nd_tot,i, j, k, idbg, ided, sent_count, ierr
		integer nbg, ned,nblocks, ntasks

		nd_tot=nbin*nbin
		allocate(reivbuffer(nd_tot))
		allocate(sentbuffer(nd_tot))
		sentbuffer=0
		reivbuffer=0
		if(mod(nbin, ntasks).ne.0)then
			print*, "error in collect data mpi!: ntasks is not times of nbin"
			stop
		end if
		sent_count=nbin*nblocks
		!print*, "rid, sent_count, nd_tot=", rid, sent_count, nd_tot
		do i=1, ntasks
			k=0
			do j=nbg, ned, ntasks
				k=k+1
				idbg=(k-1)*nbin+1+sent_count*(i-1)
				ided=idbg+nbin-1
				!write(*,fmt="(A40, 4I6, 20F10.3)") "sent:rid, j, idbg, ided=", rid, j,idbg, ided, fxy(j,1:nbin)
				sentbuffer(idbg:ided)=fxy(j,1:nbin)
			end do
		end do
		!print*, "rid,sentbuffer=", rid,  sentbuffer
		!print*
		call mpi_alltoall(sentbuffer,sent_count, MPI_DOUBLE, reivbuffer, sent_count, &
			MPI_DOUBLE, mpi_comm_world, ierr)
		!print*, "alltoall finished"
		!print*, "rid,recibuffer=", rid,  reivbuffer
		!print*
		!stop
		!do i=1, nblocks
		
		do i=1, ntasks
			k=0
			do j=i,ned, ntasks
				k=k+1
				idbg=(k-1)*nbin+1+sent_count*(i-1)
				ided=idbg+nbin-1
				!write(*,fmt="(A40, 6I6, 20F10.3)") "recv:rid, i,j,k,idbg,ided=", rid, i,j,k, idbg,ided, reivbuffer(idbg:ided)
				fxy(j, 1:nbin)=reivbuffer(idbg:ided)
			end do
		end do
		!call mpi_BARRIER(mpi_comm_world, ierr)
		!print*, "rid=", rid
		!!call dc%s2_de_110%print()

	end subroutine
	subroutine collect_data_mpi_x(fxy, nbin, nbg, ned, nblocks,ntasks)
		!use com_main_gw
		implicit none
		!type(diffuse_coeffient_type)::dc
		integer nbin
		real(8) fxy(nbin, nbin)
		real(8),allocatable::reivbuffer(:)
		real(8),allocatable::sentbuffer(:)
		integer nd_tot,i, j, idbg, ided, sent_count, ierr
		integer nbg, ned,nblocks, ntasks

		nd_tot=nbin*nbin
		allocate(reivbuffer(nd_tot))
		allocate(sentbuffer(nd_tot))
		reivbuffer=0
		if(mod(nbin, ntasks).ne.0)then
			print*, "error in collect data mpi!: ntasks is not times of nbin"
			stop
		end if
		sent_count=nbin*nblocks
		!print*, "rid, sent_count, nd_tot=", rid, sent_count, nd_tot
		do i=1, ntasks
			do j=nbg, ned
				idbg=(i-1)*sent_count+(j-nbg)*nbin+1
				ided=idbg+nbin-1
				!print*, "idbg, ided=", idbg, ided
				sentbuffer(idbg:ided)=fxy(1:nbin,j)
			end do
		end do
		!print*, "rid,sentbuffer=", rid,  sentbuffer
		!print*
		call mpi_alltoall(sentbuffer,sent_count, MPI_DOUBLE, reivbuffer, sent_count, &
			MPI_DOUBLE, mpi_comm_world, ierr)
		!print*, "alltoall finished"
		!print*, "rid,reibuffer=", rid,  reivbuffer
		!print*
		
		do i=1, nbin
			!do j=1, dc%nbin
			idbg=(i-1)*nbin+1
			ided=i*nbin
			!print*, "rid,idbg, ided=", rid, idbg, ided
			fxy(1:nbin,i)=reivbuffer(idbg:ided)
			!end do
		end do
		!call mpi_BARRIER(mpi_comm_world, ierr)
		!print*, "rid=", rid
		!!call dc%s2_de_110%print()

	end subroutine
	subroutine collect_to_root_sps_single(sps_send,sps, n)
		!use com_main_gw
		implicit none
		integer i,n
		type(particle_samples_arr_type)::sps(n)
		type(particle_samples_arr_type) sps_send
		!allocate(sps(ctl%ntasks))
		do i=1, ctl%ntasks
			if(i.ne.mpi_master_id+1)then
				!print*, "rid,i=",rid,i
				call send_particle_sample_arr_mpi(sps_send,sps(i), i-1, mpi_master_id)
			!end if
			else
				sps(mpi_master_id+1)=sps_send
			end if
		end do
	end subroutine 

	subroutine send_particle_sample_arr_mpi(sps_send, sps_recv, proc_id_source,proc_id_dest)
		!use com_main_gw
		implicit none
		type(particle_samples_arr_type)::sps_send,sps_recv
		integer ierr, proc_id_source,proc_id_dest, i, nintarr, nrealarr
		integer status(MPI_Status_size)
		integer,allocatable::intarr(:,:)
		real(8),allocatable::realarr(:,:)
		type(binary),allocatable::by(:), byini(:),bybf(:)

		if(rid.eq.proc_id_source)then
			!print*, "source:rid=",rid,proc_id_dest
			allocate(intarr(nint_particle,sps_send%n))
			allocate(realarr(nreal_particle,sps_send%n))
			allocate(by(sps_send%n),byini(sps_send%n),bybf(sps_send%n))
			do i=1, sps_send%n
				call conv_sp_int_real_arrays(sps_send%sp(i), intarr(1:nint_particle,i), realarr(1:nreal_particle,i))
				by(i)=sps_send%sp(i)%byot
				byini(i)=sps_send%sp(i)%byot_ini
				bybf(i)=sps_send%sp(i)%byot_bf
			end do
			call mpi_send(sps_send%n,1, MPI_INTEGER,proc_id_dest,0,MPI_COMM_WORLD,ierr)
			!print*, "sps_send:n=", sps_send%n, rid
			call mpi_send(intarr, nint_particle*sps_send%n, MPI_INTEGER, proc_id_dest, 0, MPI_COMM_WORLD, ierr)
			!print*, "intarr sent:rid=",rid
			call mpi_send(realarr, nreal_particle*sps_send%n, MPI_DOUBLE_PRECISION, proc_id_dest, 0, MPI_COMM_WORLD, ierr)
			!print*, "realarr sent:rid=",rid,realarr(1:2,1)
            !print*, "by sent:rid=",rid,by(1)%a_bin, proc_id_source, proc_id_dest
			call mpi_send(by(1:sps_send%n), sps_send%n, BINARY_MPI_TYPE_,proc_id_dest,0,MPI_COMM_WORLD,ierr)
			!print*, "by sent:rid=",rid,by(1)%a_bin
			call mpi_send(byini(1:sps_send%n), sps_send%n, BINARY_MPI_TYPE_,proc_id_dest,0,MPI_COMM_WORLD,ierr)
			call mpi_send(bybf(1:sps_send%n), sps_send%n, BINARY_MPI_TYPE_,proc_id_dest,0,MPI_COMM_WORLD,ierr)
			!print*, "bybf send:rid=",rid,bybf(2)%e_bin
		elseif(rid.eq.proc_id_dest)then
			!print*, "dest:rid=",rid,proc_id_source
			call mpi_recv(sps_recv%n,1, MPI_INTEGER,proc_id_source,0,MPI_COMM_WORLD,status,ierr)

			call sps_recv%init(sps_recv%n)
			!print*, "sps_recv:n=", sps_recv%n, allocated(sps_recv%sp),rid
			allocate(intarr(nint_particle,sps_recv%n))
			allocate(realarr(nreal_particle,sps_recv%n))
			allocate(by(sps_recv%n),byini(sps_recv%n),bybf(sps_recv%n))
			
			nintarr=nint_particle*sps_recv%n; nrealarr=nreal_particle*sps_recv%n

			call mpi_recv(intarr, nintarr, MPI_INTEGER, proc_id_source, 0, mpi_comm_world,status, ierr)
			!print*, "intarr recv:rid=",rid
			call mpi_recv(realarr, nrealarr, MPI_DOUBLE_PRECISION, proc_id_source, 0, mpi_comm_world, status,ierr)
			!print*, "realarr recv:rid=",rid,realarr(1:2,1)
			call mpi_recv(by(1:sps_recv%n), sps_recv%n, BINARY_MPI_TYPE_, proc_id_source, 0, mpi_comm_world, status,ierr)
			!print*, "by recv:rid=",rid,by(1)%a_bin
			call mpi_recv(byini(1:sps_recv%n), sps_recv%n, BINARY_MPI_TYPE_, proc_id_source, 0, mpi_comm_world, status,ierr)
			call mpi_recv(bybf(1:sps_recv%n), sps_recv%n, BINARY_MPI_TYPE_, proc_id_source, 0, mpi_comm_world, status,ierr)
			!print*, "bybf recv:rid=",rid,bybf(2)%e_bin
			do i=1, sps_recv%n
				call conv_int_real_arrays_sp(sps_recv%sp(i), intarr(1:nint_particle, i), realarr(1:nreal_particle,i))
				sps_recv%sp(i)%byot=by(i)
				sps_recv%sp(i)%byot_ini=byini(i)
				sps_recv%sp(i)%byot_bf=bybf(i)
				!call sps_recv%sp(i)%print("send_particle_sample_arr_mpi")
			end do

		end if
	end subroutine
 
	  
	subroutine conv_sp_int_real_arrays(sp, intarr, realarr)
		!use com_main_gw
		implicit none
		class(particle_sample_type)::sp
		integer nint, nreal
		integer intarr(nint_particle)
		real(8) realarr(nreal_particle)
		intarr(1:5)=(/sp%id,sp%obtype,sp%obidx,sp%state_flag_last,sp%exit_flag/)
		intarr(6:10)=(/sp%source,sp%within_jt, sp%rid, sp%idx, sp%length/)
		intarr(11:12)=(/sp%write_down_track,sp%N_gene/)
		realarr(1:5)=(/sp%create_time,sp%jph,sp%x,sp%jc,sp%ra/)
		realarr(6:10)=(/sp%En,sp%en0,sp%exit_time,sp%Jm,sp%weight_clone/)
		realarr(11:14)=(/sp%weight_real,sp%m,sp%r_lc,sp%jm0/)
		realarr(15:19)=(/sp%period,sp%rp, sp%tgw, sp%weight_N,sp%simu_bgtime/)	
	end subroutine
	 

	subroutine conv_int_real_arrays_sp(sp, intarr, realarr)
		!use com_main_gw
		implicit none
		class(particle_sample_type)::sp
		integer nint, nreal
		integer intarr(nint_particle)
		real(8) realarr(nreal_particle)
		sp%id=intarr(1); sp%obtype=intarr(2); sp%obidx=intarr(3)
		sp%state_flag_last=intarr(4);sp%exit_flag=intarr(5)
		!intarr(1:5)=(/sp%id,sp%obtype,sp%obidx,sp%state_flag_last,sp%exit_flag/)
		sp%source=intarr(6);sp%within_jt=intarr(7);sp%rid=intarr(8);
		sp%idx=intarr(9)
		sp%length=intarr(10); sp%write_down_track=intarr(11)
		call track_init(sp, 0)
		sp%N_gene=intarr(12)
        
		sp%create_time=realarr(1); sp%jph=realarr(2); sp%x=realarr(3)
		sp%jc=realarr(4);sp%ra=realarr(5)
		
		sp%en=realarr(6);sp%en0=realarr(7);sp%exit_time=realarr(8);
		sp%jm=realarr(9);sp%weight_clone=realarr(10)
		!realarr(6:10)=(/sp%En,sp%en0,sp%exit_time,sp%Jm,sp%weight/)
		!sp%weight_asym=realarr(11);
        sp%weight_real=realarr(11); sp%m=realarr(12)
		sp%r_lc=realarr(13);sp%jm0=realarr(14)
		!realarr(11:15)=(/sp%weight0,sp%weight_real,sp%m,sp%r_lc,sp%jm0/)
		sp%period=realarr(15);sp%rp=realarr(16);sp%tgw=realarr(17)
		sp%weight_N=realarr(18)
        sp%simu_bgtime=realarr(19)
		!call sp%print("conv_int_real_arrays_sp")
		!realarr(16:18)=(/sp%period,sp%rp, sp%tgw/)	
	end subroutine

	 

	subroutine declare_particle_mpi_type()
		!use com_main_gw
		implicit none
		type(particle)::pt
		integer i
        integer,parameter::ptcomp=9
		integer ierr, blocklen(ptcomp), vatype(ptcomp)
		!integer,dimension(ptcomp)::disp=kind(MPI_ADDRESS_KIND)
		!integer::base=kind(MPI_ADDRESS_KIND)
		integer(KIND=MPI_ADDRESS_KIND)::disp(ptcomp), base
		
		
		call mpi_get_address(pt%x, disp(1),ierr)
		call mpi_get_address(pt%vx, disp(2),ierr)
        call mpi_get_address(pt%spin, disp(3),ierr)
		call mpi_get_address(pt%M, disp(4),ierr)
		call mpi_get_address(pt%radius, disp(5),ierr)
		call mpi_get_address(pt%id, disp(6),ierr)
        call mpi_get_address(pt%obtype, disp(7),ierr)
		call mpi_get_address(pt%obidx, disp(8),ierr)
        call mpi_get_address(pt%N_gene, disp(9),ierr)

		base=disp(1)

		do i=1, ptcomp
			disp(i)=disp(i)-base
		end do
		blocklen=1
		blocklen(1)=3
		blocklen(2)=3
		blocklen(3)=3
		
		
		vatype(1:5)=MPI_DOUBLE_PRECISION
        vatype(6:9)=MPI_INTEGER
		!print*, "disp=",disp
		!stop
		call mpi_type_create_struct(ptcomp,blocklen,disp,vatype,PARTICLE_MPI_TYPE_,ierr)
		call mpi_type_commit(PARTICLE_MPI_TYPE_,ierr)

	end subroutine

	subroutine declare_binary_mpi_type()
		!use com_main_gw
		implicit none
		type(binary)::by
		integer i
		integer ierr, blocklen(20), vatype(20)
		!integer::disp(20)=kind(MPI_ADDRESS_KIND)
		!integer::base=kind(MPI_ADDRESS_KIND)
		integer(KIND=MPI_ADDRESS_KIND)::disp(20), base

		call mpi_get_address(by%ms, disp(1),ierr)
		call mpi_get_address(by%mm, disp(2),ierr)
		call mpi_get_address(by%rd, disp(3),ierr)
		call mpi_get_address(by%E, disp(4),ierr)
		call mpi_get_address(by%l, disp(5),ierr)
		call mpi_get_address(by%k, disp(6),ierr)
		call mpi_get_address(by%miu, disp(7),ierr)
		call mpi_get_address(by%mtot, disp(8),ierr)
		call mpi_get_address(by%Jc, disp(9),ierr)
		call mpi_get_address(by%a_bin, disp(10),ierr)
		call mpi_get_address(by%e_bin, disp(11),ierr)
		call mpi_get_address(by%lum, disp(12),ierr)
		call mpi_get_address(by%f0, disp(13),ierr)
		call mpi_get_address(by%Inc, disp(14),ierr)
		call mpi_get_address(by%Om, disp(15),ierr)	
		call mpi_get_address(by%Pe, disp(16),ierr)
		call mpi_get_address(by%t0, disp(17),ierr)
		call mpi_get_address(by%me, disp(18),ierr)
		call mpi_get_address(by%an_in_mode, disp(19),ierr)	
		call mpi_get_address(by%bname, disp(20),ierr)	

		base=disp(1)
		do i=1, 20
			disp(i)=disp(i)-base
		end do
		blocklen=1
		blocklen(12)=3
		blocklen(20)=100
		vatype(1:3)=PARTICLE_MPI_TYPE_
		vatype(4:18)=MPI_DOUBLE_PRECISION
		vatype(19)=MPI_INTEGER
		vatype(20)=MPI_CHARACTER

		call mpi_type_create_struct(20,blocklen,disp,vatype,BINARY_MPI_TYPE_,ierr)
		call mpi_type_commit(BINARY_MPI_TYPE_,ierr)
		
	end subroutine
  
	subroutine get_dms(dm)
		!use com_main_gw
		implicit none
		type(diffuse_mspec)::dm
		integer ierr
        call mpi_BARRIER(mpi_comm_world, ierr)

        print*, "start get diffuse coefficients", rid	
		!call dms%mb(1)%all%gxjcr%print("all%gxjjcr")
		!stop
!=================================		
        call dm_get_dc_mpi(dm)
        print*, "cal dms finished",rid
        call mpi_BARRIER(mpi_comm_world, ierr)

	end subroutine
	  

subroutine collection_fxy_int(fxy,nx,ny)
    use com_sts_type
    !class(bin_function_type)::s1d
    integer ierr, nx,ny, i,j
	integer fxy(nx,ny)
    integer sentbuffer(nx,ny,ctl%ntasks)
    integer recvbuffer(nx,ny,ctl%ntasks)

    do i=1, ctl%ntasks
        sentbuffer(1:nx,1:ny,i)=fxy(1:nx,1:ny)
    end do
    recvbuffer=0

    call mpi_alltoall(sentbuffer,nx*ny,MPI_INTEGER,recvbuffer,&
       nx*ny, MPI_INTEGER, mpi_comm_world,ierr)

    do i=1, nx
		do j=1, ny
        	fxy(i,j)=sum(recvbuffer(i,j,:))
		end do
    end do
end subroutine

subroutine collection_fx_int(fx,nx)
    use com_sts_type
    !class(bin_function_type)::s1d
    integer ierr, nx, i
	integer fx(nx)
    integer sentbuffer(nx,ctl%ntasks)
    integer recvbuffer(nx,ctl%ntasks)

    do i=1, ctl%ntasks
        sentbuffer(1:nx,i)=fx(1:nx)
    end do
    recvbuffer=0

    call mpi_alltoall(sentbuffer,nx,MPI_INTEGER,recvbuffer,&
       nx, MPI_INTEGER, mpi_comm_world,ierr)

    do i=1, nx
        fx(i)=sum(recvbuffer(i,:))
    end do
end subroutine

subroutine collection_and_avg_fx(fx,nx)
    use com_sts_type
    !class(bin_function_type)::s1d
    integer ierr, nx, ntasks, i
	real(8) fx(nx)
    real(8) sentbuffer(nx,ctl%ntasks)
    real(8) recvbuffer(nx,ctl%ntasks)

    do i=1, ctl%ntasks
        sentbuffer(1:nx,i)=fx(1:nx)
    end do
    recvbuffer=0

    call mpi_alltoall(sentbuffer,nx,MPI_double,recvbuffer,&
       nx, MPI_double, mpi_comm_world,ierr)

	do i=1, nx
		fx(i)=sum(recvbuffer(i,1:ctl%ntasks))/dble(sample_threads_number)
	end do

end subroutine

subroutine collection_and_max_fx(fx,nx)
    use com_sts_type
    !class(bin_function_type)::s1d
    integer ierr, nx, ntasks, i
	real(8) fx(nx)
    real(8) sentbuffer(nx,ctl%ntasks)
    real(8) recvbuffer(nx,ctl%ntasks)

    do i=1, ctl%ntasks
        sentbuffer(1:nx,i)=fx(1:nx)
    end do
    recvbuffer=0

    call mpi_alltoall(sentbuffer,nx,MPI_double,recvbuffer,&
       nx, MPI_double, mpi_comm_world,ierr)

	do i=1, nx
		fx(i)=maxval(recvbuffer(i,1:ctl%ntasks))
	end do

end subroutine

subroutine collection_and_min_fx(fx,nx)
    use com_sts_type
    !class(bin_function_type)::s1d
    integer ierr, nx, ntasks, i
	real(8) fx(nx)
    real(8) sentbuffer(nx,ctl%ntasks)
    real(8) recvbuffer(nx,ctl%ntasks)

    do i=1, ctl%ntasks
        sentbuffer(1:nx,i)=fx(1:nx)
    end do
    recvbuffer=0

    call mpi_alltoall(sentbuffer,nx,MPI_double,recvbuffer,&
       nx, MPI_double, mpi_comm_world,ierr)

	do i=1, nx
		fx(i)=minval(recvbuffer(i,1:ctl%ntasks))
	end do

end subroutine

subroutine collection_int(nx)
    use com_sts_type
    !class(bin_function_type)::s1d
    integer ierr,  i
	integer nx
    integer sentbuffer(ctl%ntasks)
	integer recvbuffer(ctl%ntasks)

    sentbuffer=nx
    recvbuffer=0

    call mpi_alltoall(sentbuffer,1,MPI_integer,recvbuffer,&
       1, MPI_integer, mpi_comm_world,ierr)

	nx=sum(recvbuffer(:))
end subroutine
subroutine collection_and_avg_real(nx)
    use com_sts_type
	implicit none
    !class(bin_function_type)::s1d
    integer ierr, i
	real(8) nx
    real(8) sentbuffer(ctl%ntasks)
	real(8) recvbuffer(ctl%ntasks)

	sentbuffer=nx
    recvbuffer=0

    call mpi_alltoall(sentbuffer,1,MPI_double,recvbuffer,&
       1, MPI_double, mpi_comm_world,ierr)
	nx=sum(recvbuffer(:))/dble(sample_threads_number)
end subroutine

subroutine collection_and_get_max_real(nx)
    use com_sts_type
	implicit none
    !class(bin_function_type)::s1d
    integer ierr, i
	real(8) nx
    real(8) sentbuffer(ctl%ntasks)
	real(8) recvbuffer(ctl%ntasks)

	sentbuffer=nx
    recvbuffer=0

    call mpi_alltoall(sentbuffer,1,MPI_double,recvbuffer,&
       1, MPI_double, mpi_comm_world,ierr)
	nx=maxval(recvbuffer(:))
end subroutine

subroutine collection_and_get_min_real(nx)
    use com_sts_type
	implicit none
    !class(bin_function_type)::s1d
    integer ierr, i
	real(8) nx
    real(8) sentbuffer(ctl%ntasks)
	real(8) recvbuffer(ctl%ntasks)

	sentbuffer=nx
    recvbuffer=0

    call mpi_alltoall(sentbuffer,1,MPI_double,recvbuffer,&
       1, MPI_double, mpi_comm_world,ierr)
	nx=minval(recvbuffer(:))
end subroutine

subroutine collection_and_avg_s2d(fxy, nx, ny)
    !use mpi_comu
    !use com_sts_type
    !type(s2d_type)::s2d
    integer ierr, nx, ny,  i, j
    real(8) sentbuffer(nx,ny,ctl%ntasks)
    real(8) recvbuffer(nx,ny,ctl%ntasks)
    real(8) fxy(nx, ny)

    !nx=s2d%nx
    !ny=s2d%ny
    do i=1, ctl%ntasks
        sentbuffer(1:nx,1:ny, i)=fxy(1:nx,1:ny)
        !sentbuffer(1:nx,1:ny,2, i)=s2d%fxyw(1:nx,1:ny)
    end do
    recvbuffer=0

    call mpi_alltoall(sentbuffer,nx*ny,MPI_double,recvbuffer,&
       nx*ny, MPI_double, mpi_comm_world,ierr)
    
	do i=1, nx
		do j=1, ny
			fxy(i,j)=sum(recvbuffer(i,j,:))/dble(sample_threads_number)
		end do
	end do
	
end subroutine

end module
