subroutine get_dstr_num_in_each_ird_bin(x, n, xb, xstep, nbin, fx,  n_num)
	implicit none
	integer i,j, n, nbin
	real(8) x(n), w(n), xb(nbin), xstep(nbin)
	integer fx(nbin)
	integer n_num
	integer indx
	fx=0
	n_num=0
	do i=1, n
		do j=1, nbin
			if(xb(j)-xstep(j)/2d0<x(i).and.x(i)<xb(j)+xstep(j)/2d0)then
				indx=j
				fx(indx)=fx(indx)+1
				n_num=n_num+1
			end if
		end do
	end do
end subroutine
subroutine get_dstr_num_in_each_ird_bin_weight(x, w, n, xb, xstep, nbin, fxw, n_numw)
	implicit none
	integer i,j, n, nbin
	real(8) x(n), w(n), xb(nbin), xstep(nbin)
	real(8) fxw(nbin),  n_numw
	!print*, "n, m=", n, nbin
	fxw=0
	n_numw=0
	do i=1, n
		do j=1, nbin
			!print*, "xb(j),xstep(j),x=", xb(j),xstep(j), x(i)
			if(xb(j)-xstep(j)/2d0<x(i).and.x(i)<xb(j)+xstep(j)/2d0)then
				fxw(j)=fxw(j)+w(i)
				n_numw=n_numw+w(i)
			end if
		end do
	end do
end subroutine
subroutine get_dstr_num_in_each_ird_bin_normal_kernel(x, n, xmin,xmax,xb, xstep, nbin, fx,  n_num)
	implicit none
	integer i,j, n, nbin
	real(8) x(n), w(n), xb(nbin), xstep(nbin)
	integer fx(nbin)
	integer n_num
	real(8) vx,vx1,vx2,xmin,xmax
	real(8),parameter::pi=3.14159265
	integer indx
	fx=0
	n_num=0
	do i=1, n
		do j=1, nbin
			vx=(x(i)-xb(j))/xstep(i)
			vx2=(xmax*2-xb(j)-x(i))/xstep(i)
			vx1=(xb(j)+x(i)-2*xmin)/xstep(i)
			!	vy=(y(k)-yb(j))/ystep(j)
			fx(j)=fx(j)+1d0/(2d0*pi)**0.5*exp(-0.5d0*vx**2)+1d0/(2d0*pi)**0.5*exp(-0.5d0*min(vx1,vx2)**2)
		end do
		n_num=n_num+1
	end do
end subroutine
subroutine get_dstr_num_in_each_ird_bin_delta_kernel(x, n, xmin,xmax,xb, xstep, nbin, fx,  n_num)
	implicit none
	integer i,j, n, nbin
	real(8) x(n), w(n), xb(nbin), xstep(nbin)
	integer fx(nbin)
	integer n_num
	real(8) vx,vx1,vx2,xmin,xmax
	real(8),parameter::pi=3.14159265
	integer indx
	fx=0
	n_num=0
	do i=1, n
		do j=1, nbin
			vx=(x(i)-xb(j))/xstep(j)
			vx2=(xmax*2-xb(j)-x(i))/xstep(j)
			vx1=(xb(j)+x(i)-2*xmin)/xstep(j)
			!	vy=(y(k)-yb(j))/ystep(j)
			if(vx>-0.5.and.vx<0.5)then
				fx(j)=fx(j)+1
			end if
			if(vx1>-0.5.and.vx1<0.5)then
				fx(j)=fx(j)+1
			end if
			if(vx2>-0.5.and.vx2<0.5)then
				fx(j)=fx(j)+1
			end if
		end do
		n_num=n_num+1
	end do
end subroutine
subroutine get_dstr_num_in_each_ird_bin_weight_normal_kernel(x, w, n,xmin,xmax, xb, xstep, nbin, fxw, n_numw)
	implicit none
	integer i,j, n, nbin
	real(8) x(n), w(n), xb(nbin), xstep(nbin)
	real(8) fxw(nbin),  n_numw
	real(8) vx,vx1,vx2,xmin,xmax
	real(8),parameter::pi=3.14159265
	!print*, "n, m=", n, nbin
	fxw=0
	n_numw=0
	do i=1, n
		do j=1, nbin
			!print*, "xb(j),xstep(j),x=", xb(j),xstep(j), x(i)
			!vx=(x(i)-xb(j))/xstep(j)
			!fxw(j)=fxw(j)+w(i)/(2d0*pi)**0.5*exp(-0.5d0*vx**2)
			vx=(x(i)-xb(j))/xstep(j)
			vx2=(xmax*2-xb(j)-x(i))/xstep(j)
			vx1=(xb(j)+x(i)-2*xmin)/xstep(j)
			!	vy=(y(k)-yb(j))/ystep(j)
			fxw(j)=fxw(j)+w(i)/(2d0*pi)**0.5*exp(-0.5d0*vx**2)+w(i)/(2d0*pi)**0.5*exp(-0.5d0*min(vx1,vx2)**2)
		end do
		n_numw=n_numw+w(i)
	end do
end subroutine
subroutine get_dstr_num_in_each_ird_bin_weight_delta_kernel(x, w, n,xmin,xmax, xb, xstep, nbin, fxw, n_numw)
	implicit none
	integer i,j, n, nbin
	real(8) x(n), w(n), xb(nbin), xstep(nbin)
	real(8) fxw(nbin),  n_numw
	real(8) vx,vx1,vx2,xmin,xmax
	real(8),parameter::pi=3.14159265
	!print*, "n, m=", n, nbin
	fxw=0
	n_numw=0
	do i=1, n
		do j=1, nbin
			!print*, "xb(j),xstep(j),x=", xb(j),xstep(j), x(i)
			!vx=(x(i)-xb(j))/xstep(j)
			!fxw(j)=fxw(j)+w(i)/(2d0*pi)**0.5*exp(-0.5d0*vx**2)
			vx=(x(i)-xb(j))/xstep(j)
			vx2=(xmax*2-xb(j)-x(i))/xstep(j)
			vx1=(xb(j)+x(i)-2*xmin)/xstep(j)
			if(vx>-0.5.and.vx<0.5)then
				fxw(j)=fxw(j)+1
			end if
			if(vx1>-0.5.and.vx1<0.5)then
				fxw(j)=fxw(j)+1
			end if
			if(vx2>-0.5.and.vx2<0.5)then
				fxw(j)=fxw(j)+1
			end if
		end do
		n_numw=n_numw+w(i)
	end do
end subroutine

subroutine get_dstr_num_in_each_ird_bin_weight_s2d(x,y,w, n, xb,xstep, nx, yb, ystep, ny, fxyw,  n_numw)
	implicit none
	integer i,j,k, n, nx,ny,n_numw, idx, idy
	real(8) x(n), y(n), w(n),xb(nx), yb(ny), xstep(nx), ystep(ny)
	real(8) fxyw(nx,ny)!, nxyw(nx, ny)
	real(8) xmin,xmax,ymin,ymax
	fxyw=0
	n_numw=0
	
	xmin=xb(1)-xstep(1)/2
	xmax=xb(nx)+xstep(nx)/2
	ymin=yb(1)-ystep(1)/2
	ymax=yb(ny)+ystep(ny)/2
	do k=1, n
		if(x(k)<=xmax.and.x(k)>=xmin.and.y(k)>=ymin.and.y(k)<=ymax)then
			call return_idx_ir_dstr(xb,xstep,nx,x(k),idx,1)
			
			!print*, "xb=", xb
			!print*, "x=",x(k), idx

			call return_idx_ir_dstr(yb,ystep,ny,y(k),idy,1)
			!print*, "yb=", yb
			!print*, "y=",y(k), idy
			!read(*,*)
			fxyw(idx,idy)=fxyw(idx,idy)+w(k)
			n_numw=n_numw+w(k)
		end if
	end do
end subroutine

subroutine get_dstr_num_in_each_ird_bin_weight_s2d_delta_kernel(x,y,w, n,xmin,xmax, xb,xstep, nx, yb, ystep, ny, fxyw,  n_numw)
	implicit none
	integer i,j,k, n, nx,ny,n_numw, idx, idy
	real(8) x(n), y(n), w(n),xb(nx), yb(ny), xstep(nx), ystep(ny)
	real(8) fxyw(nx,ny)!, nxyw(nx, ny)
	real(8) xmin,xmax
	real(8) vx,vy,vx2,vx1
	real(8),parameter::pi=3.14159265

	fxyw=0
	n_numw=0
	do k=1, n
		call return_idx_ir_dstr(yb,ystep,ny,y(k),idy,1)
		do i=1, nx			
			!if((y(k)<=yb(idy)+ystep(idy)/2d0).and.(y(k)>=yb(idy)-ystep(idy)/2d0))then
				vx=(x(k)-xb(i))/xstep(i)
				vx2=(xmax*2-xb(i)-x(k))/xstep(i)
				vx1=(xb(i)+x(k)-2*xmin)/xstep(i)
				if(vx>-0.5.and.vx<0.5)then
					fxyw(i,idy)=fxyw(i,idy)+w(k)
				end if
				if(vx1>-0.5.and.vx1<0.5)then
					fxyw(i,idy)=fxyw(i,idy)+w(k)
				end if
				if(vx2>-0.5.and.vx2<0.5)then
					fxyw(i,idy)=fxyw(i,idy)+w(k)
				end if
			!end if
			!end do
		end do
		n_numw=n_numw+w(k)
	end do
end subroutine
subroutine get_dstr_num_in_each_ird_bin_s2d_delta_kernel(x,y, n,xmin,xmax, xb,xstep, nx, yb, ystep, ny, fxy,  n_num)
	implicit none
	integer i,j,k, n, nx,ny,n_num,idx,idy
	real(8) x(n),y(n), xb(nx), yb(ny), xstep(nx), ystep(ny)
	integer fxy(nx,ny)!, nxy(nx, ny)
	real(8) vx, vy,vx1,vx2,xmin,xmax
	real(8),parameter::pi=3.14159265
	
	fxy=0
	n_num=0
	do k=1, n
		do i=1, nx
			call return_idx_ir_dstr(yb,ystep,ny,y(k),idy,1)
			!if((y(k)<=yb(idy)+ystep(idy)/2d0).and.(y(k)>=yb(idy)+ystep(idy)/2d0))then
				vx=(x(k)-xb(i))/xstep(i)
				vx2=(xmax*2-xb(i)-x(k))/xstep(i)
				vx1=(xb(i)+x(k)-2*xmin)/xstep(i)
				if(vx>-0.5.and.vx<0.5)then
					fxy(i,idy)=fxy(i,idy)+1
				end if
				if(vx1>-0.5.and.vx1<0.5)then
					fxy(i,idy)=fxy(i,idy)+1
				end if
				if(vx2>-0.5.and.vx2<0.5)then
					fxy(i,idy)=fxy(i,idy)+1
				end if
			!end if
			!end do
		end do
		n_num=n_num+1
	end do

end subroutine


subroutine get_dstr_num_in_each_ird_bin_weight_s2d_normal_kernel(x,y,w, n,xmin,xmax, xb,xstep, nx, yb, ystep, ny, fxyw,  n_numw)
	implicit none
	integer i,j,k, n, nx,ny,n_numw, idx, idy
	real(8) x(n), y(n), w(n),xb(nx), yb(ny), xstep(nx), ystep(ny)
	real(8) fxyw(nx,ny)!, nxyw(nx, ny)
	real(8) xmin,xmax
	real(8) vx,vy,vx2,vx1
	real(8),parameter::pi=3.14159265

	fxyw=0
	n_numw=0
	do k=1, n
		call return_idx_ir_dstr(yb,ystep,ny,y(k),idy,1)
		do i=1, nx			
			!if((y(k)<=yb(idy)+ystep(idy)/2d0).and.(y(k)>=yb(idy)+ystep(idy)/2d0))then
				vx=(x(k)-xb(i))/xstep(i)
				vx2=(xmax*2-xb(i)-x(k))/xstep(i)
				vx1=(xb(i)+x(k)-2*xmin)/xstep(i)
				!	vy=(y(k)-yb(j))/ystep(j)
				fxyw(i,idy)=fxyw(i,idy)+w(k)/(2d0*pi)**0.5*exp(-0.5d0*vx**2)+w(k)/(2d0*pi)**0.5*exp(-0.5d0*min(vx1,vx2)**2)
			!end if
			!end do
		end do
		n_numw=n_numw+w(k)
	end do
end subroutine
subroutine get_dstr_num_in_each_ird_bin_s2d_normal_kernel(x,y, n,xmin,xmax, xb,xstep, nx, yb, ystep, ny, fxy,  n_num)
	implicit none
	integer i,j,k, n, nx,ny,n_num,idx,idy
	real(8) x(n),y(n), xb(nx), yb(ny), xstep(nx), ystep(ny)
	integer fxy(nx,ny)!, nxy(nx, ny)
	real(8) vx, vy,vx1,vx2,xmin,xmax
	real(8),parameter::pi=3.14159265
	
	fxy=0
	n_num=0
	do k=1, n
		do i=1, nx
			!do j=1, ny
				call return_idx_ir_dstr(yb,ystep,ny,y(k),idy,1)
				!vx=(x(k)-xb(i))/xstep(i)
				!vy=(y(k)-yb(j))/ystep(j)
				!fxy(i,j)=fxy(i,j)+1d0/2d0/pi*exp(-0.5d0*(vx**2+vy**2))
				!if((y(k)<=yb(idy)+ystep(idy)/2d0).and.(y(k)>=yb(idy)+ystep(idy)/2d0))then
					vx=(x(k)-xb(i))/xstep(i)
					vx2=(xmax*2-xb(i)-x(k))/xstep(i)
					vx1=(xb(i)+x(k)-2*xmin)/xstep(i)
					fxy(i,idy)=fxy(i,idy)+1d0/(2d0*pi)**0.5*exp(-0.5d0*vx**2)+1d0/(2d0*pi)**0.5*exp(-0.5d0*min(vx1,vx2)**2)
				!end if
			!end do
		end do
		n_num=n_num+1
	end do

end subroutine


subroutine get_dstr_num_in_each_ird_bin_s2d(x,y, n, xb,xstep, nx, yb, ystep, ny, fxy,  n_num)
	implicit none
	integer i,j,k, n, nx,ny,n_num,idx,idy
	real(8) x(n),y(n), xb(nx), yb(ny), xstep(nx), ystep(ny)
	integer fxy(nx,ny)!, nxy(nx, ny)
	
	fxy=0
	n_num=0
	do k=1, n
		call return_idx_ir_dstr(xb,xstep,nx,x(k),idx,1)
		call return_idx_ir_dstr(yb,ystep,ny,y(k),idy,1)
		fxy(idx,idy)=fxy(idx,idy)+1
		n_num=n_num+1
	end do

end subroutine



subroutine get_dstr_num_in_each_ird_bin_weight_s2d_direct(x,y,w, n, xb,xstep, nx, yb, ystep, ny, fxyw,  n_numw)
	implicit none
	integer i,j,k, n, nx,ny,n_numw, idx, idy
	real(8) x(n), y(n), w(n),xb(nx), yb(ny), xstep(nx), ystep(ny)
	integer fxyw(nx,ny)!, nxyw(nx, ny)
	
	fxyw=0
	n_numw=0
	do k=1, n
		do i=1, nx
			if(xb(i)-xstep(i)/2d0<x(k).and.x(k)<xb(i)+xstep(i)/2d0)then
				idx=i
			end if
		end do
		do j=1, ny
			if(yb(j)-ystep(j)/2d0<x(k).and.x(k)<yb(j)+ystep(j)/2d0)then
				idy=j
			end if
		end do
		fxyw(idx,idy)=fxyw(idx,idy)+w(k)
		n_numw=n_numw+w(k)
	end do
end subroutine

subroutine get_dstr_num_in_each_ird_bin_s2d_direct(x, n, xb,xstep, nx, yb, ystep, ny, fxy,  n_num)
	implicit none
	integer i,j,k, n, nx,ny,n_num,idx,idy
	real(8) x(n), xb(nx), yb(ny), xstep(nx), ystep(ny)
	integer fxy(nx,ny)!, nxy(nx, ny)
	
	fxy=0
	n_num=0
	do k=1, n
		do i=1, nx
			if(xb(i)-xstep(i)/2d0<x(k).and.x(k)<xb(i)+xstep(i)/2d0)then
				idx=i
			end if
		end do
		do j=1, ny
			if(yb(j)-ystep(j)/2d0<x(k).and.x(k)<yb(j)+ystep(j)/2d0)then
				idy=j
			end if
		end do
		fxy(idx,idy)=fxy(idx,idy)+1
		n_num=n_num+1
	end do

end subroutine