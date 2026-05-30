#! /usr/bin/env python
import numpy as np
import myhdf5_funcs
import pandas as pd
from scipy import interpolate
import myfunc


def _weighted_variance(values, weights):
    """
    计算加权方差
    :param values: 数据数组
    :param weights: 权重数组
    :return: 加权方差
    """
    weighted_avg = np.average(values, weights=weights)
    variance = np.sum(weights * (values - weighted_avg)**2) / np.sum(weights)
    return variance**0.5

class Object:
    def __init__(self,nl):
        self.t=np.zeros((nl))
        self.dt=np.zeros((nl))
        self.logt=np.zeros((nl))
        self.size=nl
        self.rate=np.zeros((nl))
        self.nw=np.zeros((nl))
        self.dmbh=np.zeros((nl))

    def smooth_rate(self,np):
        out=Object(self.size)
        out.rate=pd.Series(self.rate).rolling(window=np,center=True).mean().values
        return out.rate
    def select(self,mask):
        out=Object(len(mask[mask==True]))
        out.t=self.t[mask]
        out.dt=self.dt[mask]
        out.logt=self.logt[mask]
        out.rate=self.rate[mask]
        out.nw=self.nw[mask]
        out.dmbh=self.dmbh[mask]
        return out
    def get_avgsct(self,dt):
        out=Object(2)

        mask=~np.isnan(self.rate)
        if(len(mask)==0):
            out.rate[0]=np.nan
            out.t[0]=np.nan
            out.dt[0]=np.nan
            out.logt[0]=np.nan
            return out
        out.rate[0]=np.average(self.rate[mask],weights=dt[mask])
        mask=~np.isnan(self.rate)
        out.t[0]=np.average(self.t[mask],weights=dt[mask])
        out.dt[0]=np.average(self.dt[mask],weights=dt[mask])
        out.logt[0]=np.average(self.logt[mask],weights=dt[mask])

        mask=~np.isnan(self.nw)
        out.nw[0]=np.average(self.nw[mask],weights=dt[mask])
        mask=~np.isnan(self.dmbh)
        out.dmbh[0]=np.average(self.dmbh[mask],weights=dt[mask])
        mask=~np.isnan(self.rate)
        out.rate[1]=_weighted_variance(self.rate[mask],weights=dt[mask])
        mask=~np.isnan(self.nw)
        out.nw[1]=_weighted_variance(self.nw[mask],weights=dt[mask])
        mask=~np.isnan(self.dmbh)
        out.dmbh[1]=_weighted_variance(self.dmbh[mask],weights=dt[mask])
        return out
    def smooth(self,npoint,method='mov'):

        
        if(method == 'mov'):
            out=Object(self.size)
            out.rate=pd.Series(self.rate).rolling(window=npoint,center=True).mean().values
            out.t=self.t
            out.dt=self.dt
            out.logt=self.logt
        elif(method == 'mov_nozero'):
            mask=self.rate!=0
            out=Object(len(self.t[mask]))
            out.t=self.t[mask]
            out.dt=self.dt[mask]
            out.logt=self.logt[mask]
            out.rate=pd.Series(self.rate[mask]).rolling(window=npoint,center=True).mean().values
            
        # print(self.rate)
        # print(out.rate)
        # print("   --------- ")
        return out
class S1d:
    def __init__(self):
        return
    def init(self,n):
       self.n=n
       self.x=np.zeros(n)
       self.y=np.zeros(n)
       self.yprepared=False
       self.xprepared=False
       self.inter_method='cubic'
       return
    def print(self):
        for i in range(self.n):
            print("i,x,y=",i, self.x[i],self.y[i]) 
    def prepare_fy(self):
        self.fy=interpolate.interp1d(self.x,self.y,kind=self.inter_method)
        self.yprepared=True
    def prepare_fx(self):
        self.fx=interpolate.interp1d(self.y,self.x,kind=self.inter_method)
        self.xprepared=True

    def get_y(self,x):
        if(not self.yprepared):
            self.prepare_fy()
        return self.fy(x)
    def get_x(self,y):
        if(not self.xprepared):
            self.prepare_fx()
        # print("y=",y)    
        return self.fx(y)   
    def readhdf5(self,fl,fsub,print_err=True):
        # print(fl, fsub)
        x,y=myhdf5_funcs.read_2d_hdf5(fl,fsub,print_err=print_err)
        self.init(len(x))
        self.x=x
        self.y=y


class Object_rates:
    
    def __init__(self,nl):
        self.typelist=["ls", "td", "emris", "emax"]
        self.size=nl
        self.alpha_mean=0
        self.fmden=S1d()
        self.fna=S1d()
        self.fma=S1d()
        self.d={}
        
        for tl in self.typelist:
            self.d[tl]=Object(nl)
    
    def readhdf5(self,fl,str_obj,str_type,i,print_err=True,_prex=''):
        da=myhdf5_funcs.read_attr_hdf5(fl,"oe_"+str_obj,str_type,print_err=print_err)
        
        if(da is not None):
            self.d[str_type].rate[i]=da[2]/1e6
            self.d[str_type].nw[i]=da[1]
        
    # def readhdf5den(self,fl,fsub,print_err=True):
    #     # print(fl, fsub)
    #     x,y=myhdf5_funcs.read_2d_hdf5(fl,fsub+"/fmden",print_err=print_err)
    #     self.fmden.init(len(x))
    #     self.fmden.x=x
    #     self.fmden.y=y

    def select(self,mask):
        out=Object_rates(len(mask[mask==True]))
        for tl in self.typelist:
            out.d[tl]=self.d[tl].select(mask)
            
        return out
    def avgsct(self,dt):
        out=Object_rates(2)
        for tl in self.typelist:
            out.d[tl]=self.d[tl].get_avgsct(dt)
        return out
    def smooth(self,np,method='mov'):
        out=Object_rates(self.size)
        for tl in self.typelist:
            out.d[tl]=self.d[tl].smooth(np,method=method)
        return out
        
class Simu_rates:
    
    def __init__(self,nl):
        self.objlist=["star","sbh","ns", "wd", "bd", "rg","all"]
        self.size=nl
        self.r0=0
        self.reff0=0
        self.rh0=0
        self.m0=0
        self.mbh0=0
        
        self.tnr=np.zeros(nl)
        self.t=np.zeros(nl)
        self.dt=np.zeros(nl)
        self.logt=np.zeros(nl)
        self.mbh=np.zeros(nl)
        self.logmbh=np.zeros(nl)
        self.reff=np.zeros(nl)
        self.gas_res=np.zeros(nl)
        self.mcl=np.zeros(nl)
        self.reff=np.zeros(nl)
        self.rh=np.zeros(nl)
        # self.dt=np.zeros(nl)
        self.mass_loss=np.zeros(nl)
        self.td_accum=np.zeros(nl)
        self.direct_swallow=np.zeros(nl)
        self.emri_accum=np.zeros(nl)
        self.emax_accum=np.zeros(nl)
        self.obj={}
        for ol in self.objlist:
            self.obj[ol]=Object_rates(nl)
    def get_avg_data_sets_given_bin(self, mbins,obtype,evtype,das='mbh'):
        num_points=len(mbins)
        x=np.zeros(num_points-1)
        xe=np.zeros(num_points-1)
        ya=np.zeros(num_points-1)
        ys=np.zeros(num_points-1)
        
        for i in range(num_points-1):
             sel=self.select(mbins[i],mbins[i+1],da=das).get_avgsct()
             x[i]=sel.mbh[0]
             xe[i]=sel.mbh[1]
             ya[i]=sel.obj[obtype].d[evtype].rate[0]
             ys[i]=sel.obj[obtype].d[evtype].rate[1]
        return x,xe,ya,ys
    
    def get_avg_data_sets(self, num_points,obtype,evtype,das='mbh'):
        mbins=np.logspace(np.log10(self.mbh[0]),np.log10(self.mbh[-1]),num_points+1)
        print("mbins=",mbins)
        return self.get_avg_data_sets_given_bin( mbins,obtype,evtype,das=das)
    
    def select_by_mask(self,mask):
        out=Simu_rates(len(mask[mask==True]))
        out.t=self.t[mask]
        out.dt=self.dt[mask]
        # print("dt=",out.dt)
        out.logt=self.logt[mask]
        out.mbh=self.mbh[mask]
        out.logmbh=self.logmbh[mask]
        out.reff=self.reff[mask]
        out.rh=self.rh[mask]
        out.gas_res=self.gas_res[mask]
        out.mass_loss=self.mass_loss[mask]
        out.direct_swallow=self.direct_swallow[mask]
        out.td_accum=self.td_accum[mask]
        out.emri_accum=self.emri_accum[mask]
        out.emax_accum=self.emax_accum[mask]
        out.mcl=self.mcl[mask]
        out.obj={}
        for ol in self.objlist:
            out.obj[ol]=self.obj[ol].select(mask)         
        return out    
    
    def select(self,xi,xf,da='time'):
        if(da=='time'):
            mask=np.logical_and(self.t>=xi,self.t<=xf)
        elif(da=='mbh'):
            mask=np.logical_and(self.mbh>=xi,self.mbh<=xf)
        # print("mask=",mask)
        # print(self.size,len(mask[mask==True]))
        return self.select_by_mask(mask)
    
    def get_avgsct(self):
        out=Simu_rates(2)
        out.obj={}
        if(len(self.dt)==0):
             out.mbh[0]=np.nan
             out.mbh[1]=np.nan
             out.t[0]=np.nan
             for ol in self.objlist:
                out.obj[ol]=Object_rates(2)
             return out
        out.t[0]=np.average(self.t,weights=self.dt)
        out.logt[0]=np.log10(out.t[0])
        out.mbh[0]=np.average(self.mbh,weights=self.dt)
        out.logmbh[0]=np.average(self.logmbh,weights=self.dt)
        out.gas_res[0]=np.average(self.gas_res,weights=self.dt)
        out.mcl[0]=np.average(self.mcl,weights=self.dt)
        out.rh[0]=np.average(self.rh,weights=self.dt)
        out.reff[0]=np.average(self.reff,weights=self.dt)
        out.mass_loss[0]=np.average(self.mass_loss,weights=self.dt)
        out.td_accum[0]=np.average(self.td_accum,weights=self.dt)
        out.direct_swallow[0]=np.average(self.direct_swallow,weights=self.dt)
        out.emri_accum[0]=np.average(self.emri_accum,weights=self.dt)
        out.emax_accum[0]=np.average(self.emax_accum,weights=self.dt)
        out.mbh[1]=_weighted_variance(self.mbh,self.dt)
        out.logmbh[1]=_weighted_variance(self.logmbh,self.dt)
        out.gas_res[1]=_weighted_variance(self.gas_res,self.dt)
        out.mcl[1]=_weighted_variance(self.mcl,self.dt)

        for ol in self.objlist:
            if (self.obj[ol] is not None):
                out.obj[ol]=self.obj[ol].avgsct(self.dt)
            else:
                print(ol+"not exist")
                out.obj[ol]=Object_rates(2)
        return out


    def get_smooth(self,np,method='mov'):
        out=Simu_rates(self.size)
        out.obj={}
        #out.objlist=[]

        out.t=self.t
        out.dt=self.dt
        out.logt=self.logt
        out.mbh=self.mbh
        out.mbh0=self.mbh0
        out.logmbh=self.logmbh
        out.gas_res=self.gas_res
        out.mcl=self.mcl
        out.reff=self.reff
        out.rh=self.rh
        out.mass_loss=self.mass_loss
        out.direct_swallow=self.direct_swallow
        out.td_accum=self.td_accum
        out.emri_accum=self.emri_accum
        out.emax_accum=self.emax_accum
        for ol in self.objlist:
            if (self.obj[ol] is not None):
                out.obj[ol]=self.obj[ol].smooth(np,method=method)
            else:
                print(ol+"not exist")
                out.obj[ol]=Object_rates(self.size)

        return out

def get_one_rates(fdir,nl,print_err=True):
    sr=Simu_rates(nl)
    sr.r0=myhdf5_funcs.read_attr_hdf5(fdir+"./output/ini/hdf5/dms_0.hdf5",".","r0(pc)")
    sr.m0=myhdf5_funcs.read_attr_hdf5(fdir+"./output/ini/hdf5/dms_0.hdf5",".","m0(msun)")
    sr.reff0=myhdf5_funcs.read_attr_hdf5(fdir+"./output/ini/hdf5/dms_0.hdf5",".","reff(pc)")
    sr.rh0=myhdf5_funcs.read_attr_hdf5(fdir+"./output/ini/hdf5/dms_0.hdf5",".","rh(pc)")
    sr.mbh0=myhdf5_funcs.read_attr_hdf5(fdir+"./output/ini/hdf5/dms_0.hdf5",".","MBH")
    # dt=np.zeros((nl))
    # t=np.zeros((nl))
    objs=["star", "sbh", "ns", "wd", "bd", "rg"]
    for i in range(nl):
        fl=fdir+"./output/ecev/dms/dms_"+str(i+1)+".hdf5"
        sr.dt[i]=myhdf5_funcs.read_attr_hdf5(fl,".","dT(Myr)")/1e3
        sr.t[i]=myhdf5_funcs.read_attr_hdf5(fl,".","Time(Myr)")/1e3
        sr.tnr[i]=myhdf5_funcs.read_attr_hdf5(fl,".", "Trlx_rh(Myr)")/1e3
        sr.logt[i]=np.log10(sr.t[i])
        sr.reff[i]=myhdf5_funcs.read_attr_hdf5(fl, ".", "reff(pc)")
        sr.mbh[i]=myhdf5_funcs.read_attr_hdf5(fl,".","MBH")
        sr.logmbh[i]=np.log10(sr.mbh[i])
        sr.gas_res[i]=myhdf5_funcs.read_attr_hdf5(fl,"dMbh","gas_reservior_left",print_err=print_err)
        sr.gas_res[sr.gas_res==0]=1e-5
        sr.mcl[i]=myhdf5_funcs.read_attr_hdf5(fl,".","Mcluster(Msun)")
        sr.rh[i]=myhdf5_funcs.read_attr_hdf5(fl,".","rh(pc)")
        sr.mass_loss[i]=myhdf5_funcs.read_attr_hdf5(fl,"dMbh","mass_loss_stellar_evolution_acum",print_err=print_err)
        sr.direct_swallow[i]=myhdf5_funcs.read_attr_hdf5(fl,"dMbh","mass_direct_swallow_acum",print_err=print_err)
        sr.td_accum[i]=np.sum(myhdf5_funcs.read_attr_hdf5(fl, "dMbh","td_disc_acum",print_err=print_err))
        sr.emri_accum[i]=np.sum(myhdf5_funcs.read_attr_hdf5(fl, "dMbh","emri_acum",print_err=print_err))
        sr.emax_accum[i]=np.sum(myhdf5_funcs.read_attr_hdf5(fl, "dMbh","emax_direct_acum",print_err=print_err))
        
        for key in sr.obj:
            sr.obj[key].readhdf5(fl,key,"ls",i,print_err=print_err)

            sr.obj[key].d["ls"].t[i]=sr.t[i]
            sr.obj[key].d["ls"].dt[i]=sr.dt[i]
            sr.obj[key].d["ls"].logt[i]=sr.logt[i]

            sr.obj[key].readhdf5(fl,key,"td",i,print_err=print_err)
            sr.obj[key].d["td"].t[i]=sr.t[i]
            sr.obj[key].d["td"].dt[i]=sr.dt[i]
            sr.obj[key].d["td"].logt[i]=sr.logt[i]

            sr.obj[key].readhdf5(fl,key,"emris",i,print_err=print_err)
            sr.obj[key].d["emris"].t[i]=sr.t[i]
            sr.obj[key].d["emris"].dt[i]=sr.dt[i]
            sr.obj[key].d["emris"].logt[i]=sr.logt[i]

            sr.obj[key].readhdf5(fl,key,"emax",i,print_err=print_err)
            sr.obj[key].d["emax"].t[i]=sr.t[i]
            sr.obj[key].d["emax"].dt[i]=sr.dt[i]
            sr.obj[key].d["emax"].logt[i]=sr.logt[i]
            
            # sr.obj[key].fmden


        td=myhdf5_funcs.read_attr_hdf5(fl,"dMbh","td_disc_acum",print_err=print_err)
        if(td is None):
            td=np.zeros(6)
        for s,id in zip(objs,range(6)):
            sr.obj[s].d["td"].dmbh[i]=td[id]

        
        lc=myhdf5_funcs.read_attr_hdf5(fl,"dMbh","lc_direct_acum",print_err=print_err)
        if(lc is None):
            lc=np.zeros(6)
        for s,id in zip(objs,range(6)):
            sr.obj[s].d["ls"].dmbh[i]=lc[id]

        
        emri=myhdf5_funcs.read_attr_hdf5(fl,"dMbh","emri_acum",print_err=print_err)
        if(emri is None):
            emri=np.zeros(2)
        sr.obj["sbh"].d["emris"].dmbh[i]=emri[1]

        emax=myhdf5_funcs.read_attr_hdf5(fl,"dMbh","emax_direct_acum",print_err=print_err)
        if(emax is None):
            emax=np.zeros(2)
        sr.obj["sbh"].d["emax"].dmbh[i]=emax[1]

        
    # sr.logt=np.log10(sr.t)    
    fl=fdir+"./output/ecev/dms/dms_"+str(nl)+".hdf5"
    for key in objs:
        ri=np.log10(sr.reff[-1]/sr.r0)
        sr.obj[key].alpha_mean=get_alpha_mean(fl,ri,key,print_err=print_err)
    
    return sr

def get_ni_ne(fdir,ti,tf):
    da=np.loadtxt(fdir+"/output/run_time",skiprows=1)
    ix=da[:,0]
    t=da[:,1]
    mask=np.logical_and(t>=ti,t<=tf)
    isel=ix[mask]
    # print(isel)
    return int(isel[0]),int(isel[-1])
def get_nf(fdir):
    nf=int(np.loadtxt(fdir+"/output/run_summary.txt"))
    return nf

def get_alpha_mean(fdir,ri,objtype,print_err=True):
    x,y=myhdf5_funcs.read_2d_hdf5(fdir, objtype+"/fslope",print_err=print_err)
    alpha_star=0
    if(x is not None):
        f=interpolate.interp1d(x,y,kind="linear")
        
        for r in np.linspace(max(-3.+ri,x[0]), min(-1.+ri,x[-1]), 100):
            #print(r, f(r))
            alpha_star=alpha_star+f(r)
        alpha_star=alpha_star/100e0
    return alpha_star

def getdstr_range(fdir,fe,ni,ne,dataname,col_n=96,col_sel=' FXW'):
    nl=ne-ni+1
    t=np.zeros((nl))
    fdstr=np.zeros((col_n,nl))
    xdstr=np.zeros((col_n,nl))
    xc=np.zeros((col_n,nl))
    for i in range(0,nl):
        #print(i)
        xds,fd=myhdf5_funcs.read_2d_hdf5(fdir+str(i+ni)+fe+".hdf5",dataname,col_sel=col_sel)
        if(xds is not None):
            # print(j, fd)
            fdstr[:,i]=fd
            xdstr[:,i]=xds
            xc[:,0]=xds
            # t[j]=myhdf5_funcs.read_attr_hdf5(fdir+"output/ecev/dms/dms_"+str(i+ni)+".hdf5",".","Time(Myr)")
    # print("t0=",t[0])
    return xc, fdstr


def get2ddstr_range(fdir,fe,ni,ne,dataname,subdir="FXYW", col_n=30, col_y=None):
    if(col_y is None):
        col_y=col_n
    nl=ne-ni+1
    t=np.zeros((nl))
    fdstr=np.zeros((col_y,col_n, nl))
    xdstr=np.zeros((col_n, nl))
    ydstr=np.zeros((col_y, nl))
    j=0
    for i in range(0,nl):
        print("get2ddstr_", i+ni)
        xds,yds,fd=myhdf5_funcs.read_contour_hdf5(fdir+str(i+ni)+fe+".hdf5",dataname,subdir)
        if(xds is not None):
            # print(j, fd)
            fdstr[:,:,j]=fd
            xdstr[:,j]=xds
            ydstr[:,j]=yds
            # t[j]=myhdf5_funcs.read_attr_hdf5(fdir+"output/ecev/dms/dms_"+str(i+ni)+".hdf5",".","Time(Myr)")
            j+=1
    # print("t0=",t[0])
    return xdstr, ydstr, fdstr
def find_nearest(array, value):
	array = np.asarray(array)
	idx = (np.abs(array - value)).argmin()
	return array[idx]

def get_rc_t(fdir):
	nl=int(np.loadtxt(fdir+"/output/run_summary.txt"))
	
	t=np.zeros(nl)
	rc=np.zeros(nl)
	rc2=np.zeros(nl)
	
	r0=myhdf5_funcs.read_attr_hdf5(fdir+"/output/ini/hdf5/dms_0.hdf5",".","r0(pc)")
	for i in range(nl):
		fl=fdir+"/output/ecev/dms/dms_"+str(i+1)+".hdf5"
		mbh=myhdf5_funcs.read_attr_hdf5(fl,".","Mbh")
		nh0=mbh/3.1**3
		
		trlx0=myhdf5_funcs.read_attr_hdf5(fl,".","Trlx_rh(Myr)")
		t[i]=myhdf5_funcs.read_attr_hdf5(fl,".","Time(Myr)")/trlx0
		r,dn=myhdf5_funcs.read_2d_hdf5(fl,"all/fmden")
		f=interpolate.interp1d(r,dn)    
		
		r1,dn1=myhdf5_funcs.read_2d_hdf5(fl,"1/star/fmden")
		r2,dn2=myhdf5_funcs.read_2d_hdf5(fl,"2/sbh/fmden")
		y=myfunc.get_a_trans(r1,dn1,dn2)
		if(i>1):
			if(np.min(y)>r1[-2]):
				#print(">",y[0],r1[-1])
				rc[i]=r1[0]-10
			else:
				rc[i]=find_nearest(y,rc[i-1])
		else:
			if(np.min(y)>r1[-2]):
				#print(">",y[0],r1[-1])
				rc[i]=r1[0]-10
			else:
				rc[i]=y[0]
		#print(rc[i],y,r1[-2])
		y=myfunc.get_a_trans(r1,dn1,dn2*10)
		if(i>1):
			if(np.min(y)>r1[-2]):
				#print(">",y[0],r1[-1])
				rc2[i]=r1[0]-10
			else:
				rc2[i]=find_nearest(y,rc2[i-1])
		else:
			if(np.min(y)>r1[-2]):
				#print(">",y[0],r1[-1])
				rc2[i]=r1[0]-10
			else:
				rc2[i]=y[0]
		#print(rc2[i],y,r1[-2])				
	rc[np.logical_or(rc<r1[1],rc>r1[-3])]=None
	rc2[np.logical_or(rc<r1[1],rc>r1[-3])]=None
	
	return t, 10**rc*r0,10**rc2*r0

def get_rc_t_one(fl,m1,m2):

	r0=myhdf5_funcs.read_attr_hdf5(fl,".","r0(pc)")
    
	mbh=myhdf5_funcs.read_attr_hdf5(fl,".","Mbh")
	
	r,dn=myhdf5_funcs.read_2d_hdf5(fl,"all/fmden")
	f=interpolate.interp1d(r,dn)    
	
	r1,dn1=myhdf5_funcs.read_2d_hdf5(fl,"1/star/fmden")
	r2,dn2=myhdf5_funcs.read_2d_hdf5(fl,"2/sbh/fmden")
	y=myfunc.get_a_trans(r1,dn1,dn2)

	rc=y[0]
	
	y=myfunc.get_a_trans(r1,dn1*m1,dn2*m2)
	
	rc2=y[0]
	
	return 10**rc*r0,10**rc2*r0

def get_Rtrans_from_data(r,fden1,m1,fden2,m2):
    return myfunc.get_a_trans(r,fden1*m1,fden2*m2)

def get_influence_radius(fma:S1d,mbh,radius_factor=2):
     
     return fma.get_x(mbh*radius_factor)