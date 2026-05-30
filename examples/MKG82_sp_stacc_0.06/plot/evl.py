#! /usr/bin/env python
import h5py
from matplotlib import pyplot as plt
import numpy as np
import math
from matplotlib import ticker,cm,colors
import myaxisfmt 
import h5py 
import myplt_funcs 
from scipy.interpolate import splrep, BSpline
import pandas as pd
import gwcom

colors=myplt_funcs.get_color_cycle()
colorstar=colors[0]
colorsbh=colors[1]
colorwd=colors[2]
colorns=colors[3]
colorbd=colors[4]
colorrg= colors[5]
 

def smooth_rate_nmin(t,nw0):
    rf_tmp=nw0.copy()
    t_tmp=t.copy()
    nw=nw0[0]; j=0; 
    tp=0
    for i in range(len(t)):
        if(i!=len(t)-1):
            if(nw0[i+1]!=0):
                tn=t[i]
                t_tmp[j]=(tp+tn)/2.
                rf_tmp[j]=nw/(tn-tp)/1e9
                j+=1
                tp=tn
                nw=nw0[i+1]
        elif(i==len(t)-1 and nw0[i]!=0):
            t_tmp[j]=(t[i]-t[i-1])/2.
            rf_tmp[j]=nw0[i]/(t[i]-t[i-1])/1e9

    rf=rf_tmp[:j].copy()
    tf=t_tmp[:j].copy()
    return tf,rf

def kernel_reg(t, rate,window_size):
    ysmooth=pd.Series(rate).rolling(window=window_size,center=True).mean()
    return ysmooth.values

def plot_one_t_mbh(sr,rawdatasmooth=1,epsilon=0):
    logt=np.log10(sr.obj["star"].d["td"].t[::rawdatasmooth])

    dmbh_edd=sr.mbh[0]*np.exp(22*10**logt*(1-epsilon))
    
    ax.plot(logt,np.log10(dmbh_edd[::rawdatasmooth]),"--k",label="MBH mass (Edd)",linewidth=2)
    ax.plot(logt,np.log10(sr.mbh[::rawdatasmooth]),"-k",markersize=0.5,label="MBH mass",linewidth=2)
    ax.plot(logt,np.log10(sr.obj["star"].d["td"].dmbh[::rawdatasmooth]),":",color=colorstar, label="TDE (MS)")
    ax.plot(logt,np.log10(sr.obj["rg"].d["td"].dmbh[::rawdatasmooth]),":",color=colorrg, label="TDE (Post-MS)")
    emris=sr.obj["sbh"].d["emris"].dmbh[::rawdatasmooth]
    emax=sr.obj["sbh"].d["emax"].dmbh[::rawdatasmooth]
    ax.plot(logt,np.log10(emris+emax),"-",color=colorsbh, label="EMRI (SBHs)")
    ax.plot(logt,np.log10(sr.obj["star"].d["ls"].dmbh[::rawdatasmooth]),"--",color=colorstar, label="LC (MS)")
    ax.plot(logt,np.log10(sr.obj["sbh"].d["ls"].dmbh[::rawdatasmooth]),"--",color=colorsbh,label="LC (SBHs)")
    ax.plot(logt,np.log10(sr.obj["ns"].d["ls"].dmbh[::rawdatasmooth]),"--",color=colorns,label="LC (NS)")
    ax.plot(logt,np.log10(sr.obj["wd"].d["ls"].dmbh[::rawdatasmooth]),"--",color=colorwd,label="LC (WD)")
    ax.plot(logt,np.log10(sr.gas_res[::rawdatasmooth]),"-",color="c",label="Gas Reservoir (left)",linewidth=2)
    ax.plot(logt,np.log10(sr.mcl[::rawdatasmooth]),"-",color="g",label="$M_{\\rm cl}$",linewidth=2)

    print("model,MHB(final)=",fdir,sr.mbh[-1])
    myaxisfmt.set_xyaxis(ymajorstep=1,xmajorstep=1)
    ax.set_ylim([1,np.log10(max(sr.mbh[-1],sr.mcl[0]))*1.05])

def get_smooth_nr(sr,objtype,evtype,rawdatasmooth=1):
    tf,rf=smooth_rate_nmin(sr.obj[objtype].d[evtype].t,sr.obj[objtype].d[evtype].nw)
    rf0=kernel_reg(np.log10(tf),np.log10(rf),rawdatasmooth)
    return np.log10(tf), rf0


def plot_rate(sr,obtype,evtype,line='-',color='k',label=None):
    emri=sr.obj[obtype].d[evtype].rate
    logt=np.log10(sr.obj[obtype].d[evtype].t)
    ax.plot(logt,np.log10(emri),linestyle=line,color=color,label=label)

def plot_one_t_rate(sr,nmv=5):
    # srl=sr.get_smooth(nmv,method="mov_nozero")

    logt,logr=get_smooth_nr(sr,"star","td", nmv)
    ax.plot(logt,logr,":",color=colorstar,label="TDE(MS)")
    logt,logr=get_smooth_nr(sr,"bd","td", nmv)
    ax.plot(logt,logr,":",color=colorbd,label="TDE(BD)")
    logt,logr=get_smooth_nr(sr,"rg","td", nmv)
    ax.plot(logt,logr,":",color=colorrg,label="TDE(Post-MS)")

    # plot_rate(srl,"star","td",line=':',color=colorstar, label="TD(MS)")
    # plot_rate(srl,"bd","td",line=':',color=colorbd, label="TD(BD)")
    # plot_rate(srl,"rg","td",line=':',color=colorrg, label="TD(Post-MS)")

    # logt,logr=get_smooth_nr(sr,"star","emris", nmv)
    # ax.plot(logt,logr,"-",color=colorstar,label="EMRI(MS)")
    # logt=np.log10(srl.obj["star"].d["emris"].t)
    # emri=srl.obj["star"].d["emris"].rate
    # ax.plot(logt,np.log10(emri),"-",color=colorstar,label="EMRI(MS)")
    logt,logr=get_smooth_nr(sr,"star","emris",nmv)
    ax.plot(logt, logr,"-",color=colorstar,label="EMRI(MS)")

    # plot_rate(srl,"star","emris",color=colorstar, label="EMRI(MS)")

    # logt,logr=get_smooth_nr(sr,"bd","td")
    # ax.plot(logt,logr,".",color="0.6",zorder=-100, \
    #             label="TD (BD)")

    logt,logr=get_smooth_nr(sr,"bd","emris",nmv)
    ax.plot(logt, logr,"-",color=colorbd,label="EMRI(BD)")

    # plot_rate(srl,"bd","emris",color=colorbd,label= "EMRI(BD)")

    # logt,logr=get_smooth_nr(sr,"sbh","emris")
    # ax.plot(logt,logr,".",color="r",zorder=-100, \
    #             label="EMRIs (SBHs)")
    logt,logr=get_smooth_nr(sr,"sbh","emris",nmv)
    ax.plot(logt, logr,"-",color=colorsbh,label="EMRI(SBHs)")
    # logt,logr=get_smooth_nr(sr,"sbh","emax",nmv)
    # ax.plot(logt, logr,"--k",label="EMAX(SBHs)")
    # emri=sr.obj["sbh"].d["emris"].rate
    # logt=np.log10(sr.obj["sbh"].d["emris"].t)
    # ax.plot(logt,np.log10(emri),".k",label="EMRI(SBH)")

    # plot_rate(srl,"sbh","emris",color=colorsbh,label= "EMRI(SBH)")
    
    # emri=srl.obj["wd"].d["emris"].rate
    # ax.plot(srl.logt, np.log10(emri),"-",color=colorwd,label="EMRI(WD)")

    logt,logr=get_smooth_nr(sr,"wd","emris",nmv)
    ax.plot(logt, logr,"-",color=colorwd,label="EMRI(WD)")

    # plot_rate(srl,"wd","emris",color=colorwd,label= "EMRI(WD)")

    # logt,logr=get_smooth_nr(sr,"ns","emris")
    # ax.plot(logt,logr,".",color="0.1",zorder=-100, \
    #             label="EMRIs (NSs)")
    logt,logr=get_smooth_nr(sr,"ns","emris",nmv)
    ax.plot(logt, logr,"-",color=colorns,label="EMRI(NS)")
    # plot_rate(srl,"ns","emris",color=colorns, label="EMRI(NS)")


    myaxisfmt.set_xaxis(xmajorstep=1)

myaxisfmt.use_tex()

plt.figure(figsize=(8,4))
plt.clf()

# models=["MKG82_sp_stacc_0.06"]
# labels=["MKG82F0.06k"]

models="../"
nmvs=10


fdir=models+"/" 
nl=int(np.loadtxt(fdir+"/output/run_summary.txt"))
sr=gwcom.get_one_rates(fdir, nl,print_err=False)

ylabel="log $M_\\bullet$ ($M_\\odot$)"
sx=1;sy=2
j=0
labelfontsize=12 
j=j+1
ax=plt.subplot(sx,sy,j)
fdir=models+"/"
nl=int(np.loadtxt(fdir+"/output/run_summary.txt"))
#nl=200
plot_one_t_mbh(sr,epsilon=0.1) 
ax.set_xlabel("log $t$ (Gyr)",fontsize=labelfontsize)
ax.set_ylabel(ylabel,fontsize=labelfontsize)
if(j==1):
	# legend_para={"handlelength":1.5}
	ax.legend(loc='lower left',ncol=3,fontsize=4,handlelength=3.2)

j=j+1
ax=plt.subplot(sx,sy,j)
plot_one_t_rate(sr,nmv=nmvs) 
ax.set_xlabel("log $t$ (Gyr)",fontsize=labelfontsize)
ax.set_ylabel("log Rate (yr$^{-1}$)",fontsize=labelfontsize) 
myaxisfmt.set_yaxis(ymajorstep=1.0)

plt.tight_layout() 
plt.savefig("fig_mbh_evl_dehnen.pdf")
