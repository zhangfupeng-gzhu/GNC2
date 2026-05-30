#! /usr/bin/env python
import h5py
from matplotlib import pyplot as plt
import numpy as np
import math
from matplotlib import ticker,cm,colors
import myaxisfmt
from scipy import interpolate, misc,optimize
import h5py
import gwcom
import pandas as pd

myaxisfmt.use_tex()

def kernel_reg(t, rate,window_size):
    ysmooth=pd.Series(rate).rolling(window=window_size,center=True).mean()
    return ysmooth.values

def get_smooth_nr(sr,objtype,evtype,rawdatasmooth=1):
    tf=sr.obj[objtype].d[evtype].t
    rate=sr.obj[objtype].d[evtype].rate
    rf0=kernel_reg(np.log10(tf),np.log10(rate),rawdatasmooth)
    return np.log10(tf), rf0

plt.figure(figsize=(5,4))
plt.clf()

cana1='#FF7F0E'

sx=1;sy=1

ax=plt.subplot(sx,sy,1)
nl=int(np.loadtxt("../output/run_summary.txt"))

fdir="../"
nl=int(np.loadtxt(fdir+"/output/run_summary.txt"))
sr=gwcom.get_one_rates(fdir, nl,print_err=False)
nmv=2
logt,logr=get_smooth_nr(sr,"star","td", nmv)
# print(logr)
ax.plot(logt,logr,"-",label="TDE(star)") 

logt,logr=get_smooth_nr(sr,"sbh","ls", nmv)
# print(logr)
ax.plot(logt,logr,"-",label="Direct Swallow (SBH)") 

logt,logr=get_smooth_nr(sr,"sbh","emris", nmv)
ax.plot(logt,logr,":",label="EMRIs(SBH)")

ax.legend(loc=[0.1,0.6] )
#ax.set_xlim([0,14])
#ax.set_ylim(1e-5,1e-3)
# ax.set_yscale("log")
ax.set_ylabel("Rates (yr$^{-1}$)",fontsize=12)
ax.set_xlabel("log $t$ (Gyr)",fontsize=12)
myaxisfmt.set_xyaxis()

plt.tight_layout()
plt.savefig("fig_rates.pdf")