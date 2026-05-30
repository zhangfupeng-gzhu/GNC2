#! /usr/bin/env python
from matplotlib import pyplot as plt
import numpy as np
import myaxisfmt 
import gwcom
import pandas as pd
import myhdf5_funcs

myaxisfmt.use_tex()
 

plt.figure(figsize=(12,4))
plt.clf()

cana1='#FF7F0E'

sx=1;sy=3


ax=plt.subplot(sx,sy,1)

fdir="../"
nl=int(np.loadtxt(fdir+"/output/run_summary.txt"))
x,fmden= myhdf5_funcs.read_2d_hdf5(fdir+"/output/ini/hdf5/dms_0.hdf5","1/star/fmden")
ax.plot(x,fmden*10**(x*2),".", label="GNC (Init)")
x,y=myhdf5_funcs.read_2d_hdf5(fdir+"/output/ini/hdf5/ini.hdf5", "frho/1_ini_frho")
ax.plot(x,y*10**(x*2),label="Dehnen")

x,fmden= myhdf5_funcs.read_2d_hdf5(fdir+"/output/ecev/dms/dms_"+str(nl)+".hdf5","1/star/fmden")
ax.plot(x,fmden*10**(x*2),"-k", label="GNC (12Gyr)")

sl=10
x,fmden= myhdf5_funcs.read_2d_hdf5(fdir+"/output/ecev/dms/dms_"+str(sl)+".hdf5","1/star/fmden")
t=myhdf5_funcs.read_attr_hdf5(fdir+"/output/ecev/dms/dms_"+str(sl)+".hdf5", ".", "Time(Myr)")
lb=format(t[0]/1000,".1f")
ax.plot(x,fmden*10**(x*2),"--", label="GNC ("+lb +"Gyr)")

ax.legend(  )

ax.set_yscale("log")
ax.set_ylabel("$\\rho r^2$",fontsize=12)
ax.set_xlabel("log $r$",fontsize=12)
myaxisfmt.set_xaxis()

ax=plt.subplot(sx,sy,2) 

sr=gwcom.get_one_rates(fdir, nl,print_err=False)

ax.plot(sr.t,sr.obj["star"].d["td"].rate,".",label="TDE(star)")  

ax.legend(  )
ax.set_yscale("log")
ax.set_ylim([1e-6,1e-3])
ax.set_ylabel("Rates (yr$^{-1}$)",fontsize=12)
ax.set_xlabel("$t$ (Gyr)",fontsize=12)
myaxisfmt.set_xaxis()

ax=plt.subplot(sx,sy,3) 

ax.plot(sr.t,sr.reff,".",label="Effective radius")  
ax.plot(sr.t,sr.rh,".",label="Influence Radius")  

ax.legend(  ) 
ax.set_ylabel("pc",fontsize=12)
ax.set_xlabel("$t$ (Gyr)",fontsize=12)
myaxisfmt.set_xyaxis()

plt.tight_layout()
plt.savefig("fig_evl.pdf")