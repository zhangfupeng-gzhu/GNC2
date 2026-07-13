from re import L
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.collections as mcoll
import matplotlib.path as mpath
import myaxisfmt
    
    
def contour_fxy(ax,x,y,z,N=100,vmax=None,vmin=None,log=False,\
    tr=False,ry=False,line=False,set_abs=False,cmap='gist_rainbow',grids=False,\
        showcb=True,**arg):
    a=np.asarray(z)
    if(set_abs):
        a=np.abs(a)
    if(log):
        a[a <= 0] =a[a>0].min()/2.
        a=np.asarray(np.log10(a))

    if(tr):
        a=np.transpose(a)
    if(ry):
        a=np.flipud(a)

    if(vmax is None):vmax=a.max()
    if(vmin is None):vmin=a.min()
    level=np.linspace(vmin,vmax,N)
    X,Y=np.meshgrid(x,y)
    #print(np.ndim(x))
    #print(np.ndim(y))
    #print(np.ndim(X))
    #print(np.ndim(Y))
    if(not line):
        cs=plt.contourf(X,Y,a,level, cmap=cmap,**arg)
    else:
        cs=plt.contour(X,Y,a,level,**arg)
    if(grids):
        plt.plot(X,Y, ".",markersize=0.5,**arg)

    if(showcb):
        cb=plt.colorbar(cs)
        myaxisfmt.set_xyaxis()   
        return cs, cb  
    else:
        myaxisfmt.set_xyaxis()   
        return cs

def get_color_cycle():
	prop_cycle = plt.rcParams['axes.prop_cycle']
	colors = prop_cycle.by_key()['color']
	return colors