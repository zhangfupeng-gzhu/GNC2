from matplotlib import pyplot as plt
from matplotlib.ticker import MaxNLocator, AutoMinorLocator,MultipleLocator, LogLocator
import matplotlib.ticker
import numpy as np
from matplotlib.ticker import ScalarFormatter

def set_yaxis_offset(OFFSET_ON=False):
	ax.get_yaxis().get_major_formatter().set_useOffset(OFFSET_ON)
	
def set_sci_format(ax,axis,scilimit=None):
	#axis='x', 'y' or 'both'
	#scilimit=(small number limit, large number limit)
	if(scilimit == None):
		scilimit=(-2,3)
	ax.ticklabel_format(style="sci",scilimits=scilimit,axis=axis)
	if(axis=='x' or 'both'):
		ax.xaxis.major.formatter._useMathText = True
	if(axis=='y' or 'both'):
		ax.yaxis.major.formatter._useMathText = True
	
def set_xaxis_offset(OFFSET_ON=False):
	ax.get_xaxis().get_major_formatter().set_useOffset(OFFSET_ON)

def set_xyaxis(ax=None, xmajorstep=None,xminorstep=None, \
			   ymajorstep=None,yminorstep=None, \
			   xmajorticks=5,ymajorticks=5,xminorticks=5,\
			   yminorticks=5,scifmt=False, xvisible=True, yvisible=True):
	plt.minorticks_on()
	if(ax==None):
		ax=plt.gca()
	if(scifmt==True):
		set_sci_format(ax,'both')
	if(xmajorstep == None):
		ax.xaxis.set_major_locator(MaxNLocator(nbins=xmajorticks))
	else:
		ax.xaxis.set_major_locator(MultipleLocator(base=xmajorstep))
	if(xminorstep == None):
		ax.xaxis.set_minor_locator(AutoMinorLocator(n=xminorticks))
	else:
		ax.xaxis.set_minor_locator(MultipleLocator(base=xminorstep))
		
	if(ymajorstep == None):	
		ax.yaxis.set_major_locator(MaxNLocator(nbins=ymajorticks))
	else:
		ax.yaxis.set_major_locator(MultipleLocator(base=ymajorstep))		
	if(yminorstep == None):
		ax.yaxis.set_minor_locator(AutoMinorLocator(n=yminorticks))
	else:	
		ax.yaxis.set_minor_locator(MultipleLocator(base=yminorstep))
	if(yvisible==False):
		ax.yaxis.set_ticklabels([])
	if(xvisible==False):
		ax.xaxis.set_ticklabels([])

def set_yaxis_log(ax=None):
	if(ax==None):
		ax=plt.gca()
	ax.yaxis.set_major_locator(LogLocator(base=10,numticks=100))
def set_xaxis_log(ax=None):
	if(ax==None):
		ax=plt.gca()
	ax.xaxis.set_major_locator(LogLocator(base=10,numticks=100))
	subs=np.linspace(0.1,0.9,9)
	ax.xaxis.set_minor_locator(LogLocator(base=10,subs=subs,numticks=100))
	
def set_yaxis(ax=None, ymajorticks=5,yminorticks=5, ymajorstep=None, yminorstep=None,scifmt=False,visible=True):
	#plt.minorticks_on()
	if(ax==None):
		ax=plt.gca()
	if(scifmt==True):
		set_sci_format(ax,'y')
	if(ymajorstep == None):	
		ax.yaxis.set_major_locator(MaxNLocator(nbins=ymajorticks))
	else:
		ax.yaxis.set_major_locator(MultipleLocator(base=ymajorstep))
	if(yminorstep == None):
		ax.yaxis.set_minor_locator(AutoMinorLocator(n=yminorticks))
	else:	
		ax.yaxis.set_minor_locator(MultipleLocator(base=yminorstep))
	if(visible==False):
		ax.yaxis.set_ticklabels([])
	
def set_xaxis(ax=None,xmajorticks=5,xminorticks=5, xmajorstep=None, xminorstep=None,scifmt=False,visible=True):
	#plt.minorticks_on()
	if(ax==None):
		ax=plt.gca()
	if(scifmt==True):
		set_sci_format(ax,'x')
	if(xmajorstep == None):
		ax.xaxis.set_major_locator(MaxNLocator(nbins=xmajorticks))
	else:
		ax.xaxis.set_major_locator(MultipleLocator(base=xmajorstep))
	if(xminorstep == None):	
		ax.xaxis.set_minor_locator(AutoMinorLocator(n=xminorticks))
	else:
		ax.xaxis.set_minor_locator(MultipleLocator(base=xminorstep))	
	if(visible==False):
		ax.xaxis.set_ticklabels([])
	
def set_subgrid_xy_label(numx, numy, i, xlim, ylim, xyoff):
	xmin, xmax=xlim
	ymin, ymax=ylim
	xoff, yoff=xyoff
	ax=plt.gca()
	nx=i//numx
	x_enlabel=False
	y_enlabel=False
	ny=i%numx 
	by=[s for s in range(1, numx*numy+1, numy)]
	bx=[s+numy*(numx-1) for s in range(1, numy+1)]
	bxhide=[s for s in range(1, numy*(numx-1)+1)]
	byhide=[(t+(s-1)*numy) for s in range(1,numx+1) for t in range(2, numy+1)]
	
#	print bxhide,byhide
#	exit()
	
	if any(i+1==s for s in by):
		if(i+1==0):
			ax.set_ylim([ymin,ymax])
		else:
			ax.set_ylim([ymin+yoff, ymax])
		y_enlabel=True
	else:
		ax.set_ylim([ymin, ymax])

	if any(i+1==s for s in bx):
		if(i+1==numx*numy):
			ax.set_ylim([ymin,ymax])
		else:
			ax.set_xlim([xmin, xmax+xoff])
		x_enlabel=True
	else:
		ax.set_xlim([xmin, xmax])

	if any(i+1==s for s in bxhide):
		ax.xaxis.set_ticklabels([])
	
	if any(i+1==s for s in byhide):
		ax.yaxis.set_ticklabels([])
	
	return x_enlabel, y_enlabel
def enforce_log_tick(ax):
	x_major = matplotlib.ticker.LogLocator(base = 10.0, numticks = 5)
	ax.xaxis.set_major_locator(x_major)
	x_minor = matplotlib.ticker.LogLocator(base = 10.0, subs = np.arange(1.0, 10.0) * 0.1, numticks = 10)
	ax.xaxis.set_minor_locator(x_minor)
	ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
	
def texts(ax,label,ps, fontsize=None, color='k',alignment='center',**keywords):
	ax.text(ps[0],ps[1],label, \
		horizontalalignment=alignment,transform = ax.transAxes, fontsize=fontsize, color=color,**keywords)	
				
def use_tex():
     plt.rcParams.update({
    "text.usetex": True})

def set_colorbar_axisformat(cbar,use_sci=False,**arg):
	if(use_sci):
		set_sci_format(cbar.ax,'y')
	set_yaxis(cbar.ax,**arg)
	