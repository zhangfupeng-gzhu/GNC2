import scipy.ndimage as ndimage
import numpy as np
from matplotlib.collections import LineCollection
from matplotlib import pyplot as plt
import matplotlib as mpl
from scipy.interpolate import interp1d
from shapely.geometry import LineString
from scipy.interpolate import CubicSpline


def get_a_trans(x,f,g,x2=None,plot=False):
	if(x2 is None):
		x2=x
    # get the cross position of two line x-f and x-g 
	first_line = LineString(np.column_stack((x, f)))
	second_line = LineString(np.column_stack((x2, g)))
	intersection = first_line.intersection(second_line)

	if intersection.geom_type == 'MultiPoint':
		if(plot): plt.plot(*LineString(intersection).xy, 'o')
		#for i in intersection.geoms:
		#	print(i)
		x,y=LineString(intersection.geoms).xy
		#return np.asarray(x),np.asarray(y)
	elif intersection.geom_type == 'Point':
		if(plot): plt.plot(*intersection.xy, 'o')
		x,y=intersection.xy
	return np.asarray(x)

def powidx_x(x,y):
	n=len(x)
	fstep=x[1]-x[0]
	xx=np.zeros([10])
	yy=np.zeros([10])
	idx=np.zeros([n])
	#print("x=", x)
	#print("y=", y)
	f=interp1d(x,y,kind='cubic')
	for i in range(n):
		f0=x[i]
		#print(f0)
		for j in range(10):
			if(i==0):
				xx[j]=fstep*float(j)/(10.-1)+f0
			elif(i==n-1):
				xx[j]=fstep*float(j)/(10.-1)-fstep+f0
			else:
				xx[j]=fstep*float(j)/(10.-1)-fstep*0.5+f0
			#print(i, j,xx[j])
			yy[j]=np.log10(f(xx[j]))
		#print(xx)
		#print(yy)
		a=np.polyfit(xx,yy,1)
		idx[i]=a[0]
	return idx

def linear_log(x, xmax, xmin):
	return np.log10(x/xmin)/np.log10(xmax/xmin)

def smooth_cont(X, Y, Z, zoom_in):
	XX=ndimage.zoom(X, zoom=zoom_in)
	YY=ndimage.zoom(Y, zoom=zoom_in)
	ZZ=ndimage.zoom(Z, zoom=zoom_in)
	return XX, YY, ZZ

def refine_errorbar(yavg, ysct,cr=1e-10):
	yr=np.copy(ysct)
	for i in range(len(yavg)):
		if(yavg[i]-ysct[i]<0): 
			yr[i]=yavg[i]-cr
	return yr

def get_log_errorbar(logyavg,logysct):
	yavg=np.copy(logyavg);
	ysctup=np.copy(logyavg);
	ysctbt=np.copy(logyavg);
	for i in range(len(logyavg)):
		yavg[i]=10**logyavg[i];
		ysctup[i]=10**(logyavg[i]+logysct[i])-yavg[i]
		ysctbt[i]=yavg[i]-10**(logyavg[i]-logysct[i])

	return yavg, ysctup, ysctbt
	

def _process_colors_by_index(self):
    """
    Color argument processing for contouring.

    The color is based in the index in the level set, not
    the actual value of the level.

    """
    self.monochrome = self.cmap.monochrome
    if self.colors is not None:
        # Generate integers for direct indexing.
        i0, i1 = 0, len(self.levels)
        if self.filled:
            i1 -= 1
        # Out of range indices for over and under:
        if self.extend in ('both', 'min'):
            i0 = -1
        if self.extend in ('both', 'max'):
            i1 += 1
        self.cvalues = list(range(i0, i1))
        self.set_norm(colors.NoNorm())
    else:
        self.cvalues = range(len(self.levels))
    self.set_array(range(len(self.levels)))
    self.autoscale_None()
    if self.extend in ('both', 'max', 'min'):
        self.norm.clip = False

    # self.tcolors are set by the "changed" method

##examples:
#def plot_contour_level_by_index():
#	orig = matplotlib.contour.ContourSet._process_colors
#	matplotlib.contour.ContourSet._process_colors = _process_colors_by_index
#	cmap = plt.cm.jet
#	figure()
#	out = plt.contourf([[.12, .2], [.8, 2]], levels=[0, .1, .3, .5, 1, 3], cmap=cmap)
#	plt.colorbar()
	# fix what we have done
#	matplotlib.contour.ContourSet._process_colors = orig


def countlines(filename):
    f = open(filename)
    try:
        lines = 1
        buf_size = 1024 * 1024
        read_f = f.read # loop optimization
        buf = read_f(buf_size)

        # Empty file
        if not buf:
            return 0

        while buf:
            lines += buf.count('\n')
            buf = read_f(buf_size)

        return lines
    finally:
        f.close()

def plot_vary_color_line(ax, myXdata, myYdata, colordata,colormap='gist_rainbow_r', **arg):
    x   = myXdata 
    y   = myYdata
    t = colordata

    # set up a list of (x,y) points
    points = np.array([x,y]).transpose().reshape(-1,1,2)
#    print points.shape  # Out: (len(x),1,2)

    # set up a list of segments
    segs = np.concatenate([points[:-1],points[1:]],axis=1)
#    print segs.shape  # Out: ( len(x)-1, 2, 2 )
                      # see what we've done here -- we've mapped our (x,y)
                      # points to an array of segment start/end coordinates.
                      # segs[i,0,:] == segs[i-1,1,:]

    # make the collection of segments
    lc = LineCollection(segs, cmap=plt.get_cmap(colormap),**arg)
    norm = mpl.colors.Normalize(vmin=0, vmax=1)
    lc.set_norm(norm)
    lc.set_array(t) # color the segments by our parameter

    # plot the collection
    ax.add_collection(lc) # add the collection to the plot
    return lc
#    plt.xlim(x.min(), x.max()) # line collections don't auto-scale the plot
#   plt.ylim(y.min(), y.max())

def compute_curve_slope(x, y, x2):
    """
    计算曲线C在x2点处的斜率，使用三次样条插值。
    
    参数:
    x  : 原始曲线的x坐标数组（允许未排序，但需无重复）
    y  : 原始曲线的y坐标数组
    x2 : 需要计算斜率的x坐标数组
    
    返回:
    slopes : x2对应点的斜率数组，超出原始范围的位置为np.nan
    """
    # 检查输入有效性
    if len(x) != len(y):
        raise ValueError("x和y长度必须相同")
    if len(x) < 2:
        raise ValueError("至少需要2个点进行插值")
    
    # 排序并去重（保留最后一个重复值）
    sorted_idx = np.argsort(x)
    x_sorted = x[sorted_idx]
    y_sorted = y[sorted_idx]
    
    # 去重处理
    mask = np.concatenate(([True], np.diff(x_sorted) != 0))
    x_clean = x_sorted[mask]
    y_clean = y_sorted[mask]
    
    if len(x_clean) < 2:
        raise ValueError("去重后有效点数不足，无法插值")
    
    # 创建三次样条插值对象
    cs = CubicSpline(x_clean, y_clean)
    
    # 计算有效范围掩码
    valid_mask = (x2 >= x_clean.min()) & (x2 <= x_clean.max())
    
    # 初始化结果数组为nan
    slopes = np.full_like(x2, np.nan, dtype=np.float64)
    
    # 计算有效点处的导数
    if np.any(valid_mask):
        slopes[valid_mask] = cs(x2[valid_mask], nu=1)  # nu=1表示一阶导数
    
    return slopes



import math

def fBPowerlawN(alpha,xb,  c, x):
    n=len(alpha)
    for i in range(n):
        if(x<=xb[i+1] and x>xb[i]):
            return c[i]*(x/xb[i])**alpha[i]


def fCBPowerLawN_prepare(alpha, xb, c1):
    """
    计算累积分段幂律函数的参数。

    参数:
        alpha (list[float]): 每个区间的幂律指数，长度为n。
        xb (list[float]): 区间边界点，长度为n+1，必须严格递增。
        c1 (float): 第一个区间的归一化常数。

    返回:
        tuple: (c, t, Q)
            c (list[float]): 各区间归一化系数，长度为n。
            t (list[float]): 累积分布函数值，长度为n+1。
            Q (float): 总积分值，用于归一化。
    
    异常:
        ValueError: 如果xb不是严格递增的。
    """
    n = len(alpha)
    if len(xb) != n + 1:
        raise ValueError("Length of xb must be n+1 where n is the length of alpha")
    
    # 初始化c数组并检查区间边界
    c = [0.0] * n
    c[0] = c1
    for i in range(n - 1):
        if xb[i + 1] <= xb[i]:
            raise ValueError("xb must be in strictly ascending order")
        c[i + 1] = c[i] * (xb[i + 1] / xb[i]) ** alpha[i]
    
    # 计算积分段和累积分布
    s = [0.0] * (n + 1)
    for i in range(n):
        x_curr = xb[i]
        x_next = xb[i + 1]
        if alpha[i] != -1.0:
            term = (x_next / x_curr) ** (alpha[i] + 1) - 1
            a_i = c[i] * x_curr / (alpha[i] + 1) * term
        else:
            a_i = c[i] * x_curr * math.log(x_next / x_curr)
        s[i + 1] = s[i] + a_i
    
    Q = s[-1]
    t = [val / Q for val in s]
    return c, t, Q