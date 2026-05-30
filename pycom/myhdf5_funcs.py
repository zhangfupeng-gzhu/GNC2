import h5py
import numpy as np
import myplt_funcs
def read_attr_hdf5(fl,daname, attrname,print_err=True):
	f=h5py.File(fl,'r')
	if( daname in f):
		da=f[daname]
	else:
		if(print_err):
			print("daname", daname, " not existed, return")
		return None
	if attrname in da.attrs:
		return da.attrs[attrname][:]
	else:
		if(print_err):
			print("attrname", attrname, " not existed, return")
		return None



def show_data_info(fl,dir):
	with h5py.File(fl, 'r') as file:
		# 指定要查看的表格路径
		dataset_path = dir  # 替换为你的表格路径
		
		try:
			# 获取表格数据集
			dataset = file[dataset_path]
			
			# 打印表格基本信息
			print(f"表格路径: {dataset_path}")
			print(f"形状: {dataset.shape}")  # (行数, 列数)
			print(f"数据类型: {dataset.dtype}")
			
			# 检查是否为结构化数组（有列名）
			if hasattr(dataset.dtype, 'names') and dataset.dtype.names:
				print(f"列名: {list(dataset.dtype.names)}")
			
			# 打印属性信息（可能包含表头信息）
			if dataset.attrs:
				print("属性信息:")
				for key, value in dataset.attrs.items():
					print(f"  {key}: {value}")
			
			# 预览前几行数据（可选）
			print("\n前5行数据预览:")
			print(dataset[:5])  # 显示前5行
		except KeyError:
			print(f"错误：找不到表格 '{dataset_path}'")
			print("可用的表格路径:")
			# 只列出所有数据集路径
			def list_datasets(name, obj):
				if isinstance(obj, h5py.Dataset):
					print(f"  - {name}")
			file.visititems(list_datasets)

def read_table_hdf5(fl,dir):
	f=h5py.File(fl,'r')
	
	if dir in f:
		return f[dir]
	else:
		return None
	
def read_2d_hdf5(fl, dir, col_sel="  FX",print_err=True):
	f=h5py.File(fl,'r')
	if dir in f:
		xx=f[dir]["   X"]
		yy=f[dir][col_sel]
	else:
		if(print_err): 
			print(dir+" doest not exist in file"+fl)
		xx=None
		yy=None
	return xx, yy

def test_path_exist_hdf5(fl,dir):
	f=h5py.File(fl,'r')
	if dir in f:
		return True
	else:
		return False

def read_contour_hdf5(fl, sdir,subdir):
	f=h5py.File(fl, 'r' )
	if sdir in f:
		z=f[sdir][subdir]
		x=f[sdir]['X']
		y=f[sdir]['Y']
		#print("shape x=", np.shape(x))
		#print(len(y),len(z))
		if(np.ndim(x)>=2):
		#	print(np.ndim(x))
			return x[:,0],y[:,0],z[:]
		else:
		#	print(np.ndim(x))
			return x,y,z
	else:
		return None,None,None

def plot_contour_hdf5(ax, fl, sdir, subdir="FXY", **arg):
	x,y,z=read_contour_hdf5(fl, sdir,subdir)
	return x,y,z, myplt_funcs.contour_fxy(ax,x,y,z,**arg)