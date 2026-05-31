#--global
LOC=zfp
COMPILE=gfortran
#COMPILE=ifort

#hdf5dirlib=/usr/lib/x86_64-linux-gnu/hdf5/serial/
hdf5dirlib=$(HDF5LIBDIR)
#hdf5dirinc=/usr/include/hdf5/serial/
hdf5dirinc=$(HDF5INCDIR)

use_new=1

ifort_gfortran_comparability=TRUE
use_hdf5=1
performance_analyze=0

OUT=libmycom2.a
libdir=$(FSHAREDIR)
 
FC=gfortran
dir:=$(libdir)glib/com/
#FLAG= -fopenmp -O3 -ffast-math -fno-cx-limited-range -fno-math-errno -funroll-loops --param max-unroll-times=2
FLAG= -fopenmp -O3 -ffast-math 
FLAG+=-fbacktrace -fcheck=all -Waliasing -Wrealloc-lhs-all -cpp -malign-double
ifeq ($(use_new), 1)
FLAG+=-fallow-invalid-boz
endif
ifeq ($(use_hdf5), 1)
FLAG+= -DUSE_HDF5
	FLAG+= -I$(hdf5dirinc) -L$(hdf5dirlib) -lhdf5_fortran -lhdf5hl_fortran -lhdf5 -lhdf5_hl
endif
ifeq ($(performance_analyze), 1)
	FLAG+= -pg -g 
endif
ifeq ($(performance_analyze), 2)
	FLAG+= -fprofile-arcs -ftest-coverage
endif

#----OBJLIST
VPATH= ingl: ingl/Gauss-Kronrod :  local_minimum : sts : math :app :file :  differential : random

OBJ_APP=APP_CONST.o  app_star.o rtnewt.o rtsafe.o app_celestial.o \
    app_physics.o 

OBJ_INGL=ingl_dopri5.o d1mach.o \
    dqawse.o dqaws.o \
    dqc25s.o dqcheb.o \
    dqk15w.o dqmomo.o \
    dqpsrt.o dqwgts.o \
    fdump.o	 i1mach.o \
    j4save.o xermsg.o \
    xerprn.o xersve.o \
    xerhlt.o xercnt.o dqagp.o dqagpe.o \
	dqags.o dqagse.o dqk21.o dqelg.o xgetua.o intgl.o  ingl_dopri853.o  
	

OBJ_FILE= FILE_STRING.o 
ifeq ($(use_hdf5),1)
	OBJ_FILE+=file_md_hdf5_table.o file_md_hdf5_file.o  file_md_hdf5_group.o file_md_hdf5.o \
	file_md_hdf5_attr.o
endif
OBJ_STS=sts_bin2.o sts_d2-1.o sts_fc.o   sts_sort.o \
		indexx.o  spear.o  betai.o crank.o gammln.o betacf.o \
		 gammp.o sts_d2-2.o md_fc.o\
		sts_type_s1d_com.o sts_type_s2d_com.o sts_type_com.o \
                 md_sts_d21.o sts_basic_type.o md_s1d_hst.o sts_type_s1d_ird_com.o \
				 sts_type_s2d_ird_com.o md_s1d_ird_hst.o sts_fc_ird.o sts_type_s1d_adv_com.o \
				 md_s2d_ird_hst.o

OBJ_SPL= spline_sub_mylib.o splin2.o splie2.o  linear_interpolate.o  

OBJ_FIT= gcf.o gser.o  

OBJ_FUN=FUN_RANDOM.o   fun_math.o sub_math.o   gaussian.o  sub_arr.o 
        
OBJ_STR=str.o

OBJ_DYDX=differential.o

OBJ_LMIN= rtbis.o   root_finding.o

OBJ_MATH=rs.o tql2.o tqlrat.o tred1.o tred2.o pythag.o  matrix.o
 

OBJS=$(OBJ_APP) $(OBJ_INGL) $(OBJ_FILE) $(OBJ_STS) $(OBJ_FUN) \
	$(OBJ_STR) $(OBJ_LMIN) $(OBJ_SPL) $(OBJ_FIT) $(OBJ_DYDX) \
	$(OBJ_MATH)   $(OBJ_INGL)


#---make target
all: $(OUT) 
 

md_s1d_hst.o:  sts_type_s1d_com.o
sts_type_s1d_com.o: sts_basic_type.o 
md_s1d_ird_hst.o: sts_type_s1d_ird_com.o
md_s2d_ird_hst.o: sts_type_s2d_ird_com.o
md_fc.o: sts_type_s1d_com.o md_s1d_hst.o
sts_type_s2d_com.o: sts_type_s1d_com.o md_s1d_hst.o
sts_type_com.o : sts_type_s1d_com.o sts_type_s2d_com.o md_s1d_hst.o\
		 md_fc.o sts_type_s2d_ird_com.o sts_type_s1d_ird_com.o md_s1d_ird_hst.o \
		 sts_type_s1d_adv_com.o md_s2d_ird_hst.o
sts_type_s2d_ird_com.o: sts_type_s1d_ird_com.o  file_md_hdf5.o
sts_type_s1d_adv_com.o: sts_basic_type.o 

file_md_hdf5.o: file_md_hdf5_table.o file_md_hdf5_file.o  file_md_hdf5_group.o  \
	file_md_hdf5_attr.o 

differential.o: sts_type_com.o


$(OUT):$(OBJS)
	ar -rc $(OUT) $(OBJS)
	cp $(OUT) $(dir)
	cp *.mod  $(dir)
%.o:%.f90 
	$(FC) -c $(FLAG) $<

%.o:%.for 
	$(FC) -c $(FLAG) $<

%.o:%.f 
	$(FC) -c $(FLAG) $<

.PHONY:clean
clean:
	rm -f *.o
	rm -f *.mod
	rm -f *.a
	rm -f *.f90~
	rm -f *.for~


