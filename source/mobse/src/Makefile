

#mobse:	clean\
#	$(OBJT1) $(LFLAGS)
#	$(CMPLR) $(FFLAGS) $(OBJT1) -o mobse.x  

USE=MPI

#--global

LOC=zfp
COMPILE=gfortran
#COMPILE=ifort

dir=../$(FSHAREDIR)
 
ifeq ($(USE),MPI)
FC=mpif90
else
FC=gfortran
endif
outdir:=$(dir)glib/mobse/
FLAG= -O3 -fopenmp -fbacktrace -std=legacy \
	-lmycom2 -fcheck=all -Waliasing  -fbounds-check   -cpp  -malign-double


OUT=libmobse.a

OBJS=comenv.o corerd.o deltat.o dgcore.o evolve.o gntage.o \
hrdiag.o instar.o kick.o mix.o mlwind.o mrenv.o ran3.o rl.o \
star.o zcnsts.o zfuncs.o pisn.o eddington.o fallback.o

hfiles = ../input/const_mobse.h ../input/zdata.h 

#---make target
all: $(OUT) 

$(OUT):$(OBJS)
	ar -rc $(OUT) $(OBJS)
	cp $(OUT)  $(outdir)
	cp $(hfiles) $(outdir)
#	sudo cp $(OUT) /usr/local/mylib/
#	sudo cp $(MODS) /usr/local/mylib/



%.o:%.f90
#   $(FC) -c -fpic $<
	$(FC) -c $(FLAG) $<

%.o:%.for
#   $(FC) -c -fpic $<
	$(FC) -c $(FLAG) $<

%.o:%.f
#   $(FC) -c -fpic $<
	$(FC) -c $(FLAG) $<

.PHONY:clean
clean:
	rm -f *.o
	rm -f *.mod
	rm -f *.a
	rm -f *.gcda
	rm -f *.gcno

