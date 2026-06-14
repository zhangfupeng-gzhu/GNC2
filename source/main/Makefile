USE=MPI

#--global

LOC=zfp
#COMPILE=ifort
COMPILE=gfortran
#hdf5dirlib=/usr/lib/x86_64-linux-gnu/hdf5/serial/
hdf5dirlib=$(HDF5LIBDIR)
#hdf5dirinc=/usr/include/hdf5/serial/
hdf5dirinc=$(HDF5INCDIR)
performance_analyze=0
use_hdf5=1
use_lapack=0
lapacklib=$(LAPACKLIBDIR)
lapackinc=$(LAPACKINCDIR) 
use_mobse=1
dir=$(FSHAREDIR)
usestack=1

ifeq ($(USE),MPI)
FC=mpif90
else
FC=gfortran
endif
comdir:=$(dir)glib/com/
libdir:=$(dir)glib/gnc/
mobsedir:=$(dir)glib/mobse/

FLAG=-I$(libdir) -L$(libdir) -lgnc2 -fstack-arrays \
-I$(comdir) -L$(comdir) -lmycom2 -fopenmp -O3 -fcheck=all -fbacktrace -cpp   -fimplicit-none \
-malign-double

ifeq ($(use_hdf5),1)
	FLAG+= -DHDF5
	FLAG+= -I$(hdf5dirinc) -L$(hdf5dirlib) -lhdf5_fortran -lhdf5hl_fortran -lhdf5 -lhdf5_hl
endif

ifeq ($(use_lapack),1)
	FLAG+=-I$(lapackinc) -L$(lapacklib) -llapack 
endif



ifeq ($(use_mobse),1)
	FLAG+=-L$(mobsedir) -I$(mobsedir) -lmobse  -DMOBSE
endif
 

ifeq ($(performance_analyze), 1)
	FLAG+=-pg -g
endif
ifeq ($(performance_analyze), 2)
	FLAG+=-fprofile-arcs -ftest-coverage
endif


libs=$(comdir)/libmycom2.a  $(mobsedir)/libmobse.a $(libdir)/libgnc2.a 

$(info === Debug Info ===)
$(info FC = $(FC))
$(info FLAG = $(FLAG))
$(info libs = $(libs))  
$(info usemobse=$(use_mobse))
$(info === End Debug Info ===)

all: main ini  get_time_range event_analysis   

%: %.f90 $(libs) | check_libs 
	$(FC) $< -o $@ $(FLAG)
	@if [ $(usestack) -eq 1 ]; then \
	    echo "Running: execstack -s $@"; \
	    execstack -s $@ || echo "Warning: execstack -s $@ failed"; \
	fi
	
check_libs:
# 	@echo "Checking required libraries..."
	@test -f $(comdir)/libmycom2.a || (echo "Error: $(comdir)/libmycom2.a not found" && false)
	@test -f $(libdir)/libgnc2.a || (echo "Error: $(libdir)/libgnc2.a not found" && false) 
	@if [ $(use_mobse) -eq 1 ]; then test -f $(mobsedir)/libmobse.a || (echo "Error: $(mobsedir)/libmobse.a not found" && false); fi

.PHONY:clean
clean:
	rm -f *.gcda
	rm -f *.gcno	