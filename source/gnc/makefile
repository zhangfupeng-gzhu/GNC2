USE=MPI

#--global 
COMPILE=gfortran 
hdf5dirlib=$(HDF5LIBDIR) 
hdf5dirinc=$(HDF5INCDIR)
performance_analyze=0
use_hdf5=1
use_mobse=1
use_bse=0 

dir=$(FSHAREDIR)
 
  
ifeq ($(USE),MPI)
FC=mpif90
else
FC=gfortran
endif
libdir:=$(dir)glib/com/
outdir:=$(dir)glib/gnc/
mobsedir:=$(dir)glib/mobse/
FLAG= -O3 -fopenmp -fbacktrace \
		-L$(libdir) -I$(libdir) -lmycom2 -fcheck=all -Waliasing -Wuninitialized -fbounds-check   -cpp \
		-fimplicit-none -ffast-math -malign-double
ifeq ($(use_mobse),1)
	FLAG+=-L$(mobsedir) -I$(mobsedir) -lmobse -DMOBSE
endif
ifeq ($(use_hdf5),1)
	FLAG+= -DHDF5
	FLAG+= -I$(hdf5dirinc) -L$(hdf5dirlib) -lhdf5_fortran -lhdf5hl_fortran -lhdf5 -lhdf5_hl
endif 
ifeq ($(performance_analyze), 1)
	FLAG+= -pg --param=max-vartrack-size=60000000 -g
endif
	ifeq ($(use_bse),1)
	FLAG+= -I$(bsedir) -L$(bsedir) -lbse -DBSE
endif
ifeq ($(use_sse),1)
	FLAG+= -I$(bsedir) -L$(bsedir) -lsse
endif

ifeq ($(performance_analyze), 2)
	FLAG+= -fprofile-arcs -ftest-coverage
endif
 

OUT=libgnc2.a

OBJS= md_binary.o md_particle_sample.o md_stellar_history.o 
OBJS+=chain_pointer_type.o chain_type.o   coeff.o  dms.o  modules.o 
ifeq ($(use_mobse), 1)
	OBJS+=mobse.o
endif 
OBJS+= md_main_gw.o  mobse_stellar.o 
OBJS+= io_txt.o  ini.o ini_single.o  Sigma_funcs.o Sigma_funcs_aniso.o 
OBJS+= EJ_funcs_ba16.o get_dc.o get_dc_2d.o
OBJS+= ac_ec_evl_single.o
OBJS+= sams_sub.o sams_sub_single.o clone.o 
OBJS+=  stellar_obj.o  sams_sts_single.o MPI_comu.o
OBJS+= read_ini_single.o  io_bin.o gen_ge.o
OBJS+=  io_hdf5.o 
OBJS+= evl_mbh_mass.o
OBJS+=star_pt_dm.o  star_pt_sub.o  coe_subs.o md_iso.o 
OBJS+=star_funcs_full_range.o dstr_subs.o self_consistent_pot.o
OBJS+= coenr.o star_aux_sub.o   ini_cond.o 
OBJS+=run_contrl.o  dstr_subs_2d.o md_dms_saving_data.o  
OBJS+=event_data.o md_mbh_evl_acc.o stellar_evolution.o
OBJS+=config.o ini_paras.o  update_sample_para.o md_bk_particle.o md_by_particle.o
  
OBJS+=RR_mpi.o  

MODS=  md_coeff.mod  model_basic.mod com_main_gw.mod
MODS+= md_chain_pointer.mod md_chain.mod
MODS+= md_binary.mod md_particle_sample.mod md_particle.mod 
MODS+= md_bk_species.mod md_dms.mod
MODS+=md_stellar_object.mod md_by_particle.mod   
MODS+=md_mobse_stellar_single.mod  
MODS+=  mpi_comu.mod
MODS+=md_mass_bins.mod
MODS+=  md_iso.mod md_star_pot.mod md_ini_data.mod 
MODS+= md_dms_saving_data.mod 
MODS+=md_event_datas.mod md_mbh_evl_acc.mod md_stellar_history.mod md_stellar_evolution.mod
MODS+=model_config.mod

#---make target
all: $(OUT) 

md_mbh_evl_acc.o: stellar_obj.o md_particle_sample.o

md_particle_collision.o: md_bk_particle.o stellar_obj.o modules.o 

md_particle_sample.o : md_binary.o  md_stellar_history.o

md_stellar_history.o: mobse_stellar.o
 

chain_pointer_type.o stellar_evolution.o:  md_particle_sample.o

md_bk_particle.o: md_particle_sample.o coeff.o 

chain_type.o: chain_pointer_type.o

md_by_particle.o: chain_type.o

stellar_obj.o: coeff.o md_bk_particle.o 

dms.o: coeff.o md_by_particle.o md_bk_particle.o stellar_obj.o  chain_type.o
 
 
modules.o: chain_type.o  dms.o md_by_particle.o  mobse_stellar.o md_iso.o   chain_type.o

MPI_comu.o : modules.o  
 
md_main_gw.o: modules.o  MPI_comu.o coeff.o 

io_txt.o track.o ini.o ini_single.o Sigma_funcs.o EJ_funcs_ba16.o \
sams_sub.o sams_sub_single.o .o   clone.o  read_ini_single.o usr_pro.o  io_bin.o io_hdf5.o   sams_sts_single.o gen_ge.o\
 RR_mpi.o ac_ec_evl_single.o   star_pt_dm.o  star_pt_sub.o\
sigma_funcs_aniso.o get_dc.o get_dc_2d.o coe_subs.o evl_mbh_mass.o \
star_funcs_full_range.o dstr_subs.o dstr_subs_2d.o self_consistent_pot.o \
coenr.o   ini_cond.o star_aux_sub.o run_contrl.o \
 update_sample_para.o : md_main_gw.o 

ini_cond.o: modules.o 

read_ini_single.o ini.o:config.o

ac_ec_evl_single.o read_ini_single.o:stellar_evolution.o

stellar_evolution.o:evl_mbh_mass.o

io_hdf5.o sams_sts_single.o ini_single.o:md_dms_saving_data.o event_data.o stellar_evolution.o

evl_mbh_mass.o ini.o read_ini_single.o: md_mbh_evl_acc.o

ini_paras.o: modules.o config.o md_mbh_evl_acc.o

libs=$(libdir)/libmycom2.a   $(mobsedir)/libmobse.a 
 

$(info === Debug Info ===)
$(info FC = $(FC))
$(info FLAG = $(FLAG))
$(info libs = $(libs))
$(info outdir = $(outdir))
$(info OBJS = $(OBJS))
$(info === End Debug Info ===)

$(OUT):$(OBJS)  
	ar -rc $(OUT) $(OBJS) 
	cp $(OUT)  $(outdir)
	cp $(MODS) $(outdir) 
 


%.o:%.f90 $(libs) 
	$(FC) -c $(FLAG) $<

%.o:%.for $(libs) 
	$(FC) -c $(FLAG) $<
 
%.o:%.f $(libs) 
	$(FC) -c $(FLAG) $<

.PHONY:clean

clean:
	rm -f *.o
	rm -f *.mod
	rm -f *.a
	rm -f *.gcda
	rm -f *.gcno

