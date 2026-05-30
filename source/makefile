SUBDIRS = com mobse/src gnc main

export COMPILE=gfortran
export FSHAREDIR=../fshare/
export use_scope=0
# 主目标
all: $(SUBDIRS)

GNCOBJS=com mobse 

gnc: $(GNCOBJS)

main: gnc 

$(SUBDIRS): 
	@echo "Compiling ... $@ "
	$(MAKE) -C $@ $(MAKECMDGOALS)

clean: $(SUBDIRS)
	@echo "clear finished"
	find ./fshare/ -name "*.o" -exec rm -f {} \;
	find ./fshare/ -name "*.mod" -exec rm -f {} \;
	find ./fshare/ -name "*.a" -exec rm -f {} \;

.PHONY: all clean $(SUBDIRS)