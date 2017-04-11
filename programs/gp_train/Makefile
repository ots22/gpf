# should contain paths to NLOPT and GPF
MOD = ../../obj   /home/raid/ots22/include
LIB = ../../lib   /home/raid/ots22/lib

F90FLAGS = -std=f2008 -Wall -Wextra -O3 $(addprefix -I,$(MOD))
F90LINKFLAGS =  $(addprefix -L,$(LIB)) -lgpf -lnlopt -llapack -lblas

gp_train: gp_train.f90 ../../lib/libgpf.a
	$(FC) $(F90FLAGS) -o gp_train gp_train.f90 $(F90LINKFLAGS)
