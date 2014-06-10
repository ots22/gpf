.SUFFIXES: 

vpath %.f90 src test
FC = gfortran-4.8
MODDIR = obj
F90FLAGS = -std=f2008 -Wall -Wextra -O3 -I /home/raid/ots22/include -I $(MODDIR) -J $(MODDIR)
F90LINKFLAGS = -lblas -llapack -lnlopt -L /home/raid/ots22/lib

TESTEXE = cov_test gp_test gp_test2 inv_test logdet_test optim_test

TESTEXE := $(addprefix bin/test/,$(TESTEXE))

RUNTESTS = cov_test gp_test gp_test2.sh inv_test logdet_test optim_test

RUNTESTS := $(addprefix bin/test/,$(RUNTESTS))
TESTOUTPUT = $(addsuffix .completed, $(RUNTESTS))

# could have additional things in tests directory
#TESTOBJ = gp_example.o
#TESTOBJ := $(addprefix obj/,$(TESTOBJ))

SRC = $(wildcard src/*.f90)
OBJ = $(SRC:.f90=.o)
OBJ := $(patsubst src/%,obj/%,$(OBJ))

.PHONY: all clean test
all: depend $(TESTEXE) lib/libgpf.a

lib/libgpf.a: $(OBJ)
	mkdir -p lib
	ar rcs lib/libgpf.a $(OBJ)

$(TESTEXE): bin/test/% : $(OBJ) $(TESTOBJ) obj/%.o depend
	mkdir -p bin/test
	$(FC) $(F90FLAGS) $(filter %.o, $^) -o $@ $(F90LINKFLAGS)
	cd bin/test; ln -sf ../../test/*.sh .

depend:
	mkdir -p obj
	makedepf90 -W -m"obj/%m.mod" -b obj src/*.f90 test/*.f90 > depend

-include depend

obj/%.o : %.f90
	mkdir -p obj
	$(FC) $(F90FLAGS) $^ -c -o $@

obj/%.mod:
	mkdir -p obj
	$(FC) $(F90FLAGS) -fsyntax-only $(filter %.f90, $^)

test: $(TESTOUTPUT)

$(TESTOUTPUT) : %.completed : %
	cd bin/test; ../../$< && touch ../../$@

clean:
	-rm -rf obj bin lib depend
