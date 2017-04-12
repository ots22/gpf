.SUFFIXES: 

vpath %.f90 src test include
vpath %.f src test include

MODDIR = obj
F90FLAGS = -std=f2008 -frealloc-lhs -Wall -Wextra -Wno-compare-reals -Wno-unused-dummy-argument -O3 -I $(MODDIR) -J $(MODDIR) -I include
F90LINKFLAGS = -lblas -llapack -lnlopt

SRC = $(wildcard src/*.f90)
OBJ = $(SRC:.f90=.o)
OBJ := $(patsubst src/%,obj/%,$(OBJ))

TESTSRC = $(wildcard test/*.f90)
TESTEXE = $(basename $(notdir $(TESTSRC)))
TESTEXE := $(addprefix bin/test/,$(TESTEXE))

TESTSCRIPTS = $(notdir $(wildcard test/*.sh))
TESTSCRIPTS := $(addprefix bin/test/,$(TESTSCRIPTS))

# The tests to run.  If a script with the same basename as a test
# executable exists, run this instead.
RUNTESTS = $(filter-out $(basename $(TESTSCRIPTS)),$(TESTEXE)) $(TESTSCRIPTS)

TESTOUTPUT = $(addsuffix .completed, $(RUNTESTS))


.PHONY: all
all: depend $(TESTEXE) lib/libgpf.a

clean:
	-rm -rf obj bin lib depend
.PHONY: clean

lib/libgpf.a: $(OBJ)
	mkdir -p lib
	ar rcs lib/libgpf.a $(OBJ)

$(TESTEXE): bin/test/% : $(OBJ) obj/%.o depend
	mkdir -p bin/test
	$(FC) $(F90FLAGS) $(filter %.o, $^) -o $@ $(F90LINKFLAGS)

$(TESTSCRIPTS) : bin/test/%.sh : test/%.sh
	cd bin/test && ln -sf ../../$< .

# makedepf90 doesn't know about some intrinsic modules
# (e.g. iso_c_binding): the -u flag supresses the warning about not
# being able to find it.
# 
# It may a warning about not finding nlopt.f, despite it being in the
# include path of the compiler
depend:
	mkdir -p obj
	makedepf90 -u iso_c_binding -I include -W -m"obj/%m.modstamp" -b obj src/*.f90 test/*.f90 > depend

-include depend

obj/%.o : %.f90
	mkdir -p obj
	$(FC) $(F90FLAGS) $(filter %.f90 %.F90, $^) -c -o $@

# gfortran (usefully) doesn't update mod files (or their timestamp) if
# nothing would change.  This can stop an unnecessary compilation
# cascade, but then the rule for a mod file deemed out of date
# (because a .f90 file changed) will then always be run (since it
# remains out of date).  One solution is a 'witness file' as a target
# (here with extension modstamp), created as witness to the recipe
# being run.
obj/%.modstamp:
	mkdir -p obj
	$(FC) $(F90FLAGS) -fsyntax-only $(filter %.f90, $^) && touch $@

.PHONY: test
test: $(TESTOUTPUT)

# convenient to allow tests to be re-run (and don't actually need to
# output '.completed' files)
.PHONY : $(TESTOUTPUT)
$(TESTOUTPUT) : %.completed : %
	cd bin/test && ./$(notdir $<)
