#
#                --- Makefile for Kiwi Core Tools ---
#
#          Override default settings in file 'Makefile.local'
#                   (see 'Makefile.local.example')
#

#### Basic settings ------------------------------------------------------------

SHELL = /bin/sh
MAKEDEPEND = ./fdepend.pl -g -d -i hdf5.mod -i omp_lib.mod
OBJDEPEND = ./objdepend.pl
FORTRANC := gfortran


#### Preset selection ----------------------------------------------------------

# Can use 'fast' or 'debug' with gfortran, g95 and ifort. You may have to create
# a custom preset when using a different compiler (see below and example in
# Makefile.local.example).

PRESET := fast


#### Installation prefixes -----------------------------------------------------

prefix = /usr/local
datarootdir = $(prefix)/share
datadir = $(datarootdir)
exec_prefix = $(prefix)
bindir = $(exec_prefix)/bin


#### Default library includes and linker settings ------------------------------

INCDUMMYOMP = -Idummy_omp_lib
LIBDUMMYOMP = dummy_omp_lib/omp_lib.o

INCMSEED = 
LIBMSEED = mseed/mseed_simple.o -lmseed

INCHDF = 
LIBHDF = -lhdf5_fortran -lhdf5 -lz

INCSAC = 
LIBSAC = -lsacio

INCFFTW = 
LIBFFTW = -lfftw3f 

LIBSMINPACK = -Lsminpack -lsminpack


#### Compiler and linker flag defaults -----------------------------------------

CFLAGS =  $(INCMSEED) $(INCHDF) $(INCSAC) $(INCFFTW) \
	      $(CFLAGS_$(FORTRANC)_$(PRESET))  #
LDFLAGS = $(LIBMSEED) $(LIBSAC)  $(LIBHDF) $(LIBSMINPACK) $(LIBFFTW) \
          $(LDFLAGS_$(FORTRANC)_$(PRESET)) #


#### Compiler specific presets  ------------------------------------------------

CFLAGS_ifort_fast   = -openmp 
LDFLAGS_ifort_fast  = -openmp

CFLAGS_ifort_debug  = -openmp -g -warn all -ftrapuv -debug all
LDFLAGS_ifort_debug = -openmp

CFLAGS_g95_fast     = $(INCDUMMYOMP) -O3
LDFLAGS_g95_fast    = $(LIBDUMMYOMP)

CFLAGS_g95_debug    = $(INCDUMMYOMP) -g -warn all -ftrapuv -debug all
LDFLAGS_g95_debug   = $(LIBDUMMYOMP) 

CFLAGS_gfortran_fast   = -fopenmp -O3
LDFLAGS_gfortran_fast  = -fopenmp

CFLAGS_gfortran_debug  = -fopenmp -g -Wall
LDFLAGS_gfortran_debug = -fopenmp


#### ---------------------------------------------------------------------------

MACHINE := $(shell ./hostinfo.pl --machine)
OS := $(shell ./hostinfo.pl --os)

-include Makefile.local

# communicate compiler settings to submake (for sminpack)
export FORTRANC

SRCS := $(shell ls *.f90)

TARGETS := eulermt source_info minimizer gfdb_build gfdb_extract gfdb_redeploy \
		  gfdb_info gfdb_specialextract gfdb_build_ahfull differential_azidist \
		  eikonal_benchmark crust ahfull

TESTS_SRCS := $(shell ls test_*.f90)
TESTS = $(TESTS_SRCS:.f90=)

.PHONY : clean clean-deps tests targets all check install uninstall

# reset make's default suffix list for implicit rules, set our own
.SUFFIXES :
.SUFFIXES : .f90 .o .d .mod

all : targets 

$(TARGETS) $(TESTS) : .sminpackdone .mseedsimple .dummyomplib .dummysacio

.sminpackdone :
	$(MAKE) -C sminpack/ && touch .sminpackdone

.mseedsimple :
	$(MAKE) -C mseed/ && touch .mseedsimple

.dummyomplib :
	$(MAKE) -C dummy_omp_lib/ && touch .dummyomplib

.dummysacio :
	$(MAKE) -C dummy_sacio/ && touch .dummysacio

targets : $(TARGETS)


install : targets
	install -d $(bindir)
	install $(TARGETS) $(bindir)
	install -d $(datadir)/kiwi
	for f in `find aux -type d -and -not -path '*/.svn*'` ; do \
	    install -d $(datadir)/kiwi/$$f ; done
	for f in `find aux -type f -and -not -path '*/.svn/*'` ; do \
	    install  $$f $(datadir)/kiwi/$$f ; done

	@echo 
	@echo '-----------------------------------------------------------------------'
	@echo '  Installation complete.'
	@echo '  Please adjust your environment variables:'
	@echo
	@echo '   * PATH should contain:'
	@echo '      ' $(bindir)
	@echo
	@echo '   * KIWI_HOME should be set to:'
	@echo '      ' $(datadir)/kiwi
	@echo '-----------------------------------------------------------------------'

uninstall :
	rm -rf -d $(datadir)/kiwi
	cd $(bindir) ; rm -f $(TARGETS)

tests : $(TESTS)

printvars :
	@echo FORTRANC = $(FORTRANC)
	@echo CFLAGS = $(CFLAGS)
	@echo LDFLAGS = $(LDFLAGS)

check : tests
	@for t in $(TESTS); do ./$$t ; done


$(TARGETS) $(TESTS) : 
	$(FORTRANC) $(filter %.o,$^) $(OMPLIB_$(FORTRANC)) $(LDFLAGS) -o $@


# implicit rules for generating depfiles
%.d : %.f90
	@$(MAKEDEPEND) $<
	@echo determining dependencies for $<...

progobjects.do : $(SRCS:.f90=.d)
	@$(OBJDEPEND) $(TARGETS) $(TESTS) -- $(SRCS:.f90=.d) > $@
	@echo determining dependencies for executables...

# implicit rule for compiling
%.o : %.f90
	$(FORTRANC) -c $(CFLAGS) $<


clean :
	rm -f *.o *.mod $(TESTS) $(TARGETS) .sminpackdone .mseedsimple .dummysacio .dummyomplib dummy_omp_lib/omp_lib.o dummy_omp_lib/omp_lib.mod
	$(MAKE) -C sminpack/ clean
	$(MAKE) -C mseed/ clean
	$(MAKE) -C dummy_omp_lib/ clean
	$(MAKE) -C dummy_sacio/ clean

    
clean-deps : clean
	rm -f *.d *.do

# include auto-created dependencies

-include progobjects.do
-include $(SRCS:.f90=.d) 
