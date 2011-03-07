# makefile for kiwi core tools

SHELL = /bin/sh
MAKEDEPEND = ./fdepend.pl -g -d -i hdf5.mod -i omp_lib.mod
OBJDEPEND = ./objdepend.pl
# compiler
FORTRANC := g95

# preset (use 'fast' or 'debug')
PROFILE := fast

# host infos
MACHINE := $(shell ./hostinfo.pl --machine)
OS := $(shell ./hostinfo.pl --os)

# installation prefix
prefix = .
datarootdir = $(prefix)/share
datadir = $(datarootdir)
exec_prefix = $(prefix)
bindir = $(exec_prefix)/bin


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
	      $(CFLAGS_$(FORTRANC)_$(PROFILE))  #
LDFLAGS = $(LIBMSEED) $(LIBSAC)  $(LIBHDF) $(LIBSMINPACK) $(LIBFFTW) \
          $(LDFLAGS_$(FORTRANC)_$(PROFILE)) #


#### Compiler specific settings ------------------------------------------------

OMPLIB_ifort        = 
CFLAGS_ifort_fast   = -openmp -fast -parallel
LDFLAGS_ifort_fast  = -openmp

CFLAGS_ifort_debug  = -openmp -g -warn all -ftrapuv -debug all
LDFLAGS_ifort_debug = -openmp


OMPLIB_g95 	        = dummy_omp_lib/omp_lib.o
CFLAGS_g95_fast     = -Idummy_omp_lib -O3
LDFLAGS_g95_fast    = 

CFLAGS_g95_debug    = -g -warn all -ftrapuv -debug all
LDFLAGS_g95_debug   = 

#### ---------------------------------------------------------------------------

-include Makefile.local

# communicate compiler settings to submake (for sminpack)
export FORTRANC

SRCS := $(shell ls *.f90)

TARGETS := eulermt source_info minimizer gfdb_build gfdb_extract gfdb_redeploy \
		  gfdb_info gfdb_specialextract gfdb_build_ahfull differential_azidist \
		  eikonal_benchmark crust

TESTS_SRCS := $(shell ls test_*.f90)
TESTS = $(TESTS_SRCS:.f90=)

.PHONY : clean clean-deps tests targets all check install uninstall

# reset make's default suffix list for implicit rules, set our own
.SUFFIXES :
.SUFFIXES : .f90 .o .d .mod

all : targets 

$(TARGETS) $(TESTS) : .sminpackdone .mseedsimple

.sminpackdone :
	$(MAKE) -C sminpack/ && touch .sminpackdone

.mseedsimple :
	$(MAKE) -C mseed/ && touch .mseedsimple

dummy_omp_lib/omp_lib.mod dummy_omp_lib/omp_lib.o : Makefile.local Makefile dummy_omp_lib/omp_lib.f90
	cd dummy_omp_lib ; $(FORTRANC) -c omp_lib.f90 -o omp_lib.o 

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
	$(FORTRANC) $^ $(OMPLIB_$(FORTRANC)) $(LDFLAGS) -o $@


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

# include auto-created dependencies

-include progobjects.do
-include $(SRCS:.f90=.d) 

clean :
	rm -f *.o *.mod $(TESTS) $(TARGETS) .sminpackdone .mseedsimple
	$(MAKE) -C sminpack/ clean
	$(MAKE) -C mseed/ clean
    
clean-deps : clean
	rm -f *.d *.do
