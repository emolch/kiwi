# makefile for kiwi core tools

SHELL = /bin/sh
MAKEDEPEND = ./fdepend.pl -g -d -i hdf5.mod

# compiler
FORTRANC := g95

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

# compiler and linker flag defaults
CFLAGS =  $(INCMSEED) $(INCHDF) $(INCSAC) $(INCFFTW)
LDFLAGS =  $(LIBMSEED) $(LIBSAC)  $(LIBHDF) $(LIBSMINPACK) $(LIBFFTW)

# use the Makefile.local if you want to append to FORTRANC, CFLAGS or LDFLAGS
# some abbreviations that can be used to append to cflags
# put CFLAGS += $(CDEBUG) or CFLAGS += $(CFAST) to Makefile.local
CFAST = -O3
CDEBUG_IFORT = -g -warn all -ftrapuv -debug all
CDEBUG_G95 = -g -Wall  -fbounds-check -ftrace=full
CDEBUG = $(if $(filter ifort, $(FORTRANC)), \
    $(CDEBUG_IFORT), \
    $(CDEBUG_G95) )

-include Makefile.local

# communicate compiler settings to submake (for sminpack)
export CFLAGS FORTRANC

SRCS := $(shell ls *.f90)
OBJECTS = $(SRCS:.f90=.o)

TARGETS = eulermt source_info minimizer \
          gfdb_build gfdb_extract gfdb_redeploy gfdb_info gfdb_specialextract \
          gfdb_build_ahfull differential_azidist eikonal_benchmark crust

TESTS_SRCS := $(shell ls test_*.f90)
TESTS = $(TESTS_SRCS:.f90=)

OBJECTS = better_varying_string.o varying_string_getarg.o constants.o util.o \
          geometry.o crust2x2.o \
          heap.o eikonal.o euler.o orthodrome.o discrete_source.o gfdb.o receiver.o \
          piecewise_linear_function.o parameterized_source.o source.o source_moment_tensor.o \
          source_bilat.o source_circular.o source_point_lp.o source_eikonal.o \
          source_all.o sparse_trace.o \
          seismogram_io.o unit.o seismogram.o read_line.o read_table.o \
          comparator.o elseis_oo.o elseis.o differentiation.o integration.o interpolation.o \
          minimizer_engine.o
          

.PHONY : clean clean-deps tests targets all check install uninstall

# reset make's default suffix list for implicit rules, set our own
.SUFFIXES :
.SUFFIXES : .f90 .o .d .mod

all : check targets 

$(TARGETS) $(TESTS) : .sminpackdone .mseedsimple

.sminpackdone :
	$(MAKE) -C sminpack/ && touch .sminpackdone

.mseedsimple :
	$(MAKE) -C mseed/ && touch .mseedsimple

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

check : tests
	@for t in $(TESTS); do ./$$t ; done

seismosizer : $(OBJECTS) seismosizer.o
	$(FORTRANC) $(OBJECTS) $@.o $(LDFLAGS) -o $@

minimizer : $(OBJECTS) minimizer.o
	$(FORTRANC) -o $@ $(OBJECTS) $@.o $(LDFLAGS)

eulermt : $(OBJECTS) eulermt.o
	$(FORTRANC) $(OBJECTS) $@.o $(LDFLAGS) -o $@

source_info : $(OBJECTS) source_info.o
	$(FORTRANC) $(OBJECTS) $@.o $(LDFLAGS) -o $@

gfdb_build : $(OBJECTS) gfdb_build.o
	$(FORTRANC) $(OBJECTS) $@.o $(LDFLAGS) -o $@

gfdb_build_ahfull : $(OBJECTS) gfdb_build_ahfull.o
	$(FORTRANC) $(OBJECTS) $@.o $(LDFLAGS) -o $@

gfdb_extract : $(OBJECTS) gfdb_extract.o
	$(FORTRANC) $(OBJECTS) $@.o $(LDFLAGS) -o $@

gfdb_specialextract : $(OBJECTS) gfdb_specialextract.o
	$(FORTRANC) $(OBJECTS) $@.o $(LDFLAGS) -o $@

gfdb_redeploy : $(OBJECTS) gfdb_redeploy.o
	$(FORTRANC) $(OBJECTS) $@.o $(LDFLAGS) -o $@

gfdb_info : $(OBJECTS) gfdb_info.o
	$(FORTRANC) $(OBJECTS) $@.o $(LDFLAGS) -o $@

differential_azidist : $(OBJECTS) differential_azidist.o
	$(FORTRANC) $(OBJECTS) $@.o $(LDFLAGS) -o $@

eikonal_benchmark : $(OBJECTS) eikonal_benchmark.o
	$(FORTRANC) $(OBJECTS) $@.o $(LDFLAGS) -o $@

crust : $(OBJECTS) crust.o
	$(FORTRANC) $(OBJECTS) $@.o $(LDFLAGS) -o $@


# simple implicit link rule for the tests
test_% : $(OBJECTS) test_%.o
	$(FORTRANC) $(OBJECTS) $@.o $(LDFLAGS) -o $@

# implicit rule for generating depfiles
%.d : %.f90
	$(MAKEDEPEND) $<
    

# implicit rule for compiling
%.o : %.f90
	$(FORTRANC) -c $(CFLAGS) $<


# include auto-created dependencies
-include $(SRCS:.f90=.d)

clean :
	rm -f *.o *.mod $(TESTS) $(TARGETS) .sminpackdone .mseedsimple
	$(MAKE) -C sminpack/ clean
	$(MAKE) -C mseed/ clean
    
clean-deps : clean
	rm -f *.d
