

SRCS := $(shell ls *.f)
OBJECTS = $(SRCS:.f=.o)
TARGET = libsminpack.a

# reset make's default suffix list for implicit rules, set our own
.SUFFIXES :
.SUFFIXES : .f90 .o .d .mod

all: libsminpack.a

$(TARGET): $(OBJECTS)
	ar cr $@ $(OBJECTS)
	ranlib $@

# implicit rule for compiling
%.o : %.f
	$(FORTRANC) -c $(CFLAGS) $<

clean :
	rm -f *.o $(TARGET)

