# the target output file
TARGET   = Betacool

# where to do the compilation
BUILDDIR = Build/UNIX/
# where the source files are
SRCDIR   = src/

# the compiler to use
CXX      = g++

# compiler flags:
#  -g -O0   adds debugging information to the executable file, without optimisation
#  -Wall turns on most, but not all, compiler warnings
#  -threads
#  -pthread
#  -fpermissive maybe not so nice.. but many things still to be improved in the code
#  -Wwrite-strings : warning on using strings attached to char*
#  -Wunused-value : warning if some defined values are not used...
CXXFLAGS  = -Wall -fopenmp -Ofast -fpermissive -Wno-write-strings -Wno-unused-value

# linker options
# Ideally we would like to compile a static binary, but on MacOSX is not so easy
uname_s := $(shell uname -s)
ifeq ($(uname_s), Darwin)
  LDFLAGS = -fopenmp
else
  LDFLAGS = -static -fopenmp
endif


# if to include something more and/or add libraries for compilation
INC     += -I.
LIBS     =

# betacool is made of several small objects:
OBJS  = Betacool.o \
        bolideU.o \
        bpData.o \
        bpIBS.o \
        bpTune.o \
        dataU.o \
        doubleU.o \
        matrixU.o \
        pellets.o \
        stdafx.o \
        vectorU.o \
        warning.o \
        xBeam.o \
        xBucket.o \
        xDistributor.o \
        xDraw.o \
        xDynamic.o \
        xEbeam.o \
        xEcool.o \
        xEffect.o \
        xForce.o \
        xHiroshi.o \
        xIBS.o \
        xLibrary.o \
        xOptics.o \
        xPowell.o \
        xRestgas.o \
        xRing.o \
        xRunge.o \
        xStoch.o \
        xTarget.o


all: $(TARGET)

$(TARGET) : $(addprefix $(BUILDDIR), $(OBJS))
	@echo "Linking $(TARGET)" 
	$(CXX) $(BUILDDIR)*.o $(LDFLAGS) -o $@

clean:
	@echo "Removing directory $(BUILDDIR) and $(TARGET)" 
	rm -rf $(BUILDDIR) $(TARGET)

$(BUILDDIR)%.o : $(SRCDIR)%.cpp
	@echo "Compiling $< into $@" 
	@mkdir -p $(BUILDDIR)
	$(CXX) $(CXXFLAGS) $(INC) -c $< -o $@;


