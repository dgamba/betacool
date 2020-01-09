# the target output file
TARGET   = Betacool_UNIX

# where to do the compilation
BUILDDIR = Build/UNIX/
# where the source files are
SRCDIR   = src/

# the compiler to use
CC       = g++

# compiler flags:
#  -g    adds debugging information to the executable file
#  -Wall turns on most, but not all, compiler warnings
#  -threads
#  -pthread
#  -fpermissive maybe not so nice.. but many things still to be improved in the code
CCFLAGS   = -g -Wall -fopenmp -O2 -fpermissive -Wwrite-strings

# the linker to use and its options
LD       = ld
LDFLAGS  = 

# for creating a static library:
DLD      = ar
DLDFLAGS = rcs

# if to include something more and/or add libraries for compilation
INC     += -I.
LIBS     =

# betacool is made of several small objects:
OBJS += Betacool.o \
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
	$(CC) $(CCFLAGS) $(INC) $(BUILDDIR)*.o -o $@

clean:
	@echo "Removing directory $(BUILDDIR)" 
	rm -rf $(BUILDDIR)
	rm $(TARGET)

$(BUILDDIR)%.o : $(SRCDIR)%.cpp
	@echo "Compiling $< into $@" 
	if [ ! -d $(BUILDDIR) ] ;     then mkdir $(BUILDDIR) ; fi;
	$(CC) $(CCFLAGS) $(INC) -c $< -o $@;


