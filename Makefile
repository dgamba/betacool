INC     += -I.
LIBS     =
CC       = gcc
CCFLAGS  = 
BUILDDIR = tmp

include ./Objects

all : $(OBJS) 
	if [ ! -d build/MACx64 ] ; then mkdir build/MACx64 ; fi;
	if [ ! -d build/MACx64/lib ] ; then mkdir build/MACx64/lib ; fi;
	cp src/*.h build/MACx64/lib
	$(DLD) $(DLDFLAGS) -o ../lib/libbetacool.a $(OBJS) $(LIBS)

	$(LD) -o ../dat/betacool  $(LDFLAGS) $(INC) $(OBJS) $(LIBS)

clean:
	rm -f $(OBJS)
#	rm -rf ../include/

../lib/obj/%.o : %.cpp
	$(CC) $(CCFLAGS) $(INC) -c $< -o $@;

OBJS += ../lib/obj/bolideu.o \
	../lib/obj/datau.o \
	../lib/obj/doubleu.o \
	../lib/obj/matrixu.o \
	../lib/obj/stdafx.o \
	../lib/obj/vectoru.o \
	../lib/obj/xbeam.o \
	../lib/obj/xdistributor.o \
	../lib/obj/xdraw.o \
	../lib/obj/xdynamic.o \
	../lib/obj/xebeam.o \
	../lib/obj/xecool.o \
	../lib/obj/xeffect.o \
	../lib/obj/xforce.o \
	../lib/obj/xhiroshi.o \
	../lib/obj/xibs.o \
	../lib/obj/xlibrary.o \
	../lib/obj/xoptics.o \
	../lib/obj/xpowell.o \
	../lib/obj/xrestgas.o \
	../lib/obj/xring.o \
	../lib/obj/xrunge.o \
	../lib/obj/xstoch.o \
	../lib/obj/xtarget.o \
	../lib/obj/betacool.o

