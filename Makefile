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

