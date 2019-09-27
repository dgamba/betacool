include ../Makefile.config

INC     += -I.
LIBS     =

include ./Objects

all : $(OBJS) 
	if [ ! -d ../lib ] ; then mkdir ../lib ; fi;
	cp *.h ../lib
	$(DLD) $(DLDFLAGS) -o ../lib/libbetacool.a $(OBJS) $(LIBS)

	$(LD) -o ../dat/betacool  $(LDFLAGS) $(INC) $(OBJS) $(LIBS)

#	$(LD) $(LDFLAGS) $(INC)-o ../lib/betacool ./betacool.cpp -L../lib -lBetacool $(LIBS)

clean:
	rm -f $(OBJS)
#	rm -rf ../include/

../lib/obj/%.o : %.cpp
	$(CC) $(CCFLAGS) $(INC) -c $< -o $@;

