#
# UTILS library.
#

L_ARCH   = $(ARCH)
LIB_NAME = d-$(L_ARCH).a

OPTFLAGS = -O

CFLAGS    = $(OPTFLAGS) -D$(ARCH)
FFLAGS    = $(OPTFLAGS)
CPPFLAGS  = -E -C -D$(ARCH)

AR = ar -rcv
CC = cc

FILES = cputime.o dpmeps.o wallclock.o

src: $(FILES)
	$(AR) $(LIB_NAME) $(FILES)

clean:
	- rm -f *.o

.SUFFIXES: .c .f .F

.c.o:
	$(CC) $(CFLAGS) -c $*.c
.f.o:
	$(FC) $(FFLAGS) -c $*.f
.F.o:
	./Fpp $*.F $(CPPFLAGS)
	$(FC) $(FFLAGS) -c $*.f
	rm $*.f


