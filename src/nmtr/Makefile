#
# NMTR directory.
#

L_ARCH   = $(ARCH)
LIB_NAME = d-$(L_ARCH).a

OPTFLAGS = -O

CFLAGS = $(OPTFLAGS) 
FFLAGS = $(OPTFLAGS)

AR = ar -rcv

FILES = destsv.o  dgqt.o  dnmtr.o

src: $(FILES) 
	$(AR) $(LIB_NAME) $(FILES)  

clean:
	- rm -f *.o
	- rm -f *.a

.c.o:
	$(CC) $(CFLAGS) -c $*.c
.f.o:
	$(FC) $(FFLAGS) -c $*.f


