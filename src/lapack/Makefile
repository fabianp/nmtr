#
# LAPACK directory.
#

L_ARCH   = $(ARCH)
LIB_NAME = d-$(L_ARCH).a

OPTFLAGS = -O

CFLAGS = $(OPTFLAGS) 
FFLAGS = $(OPTFLAGS)

AR = ar -rcv

FILES = dpotf2.o dpotrf.o ilaenv.o lsame.o xerbla.o
        
blas: $(FILES) 
	$(AR) $(LIB_NAME) $(FILES)  

clean:
	- rm -f *.o
	- rm -f *.a

.c.o:
	$(CC) $(CFLAGS) -c $*.c
.f.o:
	$(FC) $(FFLAGS) -c $*.f


