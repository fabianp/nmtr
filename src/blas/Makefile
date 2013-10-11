#
# BLAS directory.
#

L_ARCH   = $(ARCH)
LIB_NAME = d-$(L_ARCH).a

OPTFLAGS = -O

CFLAGS = $(OPTFLAGS) 
FFLAGS = $(OPTFLAGS)

AR = ar -rcv

FILES = dasum.o  dcopy.o  dgemm.o  dnrm2.o  dtrmv.o  dtrsv.o \
        daxpy.o  ddot.o   dgemv.o  dscal.o  dsyrk.o  dtrsm.o

blas: $(FILES) 
	$(AR) $(LIB_NAME) $(FILES)  

clean:
	- rm -f *.o
	- rm -f *.a

.c.o:
	$(CC) $(CFLAGS) -c $*.c
.f.o:
	$(FC) $(FFLAGS) -c $*.f


