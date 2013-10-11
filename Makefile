#
#   NMTR main directory
#

L_ARCH   = $(ARCH)
LIB_NAME = d-$(L_ARCH).a

OPTFLAGS = -O

CFLAGS   = $(OPTFLAGS) 
FFLAGS   = $(OPTFLAGS)

# Libraries.

NMTR    = src/nmtr/$(LIB_NAME) 
LAPACK  = src/lapack/$(LIB_NAME)
BLAS    = src/blas/$(LIB_NAME)
TPROBS  = src/tprobs/$(LIB_NAME)
UTILS   = src/utils/$(LIB_NAME)

LIBS = $(NMTR) $(LAPACK) $(BLAS) $(TPROBS) $(UTILS) 

install: libs exec

libs: nmtrlib lapack blas tprobs utils

nmtrlib: 
	cd src/nmtr; make

lapack:
	cd src/lapack; make

blas: 
	cd src/blas; make

tprobs:
	cd src/tprobs; make 

utils: 
	cd src/utils; make

# MINPACK-2 Newton's method for (dense) unconstrained optimization

exec : dmain.o $(LIBS)
	$(FC) $(FFLAGS) -o nmtr dmain.o $(LIBS)

clean:
	cd src/nmtr;     make clean
	cd src/lapack;   make clean
	cd src/blas;     make clean
	cd src/tprobs;   make clean
	cd src/utils;    make clean
	-rm *.o

.c.o:
	$(CC) $(CFLAGS) -c $*.c
.f.o:
	$(FC) $(FFLAGS) -c $*.f
