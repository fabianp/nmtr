#
# TPROBS directory.
#

L_ARCH   = $(ARCH)
LIB_NAME = d-$(L_ARCH).a

OPTFLAGS = -O

CFLAGS = $(OPTFLAGS) 
FFLAGS = $(OPTFLAGS)

AR = ar -rcv

FILES = deptfg.o  dgl1hs.o  dgl2hs.o  dminhs.o  dmsahs.o  dodcps.o  dsscfg.o \
        depths.o  dgl2co.o  dljcfg.o  dmsabc.o  dodcfg.o  dpjbfg.o  dsschs.o \
        dgl1fg.o  dgl2fg.o  dminfg.o  dmsafg.o  dodchs.o  dpjbhs.o

tprobs: $(FILES) 
	$(AR) $(LIB_NAME) $(FILES)  

clean:
	- rm -f *.o
	- rm -f *.a

.c.o:
	$(CC) $(CFLAGS) -c $*.c
.f.o:
	$(FC) $(FFLAGS) -c $*.f
