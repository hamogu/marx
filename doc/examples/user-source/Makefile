# -*- sh -*-
#  Makefile for example shared objects
#

#CC_SHARED = cc -G -KPIC           # Solaris with cc
#CC_SHARED = gcc -G -fPIC          # Solaris with gcc
CC_SHARED = gcc -shared -fPIC     # Linux
RPATH_OPT = -Wl,-R
# This is where we put cfitsio, etc.  What about you?
PREFIX = /nfs/cxc/a1/$(CSR_CARCH)#

# Where is your marx-dist directory?  Here is where we put ours:
#MARXDIST_ROOT = /nfs/wiwaxia/d4/ASC/src/marx-dist
MARXDIST_ROOT = /aluche/d2/src/marx-dist

#---------------------------------------------------------------------------
#CFLAGS =
CFITSIO_ROOT = /nfs/cxc/h1/davis/sys/linux-x86_64/opt/cfitsio
CFITSIO_LIB = -L$(CFITSIO_ROOT)/lib $(RPATH_OPT)$(CFITSIO_ROOT)/lib -lcfitsio
CFITSIO_INC = -I$(CFITSIO_ROOT)/include
JDMATH_INC = -I$(MARXDIST_ROOT)/jdmath/src
JDMATH_LIB = -L$(MARXDIST_ROOT)/jdmath/src/$(ARCH)objs -ljdmath
MARX_INC = -I$(MARXDIST_ROOT)/marx/libsrc
MARX_LIB = -L$(MARXDIST_ROOT)/marx/libsrc/$(ARCH)objs -lmarx
PFILE_INC = -I$(MARXDIST_ROOT)/pfile/src
PFILE_LIB = -L$(MARXDIST_ROOT)/pfile/src/$(ARCH)objs -lpfile
SLANG_INC = -I$(PREFIX)/include
SLANG_LIB = -L$(PREFIX)/lib $(RPATH_OPT)$(PREFIX)/lib -lslang

INCS = $(JDMATH_INC) $(MARX_INC) $(PFILE_INC) $(SLANG_INC)
LIBS =

XCC_SHARED = $(CC_SHARED) $(CFLAGS) $(INCS) 
XCC_TESTING = $(CC_SHARED) $(CFLAGS) $(INCS) -DTESTING

all: slang.so pnts.so point.so image.so grid.so off-axis.so eventlist.so

point.so: point.c
	$(XCC_SHARED) point.c -o point.so $(LIBS)
pnts.so: pnts.c
	$(XCC_SHARED) pnts.c -o pnts.so $(LIBS)
image.so: image.c
	$(XCC_SHARED) $(CFITSIO_INC) image.c -o image.so $(CFITSIO_LIB) $(LIBS)
grid.so: grid.c
	$(XCC_SHARED) grid.c -o grid.so $(LIBS)
off-axis.so: off-axis.c
	$(XCC_SHARED) off-axis.c -o off-axis.so $(LIBS)

slang.so: slang.c
	$(XCC_SHARED) $(SLANG_INC) slang.c -o slang.so $(LIBS) $(SLANG_LIB) -lm

eventlist.so: eventlist.c
	$(XCC_SHARED) $(CFITSIO_INC) eventlist.c -o eventlist.so $(CFITSIO_LIB) $(MARX_LIB) $(JDMATH_LIB) $(PFILE_LIB) $(LIBS)

testimage: image.c
	$(XCC_TESTING) $(CFITSIO_INC) image.c -o testimage $(CFITSIO_LIB) $(LIBS)

clean:
	/bin/rm -f *.so
