# TBD: include targets or subdir for demo and test.

AUTOMAKE_OPTIONS = subdir-objects

include_HEADERS = src/jdmath.h

lib_LIBRARIES = libjdmath.a
libjdmath_a_SOURCES = \
src/bilinear.c \
src/fftn.c \
src/rotate.c \
src/binaryio.c \
src/dinterpo.c \
src/dreadcol.c \
src/dynmem.c \
src/finterpo.c \
src/freadcol.c \
src/_jdmath.h \
src/misc.c \
src/dsort.c \
src/fsort.c \
src/vector.c \
src/hist.c \
src/random.c \
src/simplex.c \
src/complex.c \
src/qroot.c \
src/gaussrnd.c \
src/gnuplot.c \
src/init.c \
src/polyint.c \
src/loggrid.c \
src/bdata.c \
src/nan.c \
src/erf.c \
src/dgelim.c \
src/dblas.c \
src/crout.c \
src/gamma.c \
src/bisect.c \
src/machar.c

# This files is used in "include".
# List it here ot make sure it's include in the source distribution.
EXTRA_libjdmath_a_SOURCES = \
src/fftn-inc.c

# Not used
#@src/levmrq.c
#@src/romberg.c
#@src/oromberg.c
#@src/gaussj.c
