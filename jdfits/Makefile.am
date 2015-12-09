AUTOMAKE_OPTIONS = subdir-objects

lib_LIBRARIES = libjdfits.a
include_HEADERS = src/jdfits.h

libjdfits_a_SOURCES = \
src/_jdfits.h \
src/fits.c \
src/fitsfmt.c \
src/fitsbin.c \
src/fitsio.c \
src/fitsdump.c \
src/fitswrite.c \
src/fitsmisc.c \
src/fitshdr.c \
src/fitswbin.c \
src/fitscopy.c \
src/fitsheap.c \
src/fitsuser.c \
src/julian.c \
src/btread.c \
src/makefits.c \
src/ieee.c \
src/copyfits.c 

bin_PROGRAMS = fitsdump

fitsdump_SOURCES = src/fitsdump.c
fitsdump_LDADD = libjdfits.a
