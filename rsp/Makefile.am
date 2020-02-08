AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = -I build-aux

#   When a file to be included is generated during the build and not
#   part of a distribution tarball, its location is under $(builddir),
#   not under $(srcdir). This matters especially for packages that use
#   header files placed in sub-directories and want to allow builds
#   outside the source tree (see VPATH Builds). In that case we
#   recommend to use a pair of -I options, such as, e.g.,
#   ‘-Isome/subdir -I$(srcdir)/some/subdir’ or
#   ‘-I$(top_builddir)/some/subdir -I$(top_srcdir)/some/subdir’. Note
#   that the reference to the build tree should come before the
#   reference to the source tree, so that accidentally leftover
#   generated files in the source directory are ignored.

AM_CPPFLAGS = \
-I$(top_builddir)/../jdmath/src -I$(top_srcdir)/../jdmath/src \
-I$(top_builddir)/../jdfits/src -I$(top_srcdir)/../jdfits/src \
-I$(top_builddir)/../pfile/src -I$(top_srcdir)/../pfile/src \
-I$(top_builddir)/../libsrc -I$(top_srcdir)/../libsrc \
-I$(top_builddir)/../src -I$(top_srcdir)/../src

include_HEADERS = ../libsrc/marx.h src/cfits.h ../libsrc/marx-feat.h

lib_LIBRARIES = libcfits.a

# TBD: Can turn this into marx_LDADD, marx2fits_LDADD etc.?
# Or does every program need all libraries?
LDADD = ../jdmath/libjdmath.a ../jdfits/libjdfits.a ../pfile/libpfile.a \
 ../libmarx.a

libcfits_a_SOURCES = \
src/cfits.h \
src/cfits.c

bin_PROGRAMS = marxrsp rspdump

marxrsp_LDADD = $(LDADD) libcfits.a
rspdump_LDADD = $(marxrsp_LDADD)

marxrsp_SOURCES = \
src/argcargv.h \
src/argcargv.c \
src/rsp.c \
src/marxrsp.c

rspdump_SOURCES = \
src/argcargv.h \
src/argcargv.c \
src/rsp.c \
src/rspdump.c

