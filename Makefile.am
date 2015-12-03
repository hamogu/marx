AUTOMAKE_OPTIONS = subdir-objects
SUBDIRS = jdmath jdfits pfile

#   When a file to be included is generated during the build and not part of a distribution tarball, its location is under $(builddir), not under $(srcdir). This matters especially for packages that use header files placed in sub-directories and want to allow builds outside the source tree (see VPATH Builds). In that case we recommend to use a pair of -I options, such as, e.g., ‘-Isome/subdir -I$(srcdir)/some/subdir’ or ‘-I$(top_builddir)/some/subdir -I$(top_srcdir)/some/subdir’. Note that the reference to the build tree should come before the reference to the source tree, so that accidentally leftover generated files in the source directory are ignored.

AM_CPPFLAGS = \
-I$(top_builddir)/jdmath/src -I$(top_srcdir)/jdmath/src \
-I$(top_builddir)/jdfits/src -I$(top_srcdir)/jdfits/src \
-I$(top_builddir)/pfile/src -I$(top_srcdir)/pfile/src


include_HEADERS = libsrc/marx.h

lib_LIBRARIES = libmarx.a
libmarx_a_SOURCES = \
libsrc/argcargv.h \
libsrc/argcargv.c \
libsrc/marx-feat.h \
libsrc/_marx.h \
libsrc/caldb.c \
libsrc/photon.c \
libsrc/acis_geom.c \
libsrc/acis_subpix.c \
libsrc/acis-s.c \
libsrc/acis-i.c \
libsrc/marxio.c \
libsrc/hrc-i.c \
libsrc/hrc_i_geom.c \
libsrc/hrc-s.c \
libsrc/hrc_s_geom.c \
libsrc/hrcblur.c \
libsrc/acis_fef.c \
libsrc/aciscontam.c \
libsrc/spectrum.c \
libsrc/prob.c \
libsrc/data.c \
libsrc/mirror.c \
libsrc/diffract.c \
libsrc/mblur.c \
libsrc/detector.c \
libsrc/detpix.c \
libsrc/wfold.c \
libsrc/trans.c \
libsrc/marxerr.c \
libsrc/drake.c \
libsrc/reflect.c \
libsrc/hrma.c \
libsrc/source.c \
libsrc/source.def \
libsrc/util.c \
libsrc/eamirror.c \
libsrc/grating.c \
libsrc/pixlib.c \
libsrc/rdb.c \
libsrc/dither.c \
libsrc/s-saosac.c \
libsrc/s-image.c \
libsrc/s-point.c \
libsrc/s-gauss.c \
libsrc/s-beta.c \
libsrc/s-line.c \
libsrc/s-rayfile.c \
libsrc/s-disk.c \
libsrc/s-user.c \
libsrc/ixoccd.h \
libsrc/ixocatgs.c \
libsrc/ixoccd.c \
libsrc/ixoxms.c \
libsrc/ixomirror.c \
libsrc/ffield.c \
libsrc/qefuns.c


## TBD: Should be used during install to reformat data file. 
#libsrc/fix_ee.c 
# Not sure this is used at all. Can it be removed?
#libsrc/dump_ee.c 


# TBD: These are used for interpreted marx in SLANG, I assume:
#libsrc/interp.c
#libsrc/disperse.c

# TBD: What are those good for?
#libsrc/parm.h
#libsrc/parm.c


# TBD: Can turn this into marx_LDADD, marx2fits_LDADD etc.?
# Or does every program need all libraries?
LDADD = jdmath/libjdmath.a jdfits/libjdfits.a pfile/libpfile.a libmarx.a



bin_PROGRAMS = marx marx2dpde marx2fits marx2img marxasp marxcat marxpileup

marx_SOURCES = libsrc/_marx.h src/marx.c
marx2dpde_SOURCES = src/marx2dpde.c
marx2fits_SOURCES = libsrc/marx-feat.h src/marx2fits.c
marx2img_SOURCES = src/marx2img.c
marxasp_SOURCES = src/marxasp.c
marxcat_SOURCES = src/marxcat.c
marxpileup_SOURCES = \
src/henke.h \
src/henke.c \
src/acissim.h \
src/acissim.c \
src/marxpileup.c

# TBD: Missing: bestfocus detinfo geff geffdump
# TBD: are rfl and sky2chip ever used? They are not installed...
# Maybe I can remove them?
