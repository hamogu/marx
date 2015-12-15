AUTOMAKE_OPTIONS = subdir-objects

SUBDIRS = jdmath jdfits pfile

#   When a file to be included is generated during the build and not part of a distribution tarball, its location is under $(builddir), not under $(srcdir). This matters especially for packages that use header files placed in sub-directories and want to allow builds outside the source tree (see VPATH Builds). In that case we recommend to use a pair of -I options, such as, e.g., ‘-Isome/subdir -I$(srcdir)/some/subdir’ or ‘-I$(top_builddir)/some/subdir -I$(top_srcdir)/some/subdir’. Note that the reference to the build tree should come before the reference to the source tree, so that accidentally leftover generated files in the source directory are ignored.

AM_CPPFLAGS = \
-I$(top_builddir)/jdmath/src -I$(top_srcdir)/jdmath/src \
-I$(top_builddir)/jdfits/src -I$(top_srcdir)/jdfits/src \
-I$(top_builddir)/pfile/src -I$(top_srcdir)/pfile/src \
-I$(top_builddir)/libsrc -I$(top_srcdir)/libsrc

include_HEADERS = libsrc/marx.h rsp/cfits.h libsrc/marx-feat.h

lib_LIBRARIES = libmarx.a libcfits.a
libmarx_a_SOURCES = \
libsrc/argcargv.h \
libsrc/argcargv.c \
libsrc/_marx.h \
libsrc/chandra.h \
libsrc/acis.h \
libsrc/hrc.h \
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

libcfits_a_SOURCES = \
rsp/cfits.h \
rsp/cfits.c

# TBD: Can turn this into marx_LDADD, marx2fits_LDADD etc.?
# Or does every program need all libraries?
LDADD = jdmath/libjdmath.a jdfits/libjdfits.a pfile/libpfile.a libmarx.a

bin_PROGRAMS = marx marx2dpde marx2fits marx2img marxasp marxcat marxpileup \
               bestfocus detinfo geff geffdump marxrsp rspdump

marx_SOURCES = libsrc/_marx.h src/marx.c
marx2dpde_SOURCES = src/marx2dpde.c
marx2fits_SOURCES = src/marx2fits.c
marx2img_SOURCES = src/marx2img.c
marxasp_SOURCES = src/marxasp.c
marxcat_SOURCES = src/marxcat.c
marxpileup_SOURCES = \
src/henke.h \
src/henke.c \
src/acissim.h \
src/acissim.c \
src/marxpileup.c
bestfocus_SOURCES = src/bestfocus.c
detinfo_SOURCES = src/detinfo.c
geff_SOURCES = src/geff.c
geffdump_SOURCES = src/geffdump.c

marxrsp_LDADD = $(LDADD) libcfits.a
rspdump_LDADD = $(marxrsp_LDADD)

marxrsp_SOURCES = \
rsp/argcargv.h \
rsp/argcargv.c \
rsp/rsp.c \
rsp/marxrsp.c

rspdump_SOURCES = \
rsp/argcargv.h \
rsp/argcargv.c \
rsp/rsp.c \
rsp/rspdump.c


caldatarootdir = $(pkgdatadir)/data

parfilesdir = $(pkgdatadir)/pfiles
acisdatadir = $(caldatarootdir)/acis
caldbdatadir = $(caldatarootdir)/caldb
eamirrordatadir = $(caldatarootdir)/ea-mirror
gratingdatadir = $(caldatarootdir)/grating
hrcdatadir = $(caldatarootdir)/hrc
hrmadatadir = $(caldatarootdir)/hrma
obsdatadir = $(caldatarootdir)/obs
pixlibdatadir = $(caldatarootdir)/pixlib

dist_parfiles_DATA = \
par/marx.par \
par/marxasp.par \
par/marxpileup.par

dist_caldataroot_DATA = \
data/FORMAT.TXT \
data/README \
data/marxcaldb.par

dist_acisdata_DATA = \
data/acis/README

dist_caldbdata_DATA = \
data/caldb/acisD1997-04-17qeN0006.fits \
data/caldb/hetgD1996-11-01greffpr001N0007.fits \
data/caldb/acisD1999-07-22subpixN0001_marx.fits \
data/caldb/letgD1996-11-01greffpr001N0007.fits \
data/caldb/acisD1999-08-13contamN0009_marx.fits \
data/caldb/telD1999-07-23aimptsN0002.fits \
data/caldb/acisfef.fits \
data/caldb/telD1999-07-23geomN0006.fits \
data/caldb/fixfef.sl

dist_eamirrordata_DATA = \
data/ea-mirror/Ae_s1346.dat \
data/ea-mirror/README \
data/ea-mirror/Mirror_EE.bin \
data/ea-mirror/mirr-ee.bin

dist_gratingdata_DATA = \
data/grating/README \
data/grating/hetgmp1D1996-11-01marxsectorN0002.tbl \
data/grating/hetgmp3D1996-11-01marxsectorN0002.tbl \
data/grating/hetgmp4D1996-11-01marxsectorN0003.tbl \
data/grating/hetgmp6D1996-11-01marxsectorN0003.tbl \
data/grating/letgD1996-11-01greffMARXpr001N0004.dat \
data/grating/optical-constants.dat \
data/grating/optical-constants.xray

dist_hrcdata_DATA = \
data/hrc/README \
data/hrc/gold.xray \
data/hrc/carbon.dat \
data/hrc/carbon.xray \
data/hrc/chromium.dat \
data/hrc/chromium.xray \
data/hrc/gold.dat \
data/hrc/hrc_i_geom.txt \
data/hrc/hrc_s_geom.txt \
data/hrc/mcp_qe_i.v2.2.dat \
data/hrc/mcp_qe_s.v3.dat \
data/hrc/uvis_s1.v4.dat \
data/hrc/uvis_s2.v4.dat \
data/hrc/uvis_s34.v4.dat \
data/hrc/uvis_s56.v4.dat \
data/hrc/uvis_s56.v4.dat \
data/hrc/uvis_i.v4.dat

dist_hrmadata_DATA = \
data/hrma/EKCHDOS01.rdb \
data/hrma/EKCHDOS02.rdb \
data/hrma/EKCHDOS03.rdb \
data/hrma/EKCHDOS04.rdb \
data/hrma/EKCHDOS05.rdb \
data/hrma/EKCHDOS06.rdb \
data/hrma/README \
data/hrma/corr_1.dat \
data/hrma/corr_3.dat \
data/hrma/corr_4.dat \
data/hrma/corr_6.dat \
data/hrma/ideal.rdb \
data/hrma/iridium.dat \
data/hrma/iridium.xray \
data/hrma/iridium_2000.dat \
data/hrma/scat_h1_M.bin \
data/hrma/scat_h3_M.bin \
data/hrma/scat_h4_M.bin \
data/hrma/scat_h6_M.bin \
data/hrma/scat_p1_M.bin \
data/hrma/scat_p3_M.bin \
data/hrma/scat_p4_M.bin \
data/hrma/scat_p6_M.bin

dist_obsdata_DATA = \
data/obs/acis_i_obs.par \
data/obs/acis_s_obs.par \
data/obs/hrc_i_obs.par \
data/obs/hrc_s_obs.par

dist_pixlibdata_DATA = \
data/pixlib/pix_corner_lsi.par \
data/pixlib/pix_sim_table_flight.par

dist_doc_DATA = \
doc/README \
doc/dump.txt \
doc/parms.txt \
doc/tma.txt \
doc/aimpoint.txt \
doc/marx.html \
doc/source.tex

hrmadocdir = $(docdir)/hrma
examplesdocdir = $(docdir)/examples
usersrcdocdir = $(examplesdocdir)/user-source

dist_hrmadoc_DATA = \
doc/hrma/DG-19961223-1.txt \
doc/hrma/JPH-19960304.txt \
doc/hrma/JWA-19960603.txt \
doc/hrma/LVS-19960614.txt \
doc/hrma/PZ-19951128.txt \
doc/hrma/RJE-19961220.txt \
doc/hrma/TJG-19961007A.txt \
doc/hrma/TJG-19961007B.txt \
doc/hrma/TJG-19961007C.txt \
doc/hrma/WAP-19950929.txt \
doc/hrma/WAP-19960328.txt \
doc/hrma/WAP-19960410.txt \
doc/hrma/WAP-19961008.txt \
doc/hrma/WAP-19961108.txt \
doc/hrma/WAP-19970327.txt \
doc/hrma/fields.txt \
doc/hrma/hrma-OLD.txt

dist_usersrcdoc_DATA = \
doc/examples/user-source/Makefile \
doc/examples/user-source/README \
doc/examples/user-source/eventlist.c \
doc/examples/user-source/image.c \
doc/examples/user-source/point.c \
doc/examples/user-source/user.h 
