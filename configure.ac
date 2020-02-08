# Process this file with autoconf to produce a configure script.
AC_INIT([marx],[6.0],[marx-help@space.mit.edu],[],[http://space.mit.edu/cxc/marx/])

# Maybe it's easier to put things in directories named
# for GNU standards?
AC_CONFIG_SRCDIR(libsrc/diffract.c)

# Do I need this? Not sure there is anything left to keep in there.
AC_CONFIG_AUX_DIR([build-aux])


AM_INIT_AUTOMAKE([-Wall -Werror])
AM_PROG_AR
AC_PROG_CC

# Do I need this? Not suere there are any macros left I need to build.
AC_CONFIG_MACRO_DIR([build-aux])

AC_CONFIG_SUBDIRS(jdmath pfile jdfits rsp)

AC_PROG_MAKE_SET
AC_PROG_RANLIB
AC_PROG_INSTALL
AC_CANONICAL_HOST

AC_PATH_XTRA

# Make cfitsio a mandatory dependency for now
# In principle, that's optional, but lets get this build system working
# before I add options.
# Use dylib instead of static .a files?
AC_CHECK_LIB(cfitsio, ffopentest, [],
              [AC_MSG_ERROR([cfitsio library not found. Specify path to cfitsio with "./configure LDFLAGS=-L/path/to/libcfitsio.a"])])

# Check autoheader works: https://www.gnu.org/savannah-checkouts/gnu/autoconf/manual/autoconf-2.69/autoconf.html#autoheader-Invocation
AC_CHECK_HEADER([fitsio.h],
                [AC_DEFINE(HAVE_FITSIO_H, 1,
                        [Define to 1 if you have fitsio.h .])],
		[AC_MSG_ERROR([cfitsio header file not found. Specify path with 	".configure CFLAGS=-I/path/include/fitsio.h"])])

#Check these header since they cause trouble
AC_CHECK_HEADERS([stdlib.h unistd.h limits.h])


AC_CHECK_FUNCS([timegm floor memset mkdir pow sqrt strchr strstr])
AC_FUNC_MALLOC
AC_FUNC_MKTIME
AC_FUNC_REALLOC
AC_FUNC_STRTOD
AC_FUNC_ERROR_AT_LINE

dnl JD used to check header and lib with an AND condition.
AC_CHECK_HEADERS(dlfcn.h)
AC_CHECK_LIB(dl, dlopen)
AS_IF([test "x$ac_cv_header_dlfcn_h" = "xyes" -a "x$ac_cv_lib_dl_dlopen" = "xyes"],
  [AC_DEFINE([MARX_HAS_DYNAMIC_LINKING], [1], [Define to 1 if dynamic linking is available for USER source])],
 [AC_MSG_WARN("No dynamic linking. Source type USER unavailable.");AC_DEFINE([MARX_HAS_DYNAMIC_LINKING], [0], [Define to 1 if dynamic linking is available for USER source])])

AM_CONDITIONAL([HAVE_LIBDL], [test "x$ac_cv_header_dlfcn_h" = "xyes" -a "x$ac_cv_lib_dl_dlopen" = "xyes"]) 
 
dnl  Check sizeof various types
AC_CHECK_SIZEOF(short)
AC_CHECK_SIZEOF(int)
AC_CHECK_SIZEOF(long)
AC_CHECK_SIZEOF(float)
AC_CHECK_SIZEOF(double)
# Need to use SIZEOF_INT_P etc in the program for this to make sense
# I think original intention was to raise error if length is different
# then expected. 


# TBD: Is this needed here? Or enough in jdfits and pfile?
AC_SYS_LARGEFILE
AC_FUNC_FSEEKO
AC_TYPE_OFF_T
AC_CHECK_SIZEOF(off_t)

# Define preprocessor FSEEK to be either fseek or fseeko
if test "x$ac_cv_func_fseeko" != xyes; then
  AC_DEFINE([FSEEK(a,b,c)], [fseeko(a,b,c)],[Def of fseek / fseeko])
  AC_DEFINE([FTELL(a)], [ftello(a)],[Def of ftell / ftello])
else
  AC_DEFINE([FSEEK(a,b,c)], [fseek(a,b,c)],[Def of fseek / fseeko])
  AC_DEFINE([FTELL(a)], [ftell(a)],[Def of ftell / ftello])
fi


AC_CONFIG_HEADERS(libsrc/config.h:libsrc/config.hin)

AC_CONFIG_FILES([Makefile par/marx.par par/marxasp.par par/marxpileup.par])
AC_OUTPUT


echo ""
echo "You are compiling MARX with the following compiler configuration:"
echo "       CC =" "$CC"
echo "   CFLAGS =" "$CFLAGS"
echo "  LDFLAGS =" "$LDFLAGS" "$DYNAMIC_LINK_FLAGS"
echo ""
echo "MARX executables will be installed in $prefix/bin/."
echo "MARX documentation will be installed in $prefix/share/doc/marx/."
echo "MARX data files will be installed in $prefix/share/marx/data/."
echo "MARX parameter files will be installed in $prefix/share/marx/pfiles/."
echo ""
echo "To continue the build process, run 'make'."
echo ""


### Notes
# Would like to use ASCDS_LIB and ASCDS_LIB/../include as default
# LDFLAGS and CFLAGS unless given by the user.
# Here is one way (from https://lists.gnu.org/archive/html/autoconf/2006-04/msg00022.html) to do that:
# I first check to see if CFLAGS is set, then
# only muck with it if it has not been set:
#
# AC_MSG_CHECKING([whether configure should try to set CFLAGS])
# if test "x${CFLAGS+set}" = xset; then
#    enable_cflags_setting=no
# else
#    enable_cflags_setting=yes
# fi
# AC_MSG_RESULT($enable_cflags_setting)
#
# then later...
#
#        test "x$enable_cflags_setting" = xyes && CFLAGS="$CFLAGS -q64"
#
# This way, if the user sets CFLAGS, his choice is respected (though it
# might make the build crash - but that's his problem).
