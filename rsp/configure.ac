# Process this file with autoconf to produce a configure script.
AC_INIT([marxrsp],[6.0],[marx-help@space.mit.edu],[],[http://space.mit.edu/cxc/marx/])

# Maybe it's easier to put things in directories named
# for GNU standards?
AC_CONFIG_SRCDIR(src/marxrsp.c)

# Do I need this? Not sure there is anything left to keep in there.
AC_CONFIG_AUX_DIR([build-aux])


AM_INIT_AUTOMAKE([-Wall -Werror])
AM_PROG_AR
AC_PROG_CC

# Do I need this? Not suere there are any macros left I need to build.
AC_CONFIG_MACRO_DIR([build-aux])

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


AC_CONFIG_HEADERS(src/config.h:src/config.hin)

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
