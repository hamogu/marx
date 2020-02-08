AC_INIT([jdfits],[1.66],[marx-help@space.mit.edu])
# Version is for informational purposes here only.
# See jdfits.h for the version definition.

AC_CONFIG_SRCDIR(src/fitsio.c)
AC_CONFIG_AUX_DIR(build-aux)
AC_CONFIG_MACRO_DIR(build-aux)

AM_INIT_AUTOMAKE
AM_PROG_AR

AC_PROG_MAKE_SET
AC_PROG_RANLIB
AC_PROG_INSTALL

#Check these header since they cause trouble
AC_CHECK_HEADERS([stdlib.h unistd.h limits.h])

AC_CHECK_FUNCS([fseeko memset strchr])
AC_FUNC_MALLOC
AC_FUNC_REALLOC

#Check endianess
AC_CHECK_HEADERS
AC_C_BIGENDIAN

# Check sizeof various types
AC_CHECK_SIZEOF(short)
AC_CHECK_SIZEOF(int)
AC_CHECK_SIZEOF(long)
AC_CHECK_SIZEOF(float)
AC_CHECK_SIZEOF(double)
AC_CHECK_LIB(m, erf)
AC_CHECK_FUNCS( isinf isnan finite erf )


#JD_LARGE_FILE_SUPPORT
AC_SYS_LARGEFILE
AC_FUNC_FSEEKO
AC_TYPE_OFF_T
AC_CHECK_SIZEOF(off_t)

# Define preprocessor FSEEK to be other fseek or fseeko
# HMG: This seems to dublicate code in _jdfits.h
# Thus I'm commenting it out here for testing.
#if test "x$ac_cv_func_fseeko" != xyes; then
#  AC_DEFINE([FSEEK(a,b,c)], [fseeko(a,b,c)],[Def of fseek / fseeko])
#  AC_DEFINE([FTELL(a)], [ftello(a)],[Def of ftell / ftello])
#else
#  AC_DEFINE([FSEEK(a,b,c)], [fseek(a,b,c)],[Def of fseek / fseeko])
#  AC_DEFINE([FTELL(a)], [ftell(a)],[Def of ftell / ftello])
#fi

AC_CONFIG_HEADER(src/config.h:src/config.hin)
AC_CONFIG_FILES([Makefile])
AC_OUTPUT

