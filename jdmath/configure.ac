#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])

AC_INIT([jdmath],[1.83],[marx-help@space.mit.edu])
# Version is for informational purposes here only.
# See jdmath.h for the version definition.

AC_CONFIG_SRCDIR(src/vector.c)
AC_CONFIG_AUX_DIR(build-aux)
AC_CONFIG_MACRO_DIR(build-aux)

AM_INIT_AUTOMAKE

# Checks for programs.
AC_PROG_CC
AM_PROG_AR

# Which of the ones below are needed?
AC_PROG_MAKE_SET
AC_PROG_RANLIB
AC_PROG_INSTALL
#AC_CANONICAL_HOST

AC_C_BIGENDIAN

# Checks for header files.
AC_CHECK_HEADERS([limits.h stdlib.h string.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_SIZE_T
# Line above was suggested by autoscan.
# But I'm not sure the resulting symbol is currently used in the code.

# Check sizeof various types
AC_CHECK_SIZEOF(short)
AC_CHECK_SIZEOF(int)
AC_CHECK_SIZEOF(long)
AC_CHECK_SIZEOF(float)
AC_CHECK_SIZEOF(double)

AC_SYS_LARGEFILE
AC_FUNC_FSEEKO
AC_TYPE_OFF_T
AC_CHECK_SIZEOF(off_t)

# Define preprocessor FSEEK to be other fseek or fseeko
if test "x$ac_cv_func_fseeko" != xyes; then
  AC_DEFINE([FSEEK(a,b,c)], [fseeko(a,b,c)],[Def of fseek / fseeko])
  AC_DEFINE([FTELL(a)], [ftello(a)],[Def of ftell / ftello])
else
  AC_DEFINE([FSEEK(a,b,c)], [fseek(a,b,c)],[Def of fseek / fseeko])
  AC_DEFINE([FTELL(a)], [ftell(a)],[Def of ftell / ftello])
fi

# Checks for library functions.
AC_CHECK_FUNCS([memset pow sqrtisinf isnan finite erf])
AC_CHECK_LIB(m, erf)

AC_CONFIG_HEADER(src/config.h)
AC_CONFIG_FILES([Makefile demo/Makefile])
AC_OUTPUT
