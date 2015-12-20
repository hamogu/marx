AC_INIT([jdmath], [1.83], [marx-help@space.mit.edu])
# Version is for informational purposes here only.
# See jdmath.h for the version definition.

AC_CONFIG_SRCDIR(src/vector.c)
AC_CONFIG_AUX_DIR(build-aux)
AC_CONFIG_MACRO_DIR(build-aux)

AM_INIT_AUTOMAKE

AC_PROG_MAKE_SET
AC_PROG_RANLIB
AC_PROG_INSTALL
AC_CANONICAL_HOST

AC_CHECK_HEADERS
AC_C_BIGENDIAN

#Check these header since they cause trouble
AC_CHECK_HEADERS( \
stdlib.h \
unistd.h \
)

# Check sizeof various types
AC_CHECK_SIZEOF(short, 2)
AC_CHECK_SIZEOF(int, 4)
AC_CHECK_SIZEOF(long, 4)
AC_CHECK_SIZEOF(float, 4)
AC_CHECK_SIZEOF(double, 8)
AC_CHECK_LIB(m, erf)
AC_CHECK_FUNCS( isinf isnan finite erf )

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


AC_CONFIG_HEADER(src/config.h)
AC_OUTPUT(Makefile demo/Makefile)
