AC_INIT([jdmath], 1.83, [marx-help@space.mit.edu])
AC_PREFIX_DEFAULT(/usr/local)
AM_INIT_AUTOMAKE

AC_CONFIG_SRCDIR(src/vector.c)
AC_PROG_MAKE_SET
AC_PROG_RANLIB
AC_PROG_INSTALL
AC_CANONICAL_HOST

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
AC_CHECK_FUNCS( isinf isnan finite )

AC_CHECK_LIB(m, erf, [AC_DEFINE(HAVE_ERF)])

AC_CONFIG_HEADER(src/sysconf.h:src/config.hin)
AC_OUTPUT(Makefile demo/Makefile)
