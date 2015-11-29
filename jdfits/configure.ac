AC_INIT([jdfits], 1.66 [marx-help@space.mit.edu])
AM_INIT_AUTOMAKE

AC_CONFIG_SRCDIR(src/fitsio.c)
AC_PROG_MAKE_SET
AC_PROG_RANLIB
AC_PROG_INSTALL

#Check these header since they cause trouble
AC_CHECK_HEADERS( \
stdlib.h \
unistd.h \
)

AC_CHECK_FUNCS(\
fseeko \
)


#Check endianess
AC_CHECK_HEADERS
AC_C_BIGENDIAN

# Check sizeof various types
AC_CHECK_SIZEOF(short, 2)
AC_CHECK_SIZEOF(int, 4)
AC_CHECK_SIZEOF(long, 4)
AC_CHECK_SIZEOF(float, 4)
AC_CHECK_SIZEOF(double, 8)
AC_CHECK_LIB(m, erf)
AC_CHECK_FUNCS( isinf isnan finite erf )


#JD_LARGE_FILE_SUPPORT
AC_SYS_LARGEFILE
AC_FUNC_FSEEKO
AC_TYPE_OFF_T
AC_CHECK_SIZEOF(off_t)

AC_CONFIG_HEADER(src/config.h:src/config.hin)
AC_OUTPUT(Makefile)

