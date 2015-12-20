AC_INIT([pfile], [2.4.1], [marx-help@space.mit.edu])
# Version is for informational porpuses here only.
# See pfile.h for the version definition.

AC_CONFIG_SRCDIR(src/pfile.c)
AC_CONFIG_AUX_DIR(build-aux)
AC_CONFIG_MACRO_DIR(build-aux)

AM_INIT_AUTOMAKE

AC_PROG_MAKE_SET
AC_PROG_RANLIB
AC_PROG_INSTALL

#Check these header since they cause trouble
AC_CHECK_HEADERS( \
stdlib.h \
unistd.h \
limits.h \
)

#JD_LARGE_FILE_SUPPORT
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


AC_CONFIG_HEADER(src/config.h:src/config.hin)
AC_OUTPUT(Makefile)

