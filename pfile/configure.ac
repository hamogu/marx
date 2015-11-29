AC_INIT([pfile], 2.4.1, [marx-help@space.mit.edu])
AM_INIT_AUTOMAKE

AC_CONFIG_SRCDIR(src/pfile.c)
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


AC_CONFIG_HEADER(src/config.h:src/config.hin)
AC_OUTPUT(Makefile)

