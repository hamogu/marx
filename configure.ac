AC_INIT([marx], 5.2, [marx-help@space.mit.edu])
AM_INIT_AUTOMAKE

AC_CONFIG_SUBDIRS(jdmath pfile jdfits)

AC_CONFIG_SRCDIR(libsrc/diffract.c)

AC_PROG_MAKE_SET
AC_PROG_RANLIB
AC_PROG_INSTALL
AC_CANONICAL_HOST

AC_PATH_XTRA

AC_CHECK_LIB(cfitsio, ffgky, [], [AC_MSG_ERROR([cfitsio library not found. Specify path to cfitsio with "./configure LDFLAGS=-L/path/to/cfitsio.a"])])
AC_CHECK_HEADER([fitsio.h],[],[AC_MSG_ERROR([cfitsio header file not found. Specify path with ".configure CFLAGS=-I/path/to/fitsio.h"])])

#Check these header since they cause trouble
AC_CHECK_HEADERS( \
stdlib.h \
unistd.h \
)


AC_CHECK_FUNCS( timegm )

dnl JD used to check header and lib with an AND condition.
AC_CHECK_HEADERS(dlfcn.h)
AC_CHECK_LIB(dl, dlopen)
AS_IF([test "x$ac_cv_header_dlfcn_h" = "xyes" -a "x$ac_cv_lib_dl_dlopen" = "xyes"],
  [AC_DEFINE([MARX_HAS_DYNAMIC_LINKING], [], [dynamic linkin available for user source])],
 [AC_MSG_WARN("No dynamic linking. Source type USER unavailable.")])

dnl  Check sizeof various types
AC_CHECK_SIZEOF(short)
AC_CHECK_SIZEOF(int)
AC_CHECK_SIZEOF(long)
AC_CHECK_SIZEOF(float)
AC_CHECK_SIZEOF(double)

#JD_LARGE_FILE_SUPPORT
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


AC_CONFIG_HEADER(libsrc/config.h:libsrc/config.hin)


AC_OUTPUT(Makefile)


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
