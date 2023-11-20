dnl -*- sh -*-
AC_INIT
AC_CONFIG_SRCDIR([marx/par/marx.par])
JD_INIT

JD_ANSI_CC

AC_CONFIG_AUX_DIR(autoconf)
AC_PROG_MAKE_SET
AC_PROG_RANLIB
AC_PROG_INSTALL
AC_CANONICAL_HOST

AC_PATH_XTRA
JD_CHECK_FOR_LIBRARY(cfitsio, fitsio.h)

#Check these header since they cause trouble
AC_CHECK_HEADERS( \
stdlib.h \
unistd.h \
)

AC_CONFIG_SUBDIRS(jdmath pfile jdfits marx)
AC_CONFIG_FILES([Makefile:autoconf/Makefile.in marxrsp/Makefile])
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
