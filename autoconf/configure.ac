dnl -*- sh -*-
AC_INIT(marx/par/marx.par)

JD_ANSI_CC

AC_CONFIG_AUX_DIR(autoconf)
AC_PROG_MAKE_SET
AC_PROG_RANLIB

#Check these header since they cause trouble
AC_CHECK_HEADERS( \
stdlib.h \
unistd.h \
)

AC_CONFIG_SUBDIRS(jdmath pfile jdfits marx)
AC_OUTPUT(Makefile:autoconf/Makefile.in)

echo ""
echo "Now type 'make'."
echo ""
