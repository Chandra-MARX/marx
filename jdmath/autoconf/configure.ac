AC_INIT
AC_CONFIG_SRCDIR([src/vector.c])
AC_PREFIX_DEFAULT(/usr/local)
JD_INIT

AC_CONFIG_AUX_DIR(autoconf)
AC_PROG_MAKE_SET
AC_PROG_RANLIB
AC_PROG_INSTALL
AC_CANONICAL_HOST

#JD_ANSI_CC
#ELF_CFLAGS="$ELF_CFLAGS"


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

#JD_LARGE_FILE_SUPPORT

JD_SET_OBJ_SRC_DIR
JD_GCC_WARNINGS
JD_GET_MODULES(src/modules.lis)
JD_CREATE_MODULE_ORULES

AC_CONFIG_HEADERS([src/sysconf.h:src/config.hin])
AC_CONFIG_FILES([Makefile:autoconf/Makefile.in src/Makefile demo/Makefile])
AC_OUTPUT

JD_APPEND_RULES(src/Makefile)
JD_APPEND_ELFRULES(src/Makefile)
