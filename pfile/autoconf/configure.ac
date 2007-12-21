AC_INIT(src/pfile.c)
JD_INIT
JD_ANSI_CC

AC_PREFIX_DEFAULT(..)

AC_CONFIG_AUX_DIR(autoconf)
AC_PROG_MAKE_SET
AC_PROG_RANLIB
AC_PROG_INSTALL

#Check these header since they cause trouble
AC_CHECK_HEADERS( \
stdlib.h \
unistd.h \
limits.h \
)

JD_LARGE_FILE_SUPPORT

JD_SET_OBJ_SRC_DIR
JD_GCC_WARNINGS
JD_GET_MODULES(src/modules.lis)
JD_CREATE_MODULE_ORULES

AC_CONFIG_HEADER(src/sysconf.h:src/config.hin)
AC_OUTPUT(Makefile:autoconf/Makefile.in src/Makefile)

JD_CREATE_EXEC_RULE(pset)
JD_CREATE_EXEC_RULE(plist)
JD_CREATE_EXEC_RULE(pwhich)

JD_APPEND_RULES(src/Makefile)
JD_APPEND_ELFRULES(src/Makefile)
