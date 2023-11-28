AC_INIT
AC_CONFIG_SRCDIR([libsrc/diffract.c])
JD_INIT

JD_ANSI_CC

AC_CONFIG_AUX_DIR(autoconf)
AC_PROG_MAKE_SET
AC_PROG_RANLIB
AC_PROG_INSTALL

AC_CANONICAL_HOST

dnl Check these header since they cause trouble
AC_CHECK_HEADERS( \
stdlib.h \
unistd.h \
)

AC_CHECK_FUNCS( timegm )

AC_CHECK_LIB(dl,dlopen,[
	AC_CHECK_HEADER(dlfcn.h, [
		DYNAMIC_LINK_LIB="-ldl"
		AC_DEFINE(HAVE_DLFCN_H)
	       ],[
	        DYAMIC_LINK_LIB=""
		AC_MSG_WARN(cannot perform dynamic linking)
               ])],[
	DYAMIC_LINK_LIB=""
	AC_MSG_WARN(cannot perform dynamic linking)])
AC_SUBST(DYNAMIC_LINK_LIB)

case "$host_os" in
  linux )
    DYNAMIC_LINK_FLAGS="-Wl,-export-dynamic"
    ;;
  * )
    DYNAMIC_LINK_FLAGS=""
    ;;
esac
AC_SUBST(DYNAMIC_LINK_FLAGS)


dnl  Check sizeof various types
AC_CHECK_SIZEOF(short, 2)
AC_CHECK_SIZEOF(int, 4)
AC_CHECK_SIZEOF(long, 4)
AC_CHECK_SIZEOF(float, 4)
AC_CHECK_SIZEOF(double, 8)

JD_LARGE_FILE_SUPPORT

JD_SET_OBJ_SRC_DIR(libsrc)

JD_SIMPLE_LIB_DIR(jdmath)
JD_SIMPLE_LIB_DIR(jdfits)
JD_SIMPLE_LIB_DIR(pfile)

JD_GCC_WARNINGS
JD_GET_MODULES(libsrc/modules.lis)
JD_CREATE_MODULE_ORULES


AC_CONFIG_HEADERS([libsrc/sysconf.h:libsrc/config.hin])

CF_Output_Files="Makefile:autoconf/Makefile.in libsrc/Makefile"
if test -d src
then
  CF_Output_Files="$CF_Output_Files src/Makefile"
fi

AC_CONFIG_FILES([$CF_Output_Files data/Makefile doc/Makefile])
AC_OUTPUT

JD_APPEND_RULES(libsrc/Makefile)

dnl We may be configuring just marxlib and the src directory may not be
dnl part of the distribution

if test -d src
then

dnl
dnl  Now do src subdir
dnl  Note that AC_OUTPUT can only be called once.  For that reason,
dnl  src/Makefile was created above
dnl
echo "Creating Makefile in src"

JD_GET_MODULES()
JD_CREATE_EXEC_RULE(marx)
JD_CREATE_EXEC_RULE(bestfocus)
JD_CREATE_EXEC_RULE(rfl)
JD_CREATE_EXEC_RULE(geff)
JD_CREATE_EXEC_RULE(geffdump)
JD_CREATE_EXEC_RULE(dumpwfold)
JD_CREATE_EXEC_RULE(marxcat)
JD_CREATE_EXEC_RULE(marx2fits)
JD_CREATE_EXEC_RULE(marxpileup)
JD_CREATE_EXEC_RULE(acis-qe)
JD_CREATE_EXEC_RULE(testacis)
JD_CREATE_EXEC_RULE(detinfo)
JD_CREATE_EXEC_RULE(sky2chip)
JD_CREATE_EXEC_RULE(marxasp)

JD_CREATE_ORULE(acissim)
JD_CREATE_ORULE(henke)

JD_APPEND_RULES(src/Makefile)

cd src; /bin/rm -f marx-feat.h; ln -s ../libsrc/marx-feat.h marx-feat.h

fi
