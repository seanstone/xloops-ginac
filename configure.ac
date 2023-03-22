dnl Process this file with autoconf to produce a configure script.

AC_INIT
AC_CONFIG_SRCDIR([xloops/xloops.cpp])
AC_PREREQ(2.50)

dnl Configure script options
AC_ARG_WITH(ginac, [  --with-ginac=DIR        Use uninstalled copy of GiNaC])
if [[ "x$with_ginac" = xyes ]]; then
  AC_MSG_ERROR([directory must be specified for --with-ginac])
fi
AC_ARG_ENABLE(ps-doc,   [  --enable-ps-doc         build PostScript documentation [default=no]], , enable_ps_doc=no)

dnl This defines PACKAGE and VERSION
AM_CONFIG_HEADER(config.h)
AM_INIT_AUTOMAKE([xloops-ginac], [0.1.3])

dnl ===========================================================================
dnl Check for the compiler and all the utilities needed for the build:
dnl ===========================================================================
AC_PROG_CXX
AC_PROG_CXXCPP
AC_PROG_INSTALL
AM_PROG_LIBTOOL
AM_PROG_LEX
AC_PROG_YACC

dnl Switch to C++ language mode for the following libraries and headers:
AC_LANG([C++])

dnl Make sure the following libraries work by testing for symbols therein. 
dnl They are automatically added the the variable $LIBS and thus passed into
dnl the Makefile:
AC_CHECK_LIB(stdc++, cout)

dnl Check for stuff needed for building the GiNaC interactive shell (ginsh):
AC_CHECK_HEADERS(unistd.h)
AC_CHECK_HEADERS(readline/readline.h readline/history.h)
if test "x${ac_cv_header_readline_readline_h}" != "xyes" -o "x${ac_cv_header_readline_history_h}" != "xyes"; then
  GINAC_WARNING([I could not find the headers for libreadline (needed for building ginsh).])
fi
GINAC_TERMCAP
save_LIBS=$LIBS
LIBS="$LIBTERMCAP $LIBS"
AC_CHECK_LIB(readline, readline)
if test "x${ac_cv_lib_readline_readline}" = "xyes"; then
  GINAC_LIB_READLINE_VERSION
else
  GINAC_WARNING([I could not find libreadline (needed by ginsh).])
fi
GINSH_LIBS=$LIBS
LIBS=$save_LIBS
AC_SUBST(GINSH_LIBS)

dnl Make sure all the necessary new-style headers are installed on the system.
dnl If one of them cannot be found the system is probably not ANSI-conform
dnl enough so trying the .h-style headers is a waste of time.
AC_CHECK_HEADERS(iostream vector string iterator sstream stdexcept, ,
    AC_MSG_ERROR(need to have ANSI compliant headers))

dnl xloops needs the GiNaC framework up and running:
if [[ "x$with_ginac" != x ]]; then
  CPPFLAGS="$CPPFLAGS -I${with_ginac} -I${with_ginac}/ginac"
  LIBS="$LIBS ${with_ginac}/ginac/libginac.la"
  echo "Creating links for ginsh."
  (cd ginsh; \
    ln -fs ${with_ginac}/ginsh/ginsh.1 .; \
    ln -fs ${with_ginac}/ginsh/ginsh.h .; \
    ln -fs ${with_ginac}/ginsh/ginsh_fcn_help.sed .; \
    ln -fs ${with_ginac}/ginsh/ginsh_op_help.sed .; \
    ln -fs ${with_ginac}/ginsh/ginsh_lexer.ll .; \
    ln -fs ${with_ginac}/ginsh/ginsh_parser.yy . \
  )
else
  AM_PATH_GINAC([0.9.4], , AC_MSG_ERROR([need to have GiNaC installed]))
  CPPFLAGS="$CPPFLAGS $GINACLIB_CPPFLAGS"
  LIBS="$LIBS $GINACLIB_LIBS"
fi

dnl Documentation related stuff
AC_PATH_PROG(LATEX, latex, "")
AC_PATH_PROG(MAKEINDEX, makeindex, "")
AC_PATH_PROG(DVIPS, dvips, "")

DOC_TARGETS=
if test "x${enable_ps_doc}" = "xyes"; then
  DOC_TARGETS="$DOC_TARGETS xloops.ps"
fi
AC_SUBST(DOC_TARGETS)

dnl Output makefiles etc.
AC_CONFIG_FILES([
Makefile
lib/Makefile
xloops/Makefile
check/Makefile
ginsh/Makefile
doc/Makefile
])
AC_OUTPUT

echo "Configuration done. Now type \"make\"."
