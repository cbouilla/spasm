AC_INIT([spasm],[0.1],[charles.bouillaguet@lifl.fr])
AC_PREREQ([2.69])
#AC_CANONICAL_HOST
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_MACRO_DIR([m4])
LT_INIT

AM_INIT_AUTOMAKE([1.8 gnu no-dependencies foreign -Wall])
#AC_CONFIG_HEADERS([config.h])
AM_MAINTAINER_MODE

# user-definable options

# Checks for programs.
AC_PROG_CC
AC_PROG_CC_C99
AC_PROG_MAKE_SET
AC_PROG_AWK     # for TAP
AC_PROG_LIBTOOL
AC_PROG_INSTALL

AC_REQUIRE_AUX_FILE([tap-driver.sh])

# Checks for header files.
AC_HEADER_STDC
AC_HEADER_STDBOOL

# Checks for typedefs, structures, and compiler characteristics.
AC_FUNC_MALLOC

CFLAGS="$CFLAGS -Wall -Wextra -g -fPIC"

AC_SUBST(CFLAGS)
AC_SUBST(LDFLAGS)

AC_CONFIG_FILES([Makefile src/Makefile test/Makefile bench/Makefile])
AC_OUTPUT
