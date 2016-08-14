dnl Check for NUMA
dnl Charles Bouillaguet 27-03-2015
dnl Boyer Brice 07/04/11
dnl Bradford Hovinen, 2001-06-13
dnl Modified by Pascal Giorgi, 2003-12-03
dnl Inspired by gnome-bonobo-check.m4 by Miguel de Icaza, 99-04-12
dnl Stolen from Chris Lahey       99-2-5
dnl stolen from Manish Singh again
dnl stolen back from Frank Belew
dnl stolen from Manish Singh
dnl Shamelessly stolen from Owen Taylor

AC_DEFUN([SPASM_CHECK_NUMA],
[

AC_ARG_WITH(numa,
[AC_HELP_STRING([--with-numa=<path>|yes], [Use the NUMA library. If argument is yes or <empty>,
    that means the library is reachable with the standard
    search path (/usr or /usr/local). Otherwise you give
    the <path> to the directory which contains the
    library.
])],
    [if test "$withval" = yes ; then
        NUMA_HOME_PATH="${DEFAULT_CHECKING_PATH}"
        elif test "$withval" != no ; then
        NUMA_HOME_PATH="$withval ${DEFAULT_CHECKING_PATH}"
        fi],
    [NUMA_HOME_PATH="${DEFAULT_CHECKING_PATH}"])


dnl Check for existence
BACKUP_CFLAGS=${CFLAGS}
BACKUP_LIBS=${LIBS}

AC_MSG_CHECKING(for NUMA)

for NUMA_HOME in ${NUMA_HOME_PATH}
  do
#    AC_MSG_NOTICE($NUMA_HOME)
    if test -r "$NUMA_HOME/include/numa.h"; then

       if test "x$NUMA_HOME" != "x/usr"; then
           NUMA_CFLAGS="-I${NUMA_HOME}/include"
           NUMA_LIBS="-L${NUMA_HOME}/lib -lnuma"
       else
           NUMA_CFLAGS=
           NUMA_LIBS="-lnuma"
       fi

#       AC_MSG_NOTICE(found include)


       CFLAGS="${BACKUP_CFLAGS} ${NUMA_CFLAGS}"
       LIBS="${BACKUP_LIBS} ${NUMA_LIBS}"

       AC_LINK_IFELSE(
       [AC_LANG_PROGRAM([[#include <numa.h>
                          #include <stdlib.h>]],
        [[ numa_available(); ]]
	)],
	[
	  numa_found="yes"
	  break
	],
	[
	   numa_found="no (numa.h found but linking failed)"
	   unset NUMA_CFLAGS
	   unset NUMA_LIBS
	])
    else
       numa_found="no (numa.h not found)"
    fi
done

AC_MSG_RESULT($numa_found)

if test "x$numa_found" = "xyes" ; then
    AC_SUBST(NUMA_CFLAGS)
    AC_SUBST(NUMA_LIBS)
    AC_DEFINE(HAVE_NUMA,1,[Define if NUMA is installed])
    HAVE_NUMA=yes
fi

AM_CONDITIONAL(SPASM_HAVE_NUMA, test "x$HAVE_NUMA" = "xyes")

CFLAGS=${BACKUP_CFLAGS}
LIBS=${BACKUP_LIBS}

])
