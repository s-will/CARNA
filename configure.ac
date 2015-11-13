dnl -*- Autoconf -*-
dnl Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT([carna], [1.3.0], [])

AM_INIT_AUTOMAKE

dnl Checks for programs.
AC_PROG_CXX
AC_PROG_CC

AC_LANG_CPLUSPLUS

dnl Checks for header files.
AC_HEADER_STDC

dnl Set directory for m4 macro files 
AC_CONFIG_MACRO_DIR([m4])

dnl --------------------
dnl Debug option

AC_MSG_CHECKING([whether to build with debug information])
debugger=no
AC_ARG_ENABLE([debug],
    AC_HELP_STRING(
        [--enable-debug],
        [enable debug data generation (def=no)]
    ),
    debugger="$enableval"
)
AC_MSG_RESULT($debugger)
if test x"$debugger" = x"yes"; then
    AC_DEFINE([DEBUG], [], [Define to 1 for turning on debugging])
    CPPFLAGS="$CPPFLAGS -g -Wall"
else
    AC_DEFINE([NDEBUG], [], [Define to 1 for turning off debugging])
fi

dnl --------------------
dnl enable/disable graphical search with gist
AC_MSG_CHECKING([whether to build with gist support])
ENABLE_GIST=yes
AC_ARG_ENABLE([gist],
    AC_HELP_STRING(
        [--disable-gist],
        [disable gist support (def=enabled)]
    ),
    [ENABLE_GIST="$enableval"]
)
AC_MSG_RESULT([$ENABLE_GIST])

dnl turn off deprecated warnings 
CPPFLAGS="$CPPFLAGS -Wno-deprecated -Wno-unused-local-typedefs"

dnl configure for use of LocARNA library
PKG_CHECK_MODULES([DEPS], [LocARNA-1.4 >= 1.8.6])

dnl where are GECODE libs and includes

AC_MSG_CHECKING([Where are GECODE libs and headers])
if test "$GECODE_HOME" != "" ; then
   LDFLAGS="$LDFLAGS -L$GECODE_HOME/lib"
   CPPFLAGS="$CPPFLAGS -I$GECODE_HOME/include"
fi
AC_ARG_VAR(GECODE_HOME,[Home of GECODE installation (libs and headers)])
AC_MSG_RESULT([$GECODE_HOME])


GECODE_LIBFLAGS="-lgecodesearch -lgecodeminimodel -lgecodeint -lgecodeset -lgecodesupport -lgecodekernel"

dnl add GECODE libraries to list of libs for linking
LIBS="$LIBS $GECODE_LIBFLAGS"

if test "$ENABLE_GIST" = "yes" ; then
   AC_CHECK_LIB([gecodegist],[main],[HAVE_LIBGECODEGIST="yes"])
   AC_CHECK_LIB([X11],[main],[HAVE_LIBX11="yes"])
   LIBS="$LIBS -lgecodegist -lX11"
fi

if test "$ENABLE_GIST" = "yes" ; then
dnl --------------------
dnl Check for CImg library
AC_CHECK_HEADER(CImg.h,
	[CIMG_AVAILABLE="yes"],
        [CIMG_AVAILABLE="no"],
	[])

AC_MSG_CHECKING([whether CImg has correct version])
AC_COMPILE_IFELSE([AC_LANG_PROGRAM(
  [[#include <CImg.h>
  #ifndef cimg_version 
  #  error cimg not defined
  #endif
  #if cimg_version < 132
  #  error wrong version
  #endif
  ]])],
  [CIMG_CORRECT_VERSION="yes"],
  [CIMG_CORRECT_VERSION="no"])
AC_MSG_RESULT([$CIMG_CORRECT_VERSION])
fi

dnl --------------------
dnl add linker options for POSIX threads
dnl (necessary for linking Gecode) 
AX_PTHREAD
AC_MSG_NOTICE([pthread: $PTHREAD_CFLAGS, $PTHREAD_LIBS])
LIBS="$PTHREAD_LIBS $LIBS"
dnl CFLAGS="$CFLAGS $PTHREAD_CFLAGS"
LDFLAGS="$PTHREAD_CFLAGS $LDFLAGS"

dnl --------------------
dnl Static linking

AC_MSG_CHECKING([whether to link statically])
static_linkage=no
AC_ARG_ENABLE([static],
    AC_HELP_STRING(
        [--enable-static],
        [link statically (def=no)]
    ),
    static_linkage="$enableval"
)
AC_MSG_RESULT($static_linkage)
if test x"$static_linkage" = x"yes"; then
   LDFLAGS="$LDFLAGS -static"
fi


dnl Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_TYPE_SIZE_T


AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile])


dnl ------------------------------------------------------------
dnl If gist is enabled, check whether prerequisites are met.
dnl Otherwise, warn and turn off.

     
if test "$ENABLE_GIST" = "yes" ; then
   AC_MSG_NOTICE([----------------------------------------])
   GIST_OK="yes"

   if test "$CIMG_AVAILABLE" != "yes" ; then
     AC_MSG_NOTICE([CImg could not be found.])
     AC_MSG_NOTICE([Please install CImg (cimg-dev), if you need Gist.])
     GIST_OK="no"
   else
      if test $CIMG_CORRECT_VERSION != "yes" ; then
     	 AC_MSG_NOTICE([Wrong CImg version please use >=1.3.2. Please update, if you need Gist.])
	 GIST_OK="no"
      fi
   fi
   
   if test "$HAVE_LIBGECODEGIST" != "yes" ; then
     AC_MSG_NOTICE([Cannot find library gecodegist.])
     AC_MSG_NOTICE([Likely Gecode was compiled without Gist support.])
     AC_MSG_NOTICE([If you need Gist, please enable it in Gecode."])
     GIST_OK="no"
   fi

   if test "$HAVE_LIBX11" != "yes" ; then
     AC_MSG_NOTICE([Cannot find library X11, which is required for Gist.])
     AC_MSG_NOTICE([Please provide the library (libx11-dev), if you need Gist.])
     GIST_OK="no"
   fi

   dnl if GIST_OK, then define HAVE_GIST, such that we can use this
   dnl flag for conditional compilation
   if test "$GIST_OK" = "yes" ; then
     AC_DEFINE([HAVE_GIST],[],[Whether gist is available])
     AC_MSG_NOTICE([Gist enabled.])
   else
     AC_MSG_NOTICE([Gist disabled. (If you don't know it, you don't need it.)])
   fi
   AC_MSG_NOTICE([----------------------------------------])
fi


AC_OUTPUT([src/Makefile Utils/Makefile])

