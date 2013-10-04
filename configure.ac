# -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT([carna], [1.2.3], [])

AM_INIT_AUTOMAKE

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC

AC_LANG_CPLUSPLUS

# may be useful later
# AC_PATH_PROGS([PERL],[perl perl5.004 perl5.003 perl5.002 perl5.001 perl5])
# AC_PATH_PROG([HELP2MAN],[help2man],false)
# AC_PATH_PROG([POD2MAN],[pod2man],false)

# ----------------------------------------
# Vienna RNA package library path support, if not installed in standard directory/path
#AC_ARG_WITH([vrna],
#    AC_HELP_STRING(
#        [--with-vrna=PREFIX],
#        [Installation directory of the Vienna RNA library]
#      )
#) 

#AC_ARG_VAR([RNAfold],[RNAfold name (def=RNAfold)])
#AC_ARG_VAR([RNAplfold],[RNAplfold name (def=RNAplfold)])
#AC_ARG_VAR([RNAalifold],[RNAalifold name (def=RNAalifold)])


#if test "$with_vrna" = "" ; then
#   AC_PATH_PROG([RNAfold],[RNAfold],[false])
#   AC_PATH_PROG([RNAplfold],[RNAplfold],[false])
#   AC_PATH_PROG([RNAalifold],[RNAalifold],[false])
#else
#   AC_PATH_PROG([RNAfold],[RNAfold],[false],[$with_vrna/bin])
#   AC_PATH_PROG([RNAplfold],[RNAplfold],[false],[$with_vrna/bin])
#   AC_PATH_PROG([RNAalifold],[RNAalifold],[false],[$with_vrna/bin])  
#fi


# error output if ViennaRNA not found
#VRNA_OK=true;
#if test "$RNAfold" = "false" -o "$RNAplfold" = "false" -o "$RNAalifold" = "false" ; then
#    VRNA_OK=false
#fi

# end checking for Vienna package and program paths
# ----------------------------------------


# Checks for header files.
AC_HEADER_STDC


##################
# Debug option
#
AC_MSG_CHECKING([whether to build with debug information])
debuger=no
AC_ARG_ENABLE([debug],
    AC_HELP_STRING(
        [--enable-debug],
        [enable debug data generation (def=no)]
    ),
    debuger="$enableval"
)
AC_MSG_RESULT($debuger)
if test x"$debuger" = x"yes"; then
    AC_DEFINE([DEBUG], [], [Define to 1 for turning on debugging])
    CPPFLAGS="$CPPFLAGS -g -Wall"
else
    AC_DEFINE([NDEBUG], [], [Define to 1 for turning off debugging])
fi

####################
## enable/disable graphical search with gist
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

# turn off deprecated warnings 
CPPFLAGS="$CPPFLAGS -Wno-deprecated"

# configure for use of LocARNA library
PKG_CHECK_MODULES([DEPS], [LocARNA-1.2 >= 1.7.10])

## where are GECODE libs and includes
##

AC_MSG_CHECKING([Where are GECODE libs and headers])
if test "$GECODE_HOME" != "" ; then
   LDFLAGS="$LDFLAGS -L$GECODE_HOME/lib"
   CPPFLAGS="$CPPFLAGS -I$GECODE_HOME/include"
fi
AC_ARG_VAR(GECODE_HOME,[Home of GECODE installation (libs and headers)])
AC_MSG_RESULT([$GECODE_HOME])


GECODE_LIBFLAGS="-lgecodesearch -lgecodeminimodel -lgecodeint -lgecodeset -lgecodesupport -lgecodekernel"

if test $ENABLE_GIST = "yes" ; then
   GECODE_LIBFLAGS="-lgecodegist -lX11 $GECODE_LIBFLAGS"
fi

 # add GECODE libraries to list of libs for linking
LIBS="$LIBS $GECODE_LIBFLAGS"


#####################
## Check for CImg library
AC_CHECK_HEADER(CImg.h,
	[CIMG_AVAILABLE="yes"],
        [CIMG_AVAILABLE="no"],
	[])

AC_MSG_CHECKING([whether CImg has correct version])
AC_COMPILE_IFELSE([AC_LANG_PROGRAM(
  [[#include <CImg.h>
  #ifndef cimg_version 
  # error cimg not defined
  #endif
  #if cimg_version < 132
  # error wrong version
  #endif
  ]])],
  [CIMG_CORRECT_VERSION="yes"],
  [CIMG_CORRECT_VERSION="no"])
AC_MSG_RESULT([$CIMG_CORRECT_VERSION])


####################
# Static linking
#
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


# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_TYPE_SIZE_T


AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile])

#AC_FUNC_MALLOC


## ----------------------------------------
## Give information in case Vienna package missing
##

# if test "$VRNA_OK" = "false" ; then
#   AC_MSG_NOTICE([============================================================]#)
#   AC_MSG_NOTICE()
#     if test "$with_vrna" != "" ; then
#        AC_MSG_NOTICE(Cannot find Vienna RNA package in given path '$with_vrna'.)
#     else
#        AC_MSG_NOTICE(Cannot find Vienna RNA package in standard path, for non-standard path use --with-vrna=PREFIX.)
#     fi
#    AC_MSG_NOTICE()
#    AC_MSG_NOTICE(The Vienna RNA package is required for full functionality (use version >= 1.8.*).)
#    AC_MSG_NOTICE(It can be obtained from http://www.tbi.univie.ac.at/~ivo/RNA/.)
#    AC_MSG_NOTICE()
# fi


## ----------------------------------------
## Complain if gist enabled but CImg.h missing or wrong version

if test "$ENABLE_GIST" = "yes" ; then
   if test "$CIMG_AVAILABE" = "no" ; then
     AC_MSG_NOTICE([])
     AC_MSG_NOTICE([----------------------------------------])
     AC_MSG_NOTICE([ERROR: "CImg could not be found. Please install CImg library."])
     exit -1
   else
      if test $CIMG_CORRECT_VERSION = "no" ; then
         AC_MSG_NOTICE([])
     	 AC_MSG_NOTICE([----------------------------------------])
     	 AC_MSG_NOTICE([ERROR: "Wrong CImg version please use >=1.3.2"])
	 exit -1
      fi
   fi
   
   # define HAVE_GIST, such that we can use this flag for conditional
   # compilation
   AC_DEFINE([HAVE_GIST],[],[Whether gist is available])
   AC_MSG_NOTICE([])
   AC_MSG_NOTICE([----------------------------------------])
   AC_MSG_NOTICE([Gist support enabled])
else
   AC_MSG_NOTICE([])
   AC_MSG_NOTICE([----------------------------------------])
   AC_MSG_NOTICE([Gist support disabled])
fi


AC_OUTPUT([src/Makefile Utils/Makefile])

