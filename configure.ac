# -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT([ConstraintRNAAlignment], [0.1], [])

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
    AC_DEFINE([DEBUG])
    CPPFLAGS="$CPPFLAGS -g -Wall"
else
    AC_DEFINE([NDEBUG])
fi

CPPFLAGS="$CPPFLAGS -Wno-deprecated"

## where are GECODE libs and includes
##

AC_MSG_CHECKING([Where are GECODE libs])
if test "$GECODE_HOME" != "" ; then
   LDFLAGS="$LDFLAGS -L$GECODE_HOME/lib"
fi
AC_ARG_VAR(GECODE_HOME,[Home of GECODE installation (includes libs)])
AC_MSG_RESULT([$GECODE_HOME])

AC_MSG_CHECKING([Where are GECODE headers])
if test "$GECODE_SRC_HOME" != "" ; then
  CPPFLAGS="$CPPFLAGS -I$GECODE_SRC_HOME"
fi
AC_ARG_VAR(GECODE_SRC_HOME,[Home of GECODE source (contains headers)])
AC_MSG_RESULT([$GECODE_SRC_HOME])


GECODE_LIBFLAGS="-lgecodekernel -lgecodesearch -lgecodeint -lgecodeset -lgecodeminimodel  -lgecodesupport  -lgecodegist"

LDFLAGS="$LDFLAGS $GECODE_LIBFLAGS"


#####################
## Check for CImg library

AC_MSG_CHECKING([whether CImg is installed])
AC_CHECK_HEADER(CImg.h,
	[],
        [AC_MSG_ERROR([Please install CImg library])],	
	[])

AC_COMPILE_IFELSE([AC_LANG_PROGRAM(
     [[#include <CImg.h>
     #ifndef cimg_version 
     # error cimg not defined
     #endif
     #if cimg_version != 132
     # error wrong version
     #endif
     ]])],
                        [],
                        AC_MSG_ERROR([Wrong CImg version please use 132]))



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

AC_CONFIG_FILES([Makefile])

#AC_FUNC_MALLOC


AC_OUTPUT([src/Makefile])


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
