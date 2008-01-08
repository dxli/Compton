#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.61)
AC_INIT(FULL-PACKAGE-NAME, VERSION, BUG-REPORT-ADDRESS)
AC_PROG_CXX
AM_INIT_AUTOMAKE
#AUTOMAKE_OPTIONS = -foreign -Wall
#AM_CXXFLAGS = $(INTI_CFLAGS)
AC_CONFIG_SRCDIR([CPassThrough.cpp])
AC_CONFIG_HEADER([config.h])

# Checks for programs.

# Checks for libraries.

# Checks for header files.
AC_CHECK_HEADERS([sys/time.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_TIME

# Checks for library functions.
#AC_HAVE_LIBRARY(gsl,[],[echo "\tgsl package not found"])
AC_HAVE_LIBRARY(gslcblas,[],[echo "libgslcblas not found"])
AC_CHECK_LIB([gsl],[gsl_integration_qags],[],[echo "couldn't -lgsl"])
AC_CHECK_FUNCS([gettimeofday ])

#AC_HAVE_LIBRARY(mysqlpp,[],[echo "\tmysql++ package not found"])
AC_CONFIG_FILES([Makefile])
AC_OUTPUT
