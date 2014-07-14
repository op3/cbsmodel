AC_INIT([cbsmodel], [1.1], [reese@ikp.tu-darmstadt.de])
AC_PREREQ([2.64])
AM_INIT_AUTOMAKE([1.11 -Wall])

AC_CONFIG_HEADERS([config.h])
AC_ISC_POSIX
AC_PROG_CXX
AM_PROG_CC_STDC
AC_HEADER_STDC

# check for libm and GSL (GNU scientific library)
# these lines are taken from the GSL manual
AC_CHECK_LIB([m],[cos])
AC_CHECK_LIB([gslcblas],[cblas_dgemm])
AC_CHECK_LIB([gsl],[gsl_blas_dgemm])

# check for readline library
# these lines are taken from the autoconf manual
AC_CHECK_LIB([readline], [main])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
