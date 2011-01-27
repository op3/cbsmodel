AC_INIT([cbsmodel], [1.0], [reese@ikp.tu-darmstadt.de])
AM_INIT_AUTOMAKE(AC_PACKAGE_NAME, AC_PACKAGE_VERSION)

AM_MAINTAINER_MODE
AM_CONFIG_HEADER(config.h)

# I don't know exactly what these are for ...
# the AC_PROG_CXX is because we are using c++
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

# we want to have these makefiles 
AC_CONFIG_FILES([
	Makefile
	src/Makefile
])
AC_OUTPUT

