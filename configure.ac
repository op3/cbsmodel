AC_INIT([cbsmodel], [1.0], [You@that.place])
AM_INIT_AUTOMAKE(AC_PACKAGE_NAME, AC_PACKAGE_VERSION)

AM_MAINTAINER_MODE
AM_CONFIG_HEADER(config.h)

AC_ISC_POSIX
AC_PROG_CXX
AM_PROG_CC_STDC
AC_HEADER_STDC

# check for GNU scientific library
AC_CHECK_LIB([m],[cos])
AC_CHECK_LIB([gslcblas],[cblas_dgemm])
AC_CHECK_LIB([gsl],[gsl_blas_dgemm])

# check for readline library
AC_CHECK_LIB([readline], [main])

AC_CONFIG_FILES([
	Makefile
	src/Makefile
])
AC_OUTPUT

