## Additional flags to pass to aclocal when it is invoked automatically at
## make time. The ${ACLOCAL_FLAGS} variable is picked up from the environment
## to provide a way for the user to supply additional arguments.
ACLOCAL_AMFLAGS = ${ACLOCAL_FLAGS}

## Define an executable target "cbsmodel", which will be installed into the
## directory named by the predefined variable $(bindir).
bin_PROGRAMS = cbsmodel

## Define the list of source files for the "cbsmodel" target. The file extension
## .cpp is recognized by Automake, and causes it to produce rules which invoke
## the C++ compiler to produce an object file (.o) from each source file. The
## header files (.h) do not result in object files by themselves, but will be
## included in distribution archives of the project.
cbsmodel_SOURCES = src/CBSNucleus.hpp src/CBSNucleus.cpp \
                   src/Main.cpp src/ src/fit.hpp

## Define an independent executable script for inclusion in the distribution
## archive. It will not be installed on an end user's system, however.
dist_noinst_SCRIPTS = autogen.sh           \
                      doc/manual.tex       \
                      doc/r_beta_chart.pdf \
                      doc/manual.pdf
