vpath %.cpp src
vpath %.hpp include
vpath %.h   include

CBS = cbsmodel

OBJECTS = 	Main.o \
			CBSNucleus.o \

CPP = g++
LIBS = `gsl-config --libs` -lreadline
OPT = -O4
CPPFLAGS = -Wall -I./include/ \
			`gsl-config --cflags` \
			-O4

%.o: %.cpp %.hpp
	$(CPP) $(CPPFLAGS) $(OPT) -c $<

all: $(CBS)

$(CBS): $(OBJECTS)
	$(CPP) $(CPPFLAGS) $(OPT) -o $(CBS) $(OBJECTS) $(LIBS)
	

	
###############################################################################

.PHONY: clean

clean:
	rm -f $(OBJECTS)
