This project uses the GNU autotools. I'm not very familiar with them
but the way I've done it seems to work:

The project provides the files: Makefile.am, configure.ac, src/Makefile.am
These files have very few content. Except for the configure.ac every line is 
easy to understand. 

The project can be made ready for compiling by typing the following commands:

autoreconf 				# this will complain ...
automake --add-missing  # ...  this fixes partly what 'autoreconf' is complaining about
touch NEWS AUTHORS ChangeLog  # this adds the rest of the files 'autoreconf' is complaining about 
autoreconf 				# this will do the rest

I put all these commands in the 'init.sh' file.

After that, one can compile with

./configure
make

