#!/bin/sh

# Copyright 2010 Boris Daeppen <boris_daeppen@bluewin.ch>
# 
# This file is part of cbsmodel.
# 
# cbsmodel is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# cbsmodel is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with cbsmodel.  If not, see <http://www.gnu.org/licenses/>.

##################
# START BUILDING #
##################

# remove old packages
rm cbsmodel*.deb 2> /dev/null

##################
# BUILD CBSMODEL #
##################

echo 'START cbsmodel PACKAGE'

# pack manpage
mkdir -p debian/usr/share/man/man1
cp doc/cbsmodel.1 debian/usr/share/man/man1/
gzip --best debian/usr/share/man/man1/cbsmodel.1

#pack changelog
cp changelog debian/usr/share/doc/cbsmodel/
cp changelog.Debian debian/usr/share/doc/cbsmodel/
gzip --best debian/usr/share/doc/cbsmodel/changelog
gzip --best debian/usr/share/doc/cbsmodel/changelog.Debian

# add documentation to local folder
cd doc
pdflatex manual.tex
cd ..
cp doc/manual.pdf debian/usr/share/doc/cbsmodel/

# copy program
mkdir -p debian/usr/bin
cp cbsmodel debian/usr/bin

# update md5sums file of dep-tree
echo -e "\tupdate md5sums file"
rm debian/DEBIAN/md5sums
for i in $( find ./debian -path ./debian/DEBIAN -prune -o -type f -print)
do
    md5sum $i | sed -e "s/\.\/debian\///g" >> debian/DEBIAN/md5sums
done

# renew the size information
sed -i '/Installed-Size/ d' debian/DEBIAN/control # delete
echo "Installed-Size: $(du -s --exclude DEBIAN debian/ | cut -f1)" >> debian/DEBIAN/control

# create deb package
echo -e "\tbuild package"
fakeroot dpkg-deb --build debian \
$( grep Package debian/DEBIAN/control | cut -d" " -f2 )_\
$( grep Version debian/DEBIAN/control | cut -d" " -f2 )_\
$( grep Architecture debian/DEBIAN/control | cut -d" " -f2 )\
.deb

# remove packed things,
# I don't need it in src
rm debian/usr/share/man/man1/cbsmodel.1.gz
rm debian/usr/share/doc/cbsmodel/changelog.gz
rm debian/usr/share/doc/cbsmodel/changelog.Debian.gz
rm doc/manual.pdf


echo 'DONE'
echo "don't forget to check the package with lintian!"

