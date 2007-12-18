#! /bin/sh

aclocal
autoconf
automake --add-missing --foreign
echo "to build, run ./configure && make"
