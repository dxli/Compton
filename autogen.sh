#! /bin/sh

aclocal
autoconf
autoheader
automake --add-missing --foreign
echo "to build, run ./configure && make"
