#! /bin/sh

aclocal
autoconf
automake --add-missing --foreign
./configure
make
