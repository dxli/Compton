#!/bin/bash

aclocal
libtoolize --copy --force --automake
autoconf
aclocal
autoheader
automake --add-missing --foreign
echo "to build, run: ./configure && make"
