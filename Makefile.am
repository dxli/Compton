#general c++ compiling rule
AM_CPPFLAGS = -fomit-frame-pointer -fexpensive-optimizations
bin_PROGRAMS = CPassThrough CPassThrough2D CRefBuried CPTBrutal binMatrix
CRefBuried_SOURCES = CRefBuried.cpp
CPTBrutal_SOURCES = CPTBrutal.cpp
CPTBrutal_LDADD =  compton.o blockBrutal.o randomInit.o
LDADD = compton.o block.o randomInit.o
include_HEADERS=compton.h block.h
