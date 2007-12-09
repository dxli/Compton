#general c++ compiling rule


EXECS=CPassThrough

all : ${EXECS}
.PHONY : all

CXX=g++
CXX3=g++
CC=gcc
CFLAGS =-O3 -march=nocona -fomit-frame-pointer -m64 -fopenmp -Wall
CXXFLAGS=${CFLAGS}
LDFLAGS=-Wl,-O1 -Wl,--as-needed -lgrace_np -lgsl -lgslcblas -lgomp
OBJS=compton.o block.o randomInit.o
DEPS=compton.h block.h

$(OBJS):%.o: %.cpp
		${CXX} $(CXXFLAGS) -c -o $@ $<

$(EXECS):%:%.o ${OBJS} ${DEPS}
		${CXX} -o $* $< ${OBJS} ${LDFLAGS}

clean:
		rm -f *.o ${EXECS}
