# Makefile to compile compton fitting code
#       Joshua Hoskins 
#         April 2017                                                                                                                             
#

ROOTLIBS   = $(shell root-config --libs ) -lSpectrum
ROOTGLIBS  = $(shell root-config --glibs)
LIB        = -L/usr/lib64/ -lboost_system -lboost_filesystem
INCLUDES   = -I$(shell root-config --incdir) -Iinclude/ -I/usr/include/
CC         = g++ ${INCLUDES}
SRC        = src
CFLAGS     = -O -Wall ${INCLUDES} ${LIB}

all: nmranalysis

%.o: %.cc
	${CC} ${CFLAGS} -c -o $@ $< 
nmranalysis : nmranalysis.o ${SRC}/NMRAnalysis.o 
	${CC} ${INCLUDES} -o $@  ${CFLAGS} $^ ${ROOTLIBS} ${ROOTGLIBS} ${LIB}
gainanalysis : gainanalysis.o ${SRC}/NMRAnalysis.o 
	${CC} ${INCLUDES} -o $@  ${CFLAGS} $^ ${ROOTLIBS} ${ROOTGLIBS} ${LIB}
companalysis : companalysis.o ${SRC}/NMRAnalysis.o 
	${CC} ${INCLUDES} -o $@  ${CFLAGS} $^ ${ROOTLIBS} ${ROOTGLIBS} ${LIB}
clean:
	rm -f *.o *~ src/*.o src/*~ include/*~
