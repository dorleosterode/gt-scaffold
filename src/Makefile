.PHONY: all clean

CC=gcc
LD=${CC}

GT_INCLUDE:=-DWITHOUT_CAIRO -I${GTDIR}/src
CFLAGS=-g -Wall -Werror -O3 -Wunused-parameter -Wpointer-arith ${GT_INCLUDE}
LDFLAGS=-lm
LIBS=$(GTDIR)/lib/libgenometools.a -lpthread -lm
OBJ=gt_scaffolder_graph.o test.o

all:test.x

gt_scaffolder_graph.o: gt_scaffolder_graph.h 

test.o: gt_scaffolder_graph.h

test.x:${OBJ}
	${LD} ${OBJ} -o $@ ${LDFLAGS} $(LIBS)

clean:
	${RM} *.x *.o *~
