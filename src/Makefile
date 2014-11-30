.PHONY: all clean

CC=gcc
LD=${CC}

GT_INCLUDE:=-DWITHOUT_CAIRO -I${GTDIR}/src
CFLAGS=-g -Wall -Werror -O3 -Wunused-parameter -Wpointer-arith ${GT_INCLUDE}
LDFLAGS=-lm
LIBS=$(GTDIR)/lib/libgenometools.a -lpthread -lm
OBJ=gt_scaffolder_graph.o test.o
OBJ2=gt_scaffolder_graph.o test_2.o

all:test.x test_2.x

gt_scaffolder_graph.o: gt_scaffolder_graph.h

test.o: gt_scaffolder_graph.h

test_2.o: gt_scaffolder_graph.h

test.x:${OBJ}
	${LD} ${OBJ} -o $@ ${LDFLAGS} $(LIBS)

test_2.x:${OBJ2}
	${LD} ${OBJ2} -o $@ ${LDFLAGS} $(LIBS)

src_check:
	${GTDIR}/scripts/src_check *.c

test: test.x
	../testsuite/gt_scaffolder.sh

clean:
	${RM} *.x *.o *~ *.dot
