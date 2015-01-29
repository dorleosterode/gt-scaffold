.PHONY: all clean

CC?=gcc
LD=${CC}

GT_INCLUDE:=-DWITHOUT_CAIRO -I${GTDIR}/src
CFLAGS=-g -Wall -Werror -O3 -Wunused-parameter -Wpointer-arith ${GT_INCLUDE}
LDFLAGS=-lm
LIBS=$(GTDIR)/lib/libgenometools.a -lpthread -lm
OBJ=gt_scaffolder_graph.o gt_scaffolder_parser.o gt_scaffolder_algorithms.o gt_scaffolder_bamparser.o test.o

all: src_check test.x

# generic compilation rule which creates dependency file on the fly
%.o: %.c
	$(CC) -c $< -o $@ $(CFLAGS) -MT $@ -MMD -MP -MF $(@:.o=.d)
# read dependencies
-include $(wildcard *.d)

test.x:${OBJ}
	${LD} ${OBJ} -o $@ ${LDFLAGS} $(LIBS)

src_check:
	${GTDIR}/scripts/src_check *.c

src_clean:
	${GTDIR}/scripts/src_clean *.c

test: test.x
	../testsuite/testsuite.rb

testall: test.x test
	../testsuite/gt_scaffolder_testall.sh

clean:
	${RM} *.d *.dot *.o *.x *~ *.de *.scaf *.fa
