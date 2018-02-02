CC=gcc
I=../include
# LIB=-L../lib
CFLAGS=-Wall -O3
BINDIR=../../bin

tephra2:	new_tephra.o parameters_strat.h tephra_calc.o
		$(CC) $(CFLAGS) $(LIB) -o tephra2  new_tephra.o tephra_calc.o -lgc -lm
		install tephra2 $(BINDIR)

new_tephra.o:	new_tephra.c common_structures_strat.h makefile
		$(CC) $(CFLAGS) -I$(I) -c new_tephra.c

tephra_calc.o:		tephra_calc.c common_structures_strat.h makefile
		$(CC) $(CFLSGS) -I$(I) -c tephra_calc.c

clean:
		rm -fv *.o *.bak tephra2
