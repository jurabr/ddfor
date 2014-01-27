CC=gcc
DEBUG=-O3 -Wall -pedantic -DLARGE -DUI 
DEBUG=-O0 -Wall -pedantic -DLARGE -DUI -g
CFLAGS=$(DEBUG)
LIBS=-lm

all: ddfor aaem

ddfor: ddfor.c
	$(CC) $(CFLAGS) ddfor.c $(LIBS) -o $(@)

ddfor_a.o: ddfor.c
	$(CC) $(CFLAGS) -DGR ddfor.c -c -o ddfor_a.o

aaem.o: aaem.c 
	$(CC) $(CFLAGS) aaem.c -c

aaem: aaem.o ddfor_a.o
	$(CC) $(CFLAGS) aaem.o ddfor_a.o $(LIBS) -o $(@)

clean:
	-rm -f ddfor aaem *.o core
