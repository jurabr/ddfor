CC=gcc
DEBUG=-O3 -Wall -pedantic -DLARGE -DUI 
DEBUG=-O0 -Wall -pedantic -DLARGE -DUI -g
CFLAGS=$(DEBUG)
LIBS=-lm

all: ddfor aaem

ddfor: ddfor.c
	$(CC) $(CFLAGS) ddfor.c $(LIBS) -o $(@)

aaem: aaem.c
	$(CC) $(CFLAGS) aaem.c $(LIBS) -o $(@)

clean:
	-rm -f ddfor aaem *.o core
