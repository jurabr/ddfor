CC=gcc
DEBUG=-O3 -Wall -pedantic -DLARGE -DUI -DFREEROT
DEBUG=-O0 -Wall -pedantic -DLARGE -DUI -DFREEROT -g
CFLAGS=$(DEBUG)
LIBS=-lm

all: ddfor

ddfor: ddfor.c
	$(CC) $(CFLAGS) ddfor.c $(LIBS) -o $(@)

clean:
	-rm -f ddfor *.o core
