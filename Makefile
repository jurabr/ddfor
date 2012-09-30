CC=gcc
#DEBUG=-O0 -Wall -pedantic -g -DLARGE
DEBUG=-O3 -Wall -pedantic -DLARGE
CFLAGS=$(DEBUG)
LIBS=-lm

all: ddfor

ddfor: ddfor.c
	$(CC) $(CFLAGS) ddfor.c $(LIBS) -o $(@)

clean:
	-rm -f ddfor *.o core
