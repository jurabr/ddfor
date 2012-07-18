
CC=gcc
DEBUG=-O0 -Wall -pedantic -g 
CFLAGS=$(DEBUG)
LIBS=-lm

all: ddfor

ddfor: ddfor.c
	$(CC) $(CFLAGS) ddfor.c $(LIBS) -o $(@)

clean:
	-rm -f ddfor *.o core
