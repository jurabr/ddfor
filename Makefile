CC=gcc
DEBUG=-O3 -Wall -pedantic -DLARGE -DUI 
DEBUG=-O0 -Wall -pedantic -DLARGE -DUI -g
CFLAGS=$(DEBUG)
LIBS=-lm

all: ddfor 

ddfor: ddfor.c
	$(CC) $(CFLAGS) ddfor.c $(LIBS) -o $(@)

dslab: dslab.c
	$(CC) $(CFLAGS) dslab.c $(LIBS) -o $(@)

dgen: dgen.c
	$(CC) $(CFLAGS) dgen.c $(LIBS) -o $(@)

ddfor_a.o: ddfor.c
	$(CC) $(CFLAGS) -DGR ddfor.c -c -o ddfor_a.o

aaem.o: aaem.c 
	$(CC) $(CFLAGS) aaem.c -c

aaem: aaem.o ddfor_a.o
	$(CC) $(CFLAGS) aaem.o ddfor_a.o $(LIBS) -o $(@)

bipe: bipe.c
	$(CC) $(CFLAGS) bipe.c $(LIBS) -o $(@)

bbet: bbet.c
	$(CC) $(CFLAGS) bbet.c $(LIBS) -o $(@)

btor: btor.c
	$(CC) $(CFLAGS) btor.c $(LIBS) -o $(@)



clean:
	-rm -f ddfor dslab aaem bipe bbet btor *.o core
