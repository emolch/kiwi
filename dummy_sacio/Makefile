CFLAGS = -std=c99 -Wall
CC = gcc

libsacio.a : sacio.o
	ar cr $@ sacio.o
	ranlib $@

sacio.o : sacio.c
	$(CC) -c $(CFLAGS)  sacio.c

clean :
	rm -f *.o libsacio.a
