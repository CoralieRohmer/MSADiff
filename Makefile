CC=g++
CFLAGS=-Wall -c -O3

all: MTool

MTool: main.o MSA.o
	$(CC) -o MTool main.o MSA.o


MSA.o: MSA.cpp MSA.h
	$(CC) $(CFLAGS) MSA.cpp

main.o: main.cpp MSA.h
	$(CC) $(CFLAGS) main.cpp

clean:
	rm *o

mr proper:
	rm *o main
