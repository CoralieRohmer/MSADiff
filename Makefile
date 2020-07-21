CC=g++
CFLAGS=-Wall -c

all: main

main: main.o MSA.o
	$(CC) -o main main.o MSA.o


MSA.o: MSA.cpp MSA.h
	$(CC) $(CFLAGS) MSA.cpp

main.o: main.cpp MSA.h
	$(CC) $(CFLAGS) main.cpp

clean:
	rm *o

mr proper:
	rm *o main
