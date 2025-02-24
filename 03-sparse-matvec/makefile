# Makefile to compile mat-vec.c 
CC = gcc

all: mat-vec mat-vec-dist

mmio.o: mmio.c mmio.h
	$(CC) -c mmio.c

mat-vec: mat-vec.c mmio.o mmio.h
	$(CC) -o mat-vec mat-vec.c mmio.o -lm

mat-vec-dist: mat-vec-dist.c mmio.o mmio.h
	$(CC) -o mat-vec-dist mat-vec-dist.c mmio.o -lm

clean:
	rm -f mat-vec mmio.o
