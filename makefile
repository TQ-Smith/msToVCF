
# File: makefile
# Date: 18 December 2024
# Author: T. Quinn Smith
# Principal Investigator: Dr. Zachary A. Szpiech
# Purpose: Build msToVCF.

CC?=gcc
CFLAGS = -c -Wall -g
LFLAGS = -g -o

bin/msToVCF: src/Main.o
	mkdir -p bin
	$(CC) $(LFLAGS) bin/msToVCF src/Main.o -lz -lm

src/Main.o:
	$(CC) $(CFLAGS) src/Main.c -o src/Main.o

.PHONY: clean
clean:
	rm src/Main.o bin/msToVCF