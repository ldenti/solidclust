CC=gcc
CFLAGS=-g -O3 -Wall
LDFLAGS=-lm -lz

.PHONY: all clean

all: ioc

ioc: kmer.o rsort.o main.o
	@echo "* Linking $<"
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.c
	@echo '* Compiling $<'
	$(CC) $(CFLAGS) -o $@ -c $<

clean:
	rm -rf *.o
