CC = gcc-mp-8

CFLAGS = -Wall -O3 -lm -std=c99 -fopenmp

all: lsystems 

%: %.c 
	$(CC) $(CFLAGS) -o $@ $^

%_debug: %.c 
	$(CC) -g $(CFLAGS) -o $@ $^

.PHONY: clean

debug: lsystems_debug

clean:
	rm -f lsystems lsystems_debug
