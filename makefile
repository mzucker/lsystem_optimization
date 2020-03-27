CC = gcc-mp-8
CFLAGS = -Wall -g -O3 -lm -std=c99 -fopenmp

#CC = cc
#CFLAGS = -Wall -O3 -lm -std=c99 -lpthread

all: lsystems 

%: %.c makefile
	$(CC) $(CFLAGS) -o $@ $<

%_debug: %.c 
	$(CC) -g $(CFLAGS) -o $@ $^

.PHONY: clean

debug: lsystems_debug

clean:
	rm -f lsystems lsystems_debug
