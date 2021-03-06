CC = gcc-mp-8
CFLAGS = -Wall -O3 -lm -std=c99 -fopenmp

#CC = cc
#CFLAGS = -Wall -O3 -lm -std=c99 -lpthread

PROGS = lsystems_v2 lsystems_v3 lsystems_v4 lsystems_v5 lsystems_v6 lsystems_v7
DEBUG_PROGS = $(PROGS:=_debug)

all: $(PROGS)

%: %.c makefile
	$(CC) $(CFLAGS) -o $@ $<

%_debug: %.c makefile
	$(CC) -g $(CFLAGS) -o $@ $<

.PHONY: clean

debug: $(DEBUG_PROGS)

clean:
	rm -f $(PROGS) $(DEBUG_PROGS)
