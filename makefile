CFLAGS = -Wall -O3 -lm -std=c99

all: lsystems 

%: %.c 
	$(CC) $(CFLAGS) -o $@ $^

%_debug: %.c 
	$(CC) $(CFLAGS) -o $@ $^

.PHONY: clean

clean:
	rm -f lsystems lsystems_debug
