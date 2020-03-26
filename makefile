CFLAGS = -Wall -O3 -lm -std=c99

all: lsystems 

%: %.c 
	$(CC) $(CFLAGS) -o $@ $^

%_debug: %.c 
	$(CC) -g $(CFLAGS) -o $@ $^

.PHONY: clean
.PHONY: debug

debug: lsystems_debug

clean:
	rm -f lsystems lsystems_debug
