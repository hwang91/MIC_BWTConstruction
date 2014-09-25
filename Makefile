CC = gcc

setBWT: 
	$(CC) -O3 -o setBWT parallel.0.0.c Timing.c
 
.PHONY: clean
clean:
	rm setBWT
