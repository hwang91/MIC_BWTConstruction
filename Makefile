CC = clang

CFLAGS=	-O3

setBWT: 
	$(CC) $(CFLAGS) -o setBWT setBWT.c Timing.c
 
.PHONY: clean
clean:
	rm setBWT
