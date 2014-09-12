CC = gcc

radixSort: 
	$(CC) -o radixSort radixSort.c Timing.c
 
.PHONY: clean
clean:
	rm radixSort
