CC = gcc

radixSort: 
	$(CC) -o radixSort radixSort.c
 
.PHONY: clean
clean:
	rm radixSort