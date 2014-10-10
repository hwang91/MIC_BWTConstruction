//
//	main.c
//	Construct BWT of reads collection
//	
//	Heng Wang
//	Software Research Institute 
//	National University of Defense Technology	
//	2014.10.06
//
//	Changes:
//
//
////////////////////////////////////////////////////////////////

#include "setBWT.h"

int main(int argc, char** argv){

	ParseArguments(Arguments arguments);
	LoadReadsToPack(char * readsFile, uint8_t * readsPack);

	for(prefix = 0; prefix < 256; prefix++){
		AllocMemory();
		Partition();
		AmericanFlagSort();
		OutputBWT();
		FreeMemory();
	}

	FreeReadsPack();

	return 0;
}