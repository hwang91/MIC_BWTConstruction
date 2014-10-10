//	readsGenerator.c
//	Generate random reads
//	Version 0.0
//
//	Heng Wang
//	National University of Defense Technology
//	2014/10/09
//
/////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define READ_LENGTH		100
#define ConvertInt2char(i) (i) == 0 ? 'A' : ((i) == 1 ? 'C' : ((i) == 2 ? 'G' : 'T'))

char singleRead[READ_LENGTH];

//unsigned baseCount[4] = {0};

void GenerateReads(int readLen, unsigned long long readsNum);


int main(int argc, char ** argv){

	//int readLen = atoi(argv[1]);
	int readLen = READ_LENGTH;
	unsigned long long readsNum = atoi(argv[1]);

	GenerateReads(readLen, readsNum);
	return 0;
}

void GenerateReads(int readLen, unsigned long long readsNum){
	unsigned long long readIndex;
	for(readIndex = 0; readIndex < readsNum; readIndex++){
		srand(readIndex);
		int baseId;
		for(baseId = 0; baseId < readLen; baseId++){
			int tmp = rand()%4;
			singleRead[baseId] = ConvertInt2char(tmp);
			//baseCount[tmp]++; 
		}
		printf("%s\n", singleRead);
	}
}


