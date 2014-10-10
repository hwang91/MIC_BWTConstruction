#include <stdio.h>
#include <inttypes.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "Timing.h"


#define MAX_LINE_LENGTH 256
#define ALPHABETA_SIZE  4
#define PRESORT_LEN     4   //Fixed value
//#define READ_LEN        100 //Should be multiple of 4
//#define READ_CODE_LEN   READ_LEN/4
//#define PARTION_COUNT   READ_NUM*READ_LEN/256 * 8
//#define READ_NUM        10000
#define PARTION_NUM     pow(ALPHABETA_SIZE, PRESORT_LEN)

#define THRESHOLD   20

typedef uint8_t * string;

typedef struct List{
    string * sa;
    int sn;
    int si;
} List;

typedef struct Arguments{
	int READ_LEN;
	unsigned long long READ_NUM;
	int READ_CODE_LEN;
	unsigned long long PARTION_COUNT;
} Arguments;


