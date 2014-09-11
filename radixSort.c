//radixSort.c
// Wang Heng
// Sept 3, 2014

//TODO
//  1) Change filename to BWTConstruct.c
//  2) Get Elapsed time for each module
//  3) Move some codes out of the main fuction
//  4) Parallel on CPU
//  5) MIC version
//  6) Add the last column of original sequences to the BWT
//  7) SSE instructors

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

//#define PRINT_PASS
//#define PRINT_COUNT
//#define DEBUG

#define ALPHABETA_SIZE  4
#define MAX_LINE_LENGTH 256
#define WRITE_BUFFER_SIZE   1024
//#define READ_LEN        4
//#define READ_NUM        3


//char s[READ_NUM][READ_LEN];
//int  Count[READ_LEN][ALPHABETA_SIZE + 1] = {0}; 
char ConvertBase2Int(char c); //Convert ACGT to 1234;

int main(int argc, char** argv){

    int READ_LEN = atoi(argv[2]);
    int READ_NUM = atoi(argv[3]);

    char* s = (char*)malloc(sizeof(char) * READ_LEN * READ_NUM);
    if(s == NULL){
        printf("Fail to alloc memory for sequences. Exiting ... \n");
        exit(1);
    }
    //memset(s, 32, READ_LEN * READ_NUM);

    int * Count = (int*)malloc(sizeof(int) * READ_LEN * (ALPHABETA_SIZE + 1));
    if(Count == NULL){
        printf("Fail to alloc memory for Count. Exiting ... \n");
        exit(1);
    }


    int i = 0;
    int j = 0;
    FILE *fp, *fp_out;
    char OutputFileName[MAX_LINE_LENGTH];
    strcpy(OutputFileName, argv[1]);
    strcat(OutputFileName,".bwt");

    if((fp = fopen(argv[1], "rt")) == NULL){
        printf("Cannot open file strike any key exit!");
        exit(1);
    }

    if((fp_out = fopen(OutputFileName, "w")) == NULL){
        printf("Cannot open Output file strike any key exit!");
        exit(1);
    }


    char readLine[MAX_LINE_LENGTH];
    //fgets(readLine, MAX_LINE_LENGTH, fp);
    while(fgets(readLine, MAX_LINE_LENGTH, fp)){
        if(readLine[0] == '>') continue;
        for(j = 0; j < READ_LEN; j++){
            char c = readLine[j];
            if(c != 'A' && c != 'C' && c != 'G' && c != 'T'){
                fprintf(stderr, "%c\t%d\t%d\n", c, i, j);
            }
            char tmp = ConvertBase2Int(c);
            s[i * READ_LEN + j] = tmp;
            Count[j * (ALPHABETA_SIZE + 1) + (tmp - '0')]++;
        }
        i++;
        //fgets(readLine, MAX_LINE_LENGTH, fp);
    }

/*
    //Print sequence code
    for(i = 0; i < READ_NUM; i++){
        for(j = 0; j < READ_LEN; j++)
            putchar(s[i * READ_LEN + j]);
        putchar('\n');
    }
*/
    for(i = READ_LEN - 2; i >= 0; i--){
        for(j = 0; j < ALPHABETA_SIZE + 1; j++)
            Count[i * (ALPHABETA_SIZE + 1) + j] += Count[(i + 1) * (ALPHABETA_SIZE + 1) + j];
    }

    // Calculate # of 0 
    for(j = 0; j < READ_LEN; j++)
        Count[j * (ALPHABETA_SIZE + 1)] = READ_LEN * READ_NUM  - Count[j * (ALPHABETA_SIZE + 1) + 1] \
         - Count[j * (ALPHABETA_SIZE + 1) + 2] - Count[j * (ALPHABETA_SIZE + 1) + 3] - Count[j * (ALPHABETA_SIZE + 1) + 4];

    //Accumulation
    for(i = 0; i < READ_LEN; i++){
        for(j = 1; j < ALPHABETA_SIZE + 1; j++)
            Count[i * (ALPHABETA_SIZE + 1) + j] += Count[i * (ALPHABETA_SIZE + 1) + j - 1];
    }

#ifdef PRINT_COUNT
    for(i = 0; i < READ_LEN; i++){
        printf("%d\t", i );
        for(j = 0; j < ALPHABETA_SIZE + 1; j++)
            printf("%d\t", Count[i * (ALPHABETA_SIZE + 1) + j]);
        putchar('\n');
    }
#endif //PRINT_COUNT


    int * rank_before = (int*)malloc(sizeof(int) * READ_NUM * READ_LEN);
    if(rank_before == NULL){
        printf("Fail to alloc memory for rank_before. Exiting ... \n");
        exit(1);
    }
 
    int * rank_after  = (int*)malloc(sizeof(int) * READ_NUM * READ_LEN);
    if(rank_after == NULL){
        printf("Fail to alloc memory for rank_after. Exiting ... \n");
        exit(1);
    }

    for(i = 0; i < READ_LEN * READ_NUM; i++){
        rank_after[i] = i;
    }

    //LSD Radix Sort
    for(i = READ_LEN - 1; i >= 0; i--){
        //Count sort
        int k, value, pos;

        // set rank_before to rank_after
        for(k = 0; k < READ_LEN * READ_NUM ; k++) 
            rank_before[k] = rank_after[k];

        //Counting Sort
        for(k = READ_NUM * READ_LEN - 1; k >= 0; k--){
            #ifdef DEBUG
                printf("%s %d\n", "Debug", k);
            #endif
            value = ((rank_before[k] % READ_LEN + i) >= READ_LEN )? 0 : s[rank_before[k] + i] - '0'; 
            if(value > 4 || value < 0){
                printf("Error. %d  %d  %d\n", value, i, k);
                exit(0);
            }
            pos   = Count[i * (ALPHABETA_SIZE + 1) + value]; 
            rank_after[pos - 1] = rank_before[k]; 
            Count[i * (ALPHABETA_SIZE + 1) + value]--; 
        }
        
    #ifdef PRINT_PASS
        for(k = 0; k < READ_LEN * READ_NUM; k++)
            printf("%d\t", rank_after[k]);
        printf("\n");
    #endif

    }

    free(Count);
    free(rank_before); 

    fprintf(stderr, "%s", "Printing BWT ... ");
    
    //Print BWT
    char * BWT = (char*)malloc(sizeof(char) * READ_LEN * READ_NUM);

    for(i = 0; i < READ_LEN * READ_NUM; i++){
        if(rank_after[i] % READ_LEN == 0) {
            BWT[i] = '0';
            continue;
        }
        BWT[i] = s[rank_after[i] - 1];
    }

    for(i = 0; i < READ_NUM; i++){
        fwrite(BWT + i * READ_LEN, 1, READ_LEN, fp_out);
        fwrite("\n",1,1,fp_out);
    }

    printf("Done\n");
    free(BWT);
    

 
    free(rank_after);
    free(s);
    return 0;
}
char ConvertBase2Int(char c){
    char tmp;
    switch(c){
        case 'A': tmp = '1'; break;
    //    case 'a': tmp = '1'; break;
        case 'C': tmp = '2'; break;
    //    case 'c': tmp = '2'; break;
        case 'G': tmp = '3'; break;
    //    case 'g': tmp = '3'; break;
        case 'T': tmp = '4'; break;
    //    case 't': tmp = '4'; break;
        case 'N': tmp = '1'; break;
        default : tmp = '1'; 
    }
    return tmp;
}
