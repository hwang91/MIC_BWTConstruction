// SetBWTParallel.c
// Construct BWT for a large collection of short-reads in parallel
// Version 0.0
// Heng Wang
// National University of Defense Technology
// 2014/09/21

#include <stdio.h>
#include <inttypes.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
//#include <omp.h>

#define MAX_LINE_LENGTH 256
#define ALPHABETA_SIZE  4
#define PRESORT_LEN     4
#define READ_LEN        100
#define READ_CODE_LEN   READ_LEN/4
#define PARTION_COUNT   50000
#define READ_NUM        100
#define PARTION_NUM     pow(ALPHABETA_SIZE, PRESORT_LEN)

uint8_t Convert4BaseToOneUint8(char s1, char s2, char s3, char s4);
void  GetFragment(uint8_t * destination, uint8_t * source, uint8_t starPosition);

int main(int argc, char ** argv){

/*
    int READ_LEN = atoi(argv[2]);
    int READ_NUM = atoi(argv[3]);
*/

    int i = 0, j;

    uint8_t * s = (uint8_t*)malloc(sizeof(uint8_t) * READ_CODE_LEN * READ_NUM); //4: to be replaced by formula of ALPHABETA_SIZE
    if(s == NULL){
        printf("Fail to alloc memory for sequences. Exiting ... \n");
        exit(1);
    }

    ///////////////////////////////////////////////////////////////////
    //  Read seuquence and convert ACGT to 0123 respetively. 
    ///////////////////////////////////////////////////////////////////
    FILE *fp = fopen(argv[1], "r");
    if(fp == NULL) {
        fprintf(stderr, "%s\n", "Fail to open file");
        exit(0);
    }

    char readLine[MAX_LINE_LENGTH];
    //fgets(readLine, MAX_LINE_LENGTH, fp);
    while(fgets(readLine, MAX_LINE_LENGTH, fp)){
        if(readLine[0] == '>') continue;
        for(j = 0; j < READ_CODE_LEN; j++){
            s[i * READ_CODE_LEN + j] = Convert4BaseToOneUint8(readLine[j*4], readLine[j*4+1], readLine[j*4+2], readLine[j*4+3]);                 
        }
        i++;
    }
    fclose(fp);

    ///////////////////////////////////////////////////////////////////
    //  Print s to screen for check
    ///////////////////////////////////////////////////////////////////
  
    for(i = 0; i < READ_NUM * READ_LEN/4; i++){
        printf("%4X", s[i]);
        if(i%25 == 24) putchar('\n');
    }

    ///////////////////////////////////////////////////////////////////
    //  List star position of suffixes with particular 4-prefix
    ///////////////////////////////////////////////////////////////////
    int prefix;
    for(prefix = 0; prefix < 256; prefix++){ // For every possible 4-prefix

        uint8_t * S_Prefix = (uint8_t*)malloc(sizeof(uint8_t)*(READ_CODE_LEN)*PARTION_COUNT);
        if(S_Prefix == NULL) {
            fprintf(stderr, "%s\n", "Fail to alloc memory: S_Prefix");
            exit(0);
        }
    
        memset(S_Prefix, 0, READ_CODE_LEN*PARTION_COUNT);
    

        int S_Prefix_index = 0;

        int readPos;

        for(readPos = READ_LEN - 4; readPos >= 0; readPos--){
            if(!(readPos%4)) {
                for(j = 0; j < READ_NUM; j++){
                    if(s[j * READ_CODE_LEN + readPos/4] == prefix) {
                        GetFragment(S_Prefix + (READ_CODE_LEN)*S_Prefix_index, s + j * READ_CODE_LEN, readPos);
                        S_Prefix[READ_CODE_LEN*S_Prefix_index] = (readPos == 0 ? 4 : s[j * READ_CODE_LEN + readPos/4 - 1] & 0x03);
                        S_Prefix_index++;
                    }
                } 
            } else {
                for(j = 0; j < READ_NUM; j++){
                    if((uint8_t)(s[j * READ_CODE_LEN + readPos/4] << (readPos%4 * 2)) + (uint8_t)(s[j * READ_CODE_LEN + \
                                                                        readPos/4 + 1] >> (8 - readPos%4 * 2)) == prefix){
                        GetFragment(S_Prefix + (READ_CODE_LEN)*S_Prefix_index, s + j * READ_CODE_LEN, readPos);
                        uint8_t sw = readPos % 4 - 1;
                        switch (sw){
                            case 0 : S_Prefix[READ_CODE_LEN*S_Prefix_index] = (uint8_t)(s[j * READ_CODE_LEN + readPos/4] & 0xC0 >> 6); break;
                            case 1 : S_Prefix[READ_CODE_LEN*S_Prefix_index] = (uint8_t)(s[j * READ_CODE_LEN + readPos/4] & 0x30 >> 4); break;
                            case 2 : S_Prefix[READ_CODE_LEN*S_Prefix_index] = (uint8_t)(s[j * READ_CODE_LEN + readPos/4] & 0x0C >> 2); break;
                        }
                        S_Prefix_index++;
                    }
                }
            }  
        }

        //uint8_t * S_Prefix_tmp = realloc(S_Prefix, sizeof(uint8_t) * S_Prefix_index * READ_CODE_LEN);
        //if(!S_Prefix_tmp) S_Prefix = S_Prefix_tmp;

        for(i = 0; i < S_Prefix_index * (READ_CODE_LEN); i++){
            printf("%4X", S_Prefix[i]);
            if(i%(READ_CODE_LEN) == READ_CODE_LEN-1) putchar('\n');
        }
        putchar('\n');
        free(S_Prefix);

    }

    ///////////////////////////////////////////////////////////////////
    //  Free memory allocated
    ///////////////////////////////////////////////////////////////////
    //free(S_Prefix);

    fprintf(stderr, "%s", "Free s ...");
    free(s);
    fprintf(stderr, "%s\n", "Done.");
    return 0;
}

uint8_t Convert4BaseToOneUint8(char s1, char s2, char s3, char s4){
    uint8_t i1 = (s1 == 'A' ? 0 : (s1 == 'C' ? 1 : (s1 == 'G'? 2 : 3)));
    uint8_t i2 = (s2 == 'A' ? 0 : (s2 == 'C' ? 1 : (s2 == 'G'? 2 : 3)));
    uint8_t i3 = (s3 == 'A' ? 0 : (s3 == 'C' ? 1 : (s3 == 'G'? 2 : 3)));
    uint8_t i4 = (s4 == 'A' ? 0 : (s4 == 'C' ? 1 : (s4 == 'G'? 2 : 3)));
    return i1 * 64 + i2 * 16 + i3 * 4 + i4;
}

void GetFragment(uint8_t* destination, uint8_t * source, uint8_t starPosition){
    uint8_t i,j;
    if(starPosition%4 == 0){
        for(i = 0, j = starPosition/4; j < READ_CODE_LEN; i++,j++){
           destination[i] = source[j];
        }
        //strncpy(destination, source + starPosition/4, READ_CODE_LEN - starPosition/4);
    } else {
        for(i = 0, j = starPosition/4; j < READ_CODE_LEN - 1; i++,j++){
            destination[i] = (uint8_t)(source[j]<<(starPosition%4 * 2)) + (uint8_t)(source[j+1] >>(8 - starPosition%4 * 2));
        }
        destination[i] = (uint8_t)(source[READ_CODE_LEN-1] << (starPosition%4 * 2));
    }
    return;
}