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

#define MAX_LINE_LENGTH 256
#define ALPHABETA_SIZE  4
#define PRESORT_LEN     4   //Fixed value
//#define READ_LEN        100 //Should be multiple of 4
//#define READ_CODE_LEN   READ_LEN/4
//#define PARTION_COUNT   READ_NUM*READ_LEN/256 * 8
//#define READ_NUM        10000
#define PARTION_NUM     pow(ALPHABETA_SIZE, PRESORT_LEN)

#define THRESHOLD   15

int READ_LEN;
int READ_NUM;
int READ_CODE_LEN;
int PARTION_COUNT;

typedef uint8_t * string;
typedef struct List{
    string * sa;
    int sn;
    int si;
} List;

void rsort(string *, int);
void SelectionSort(string *a, int n, int b);
uint8_t Convert4BaseToOneUint8(char s1, char s2, char s3, char s4);
void  GetFragment(uint8_t * destination, uint8_t * source, uint8_t starPosition);

int main(int argc, char ** argv){


    READ_LEN = atoi(argv[2]);
    READ_NUM = atoi(argv[3]);
    READ_CODE_LEN = READ_LEN/4;
    PARTION_COUNT = READ_LEN * READ_NUM / 256 * 8; //READ_LEN * READ_NUM should be larger than 256


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
    //  Print s in hexadecimal way
    ///////////////////////////////////////////////////////////////////
    /*
    for(i = 0; i < READ_NUM * READ_LEN/4; i++){
        printf("%4X", s[i]);
        if(i%25 == 24) putchar('\n');
    }
    */
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

        ///////////////////////////////////////////////////////////////////
        //  List the suffixes with specific prefix
        ///////////////////////////////////////////////////////////////////        
        // TODO: Merge the 3 cases
        //Last 3 bases considered specially
        for(readPos = READ_LEN -1; readPos > READ_LEN -4; readPos--){
            for(j = 0; j < READ_NUM; j++){
                uint8_t value = (uint8_t)(s[j*READ_CODE_LEN + READ_CODE_LEN-1] << (readPos%4 * 2));
                if(value == prefix){
                    //S_Prefix[READ_CODE_LEN*S_Prefix_index] = value;
                    
                    uint8_t sw = readPos % 4 - 1;
                    switch (sw){
                        case 0 : S_Prefix[READ_CODE_LEN*S_Prefix_index] = ((uint8_t)(s[j * READ_CODE_LEN + readPos/4] & 0xC0) >> 6); break;
                        case 1 : S_Prefix[READ_CODE_LEN*S_Prefix_index] = ((uint8_t)(s[j * READ_CODE_LEN + readPos/4] & 0x30) >> 4); break;
                        case 2 : S_Prefix[READ_CODE_LEN*S_Prefix_index] = ((uint8_t)(s[j * READ_CODE_LEN + readPos/4] & 0x0C) >> 2); break;
                    }
                    
                    S_Prefix_index++;
                    if(S_Prefix_index > PARTION_COUNT){ // TODO: realloc
                        fprintf(stderr, "[ERROR] Number of %4d-suffixes overflows\n", prefix);
                        exit(0);
                    }
                }

            }
        }
        
        for(readPos = READ_LEN - 4; readPos >= 0; readPos--){
            if(!(readPos%4)) {
                for(j = 0; j < READ_NUM; j++){
                    if(s[j * READ_CODE_LEN + readPos/4] == prefix) {
                        GetFragment(S_Prefix + (READ_CODE_LEN)*S_Prefix_index, s + j * READ_CODE_LEN, readPos);
                        S_Prefix[READ_CODE_LEN*S_Prefix_index] = (readPos == 0 ? 4 : s[j * READ_CODE_LEN + readPos/4 - 1] & 0x03);
                        S_Prefix_index++;
                        if(S_Prefix_index > PARTION_COUNT){
                            fprintf(stderr, "[ERROR] Number of %4d-suffixes overflows\n", prefix);
                        }
                    }
                } 
            } else {
                for(j = 0; j < READ_NUM; j++){
                    if((uint8_t)(s[j * READ_CODE_LEN + readPos/4] << (readPos%4 * 2)) + (uint8_t)(s[j * READ_CODE_LEN + \
                                                                        readPos/4 + 1] >> (8 - readPos%4 * 2)) == prefix){
                        GetFragment(S_Prefix + (READ_CODE_LEN)*S_Prefix_index, s + j * READ_CODE_LEN, readPos);
                        
                        uint8_t sw = readPos % 4 - 1;
                        switch (sw){
                            case 0 : S_Prefix[READ_CODE_LEN*S_Prefix_index] = ((uint8_t)(s[j * READ_CODE_LEN + readPos/4] & 0xC0) >> 6); break;
                            case 1 : S_Prefix[READ_CODE_LEN*S_Prefix_index] = ((uint8_t)(s[j * READ_CODE_LEN + readPos/4] & 0x30) >> 4); break;
                            case 2 : S_Prefix[READ_CODE_LEN*S_Prefix_index] = ((uint8_t)(s[j * READ_CODE_LEN + readPos/4] & 0x0C) >> 2); break;
                        }
                        
                        S_Prefix_index++;
                        if(S_Prefix_index > PARTION_COUNT){
                            fprintf(stderr, "[ERROR] Number of %4d-suffixes overflows\n", prefix);
                        }                        
                    }
                }
            }  
        }

        //fprintf(stderr,"%d\t%d\n", prefix, S_Prefix_index);

        ///////////////////////////////////////////////////////////////////
        //  Print all suffixes, in order of 0000 to 3333
        ///////////////////////////////////////////////////////////////////
        /*
        for(i = 0; i < S_Prefix_index * (READ_CODE_LEN); i++){
            printf("%4X", S_Prefix[i]);
            if(i%(READ_CODE_LEN) == READ_CODE_LEN-1) putchar('\n');
        }
        */
        string * aa = (string*)malloc(sizeof(string)*S_Prefix_index);
        for(i = 0; i < S_Prefix_index; i++)
            aa[i] = S_Prefix + READ_CODE_LEN * i;

        rsort(aa, S_Prefix_index);
        
        ///////////////////////////////////////////////////////////////////
        //  Print sorted suffixes
        ///////////////////////////////////////////////////////////////////
        /*
        for(i = 0; i < S_Prefix_index; i++){
            for(j = 0; j < READ_CODE_LEN; j++)
                printf("%4X", *(*(aa + i) + j));
            printf("\n");
        }
        */

        ///////////////////////////////////////////////////////////////////
        //  Print BWT
        ///////////////////////////////////////////////////////////////////
        for(i = 0; i < S_Prefix_index; i++){
            printf("%u", *(aa+i)[0]);
        }
        

        free(aa);
        free(S_Prefix);

    }
    printf("\n");

    ///////////////////////////////////////////////////////////////////
    //  Free memory allocated
    ///////////////////////////////////////////////////////////////////
    free(s);
    return 0;
}

//Convert 4 bases to an uint8_t
uint8_t Convert4BaseToOneUint8(char s1, char s2, char s3, char s4){
    uint8_t i1 = (s1 == 'A' ? 0 : (s1 == 'C' ? 1 : (s1 == 'G'? 2 : 3)));
    uint8_t i2 = (s2 == 'A' ? 0 : (s2 == 'C' ? 1 : (s2 == 'G'? 2 : 3)));
    uint8_t i3 = (s3 == 'A' ? 0 : (s3 == 'C' ? 1 : (s3 == 'G'? 2 : 3)));
    uint8_t i4 = (s4 == 'A' ? 0 : (s4 == 'C' ? 1 : (s4 == 'G'? 2 : 3)));
    return i1 * 64 + i2 * 16 + i3 * 4 + i4;
}

void GetFragment(uint8_t* destination, uint8_t * source, uint8_t starPosition){
    uint8_t i,j;
    if(starPosition > READ_LEN - 4){
        destination[0] =  (uint8_t)(source[READ_CODE_LEN-1] << (starPosition%4 * 2));
    }
    else if(starPosition%4 == 0){
        for(i = 0, j = starPosition/4; j < READ_CODE_LEN; i++,j++){
           destination[i] = source[j];
        }
        //strncpy(destination, source + starPosition/4, READ_CODE_LEN - starPosition/4);
    } 
    else {
        for(i = 0, j = starPosition/4; j < READ_CODE_LEN - 1; i++,j++){
            destination[i] = (uint8_t)(source[j]<<(starPosition%4 * 2)) + (uint8_t)(source[j+1] >>(8 - starPosition%4 * 2));
        }
        destination[i] = (uint8_t)(source[READ_CODE_LEN-1] << (starPosition%4 * 2));
    }
}

void rsort(string * a, int n){ //Sort n strings

    List *stack = (List*)malloc(sizeof(List) * READ_NUM);
    List *sp = stack; 
    string          *pile[256], *ai, *ak, *ta;
    static int      count[256] = {0};
    int             b = 1, c, cmin, *cp, nc = 0; // nc: number of unempty buckets

#define push(a, n, i)   sp->sa = a, sp->sn = n, (sp++)->si = i
#define pop(a, n, i)    a = (--sp)->sa, n = sp->sn, i = sp->si
#define stackempty()    (sp <= stack)


    ta = malloc(n*sizeof(string)); // Total size: n
    push(a, n, 1); // sp->sa = a, sp->sn = n, sp->si = 0, sp++; 
    while(sp > stack) {
        pop(a, n, b); // sp--, a = sp->sa, n = sp->sn, b = sp->si;
        //Ignore the bytes before b
    
        //When the total num is less thah THRESHOLD, sort them by means of comparison.
        if(n < THRESHOLD) {
            SelectionSort(a, n, b);
            continue;
        }
    
        cmin = 255; //Minimum value in the bytes
        for(ak = a+n; --ak >= a; ) {
            c = (*ak)[b]; //Fetch the b-th byte of ak-th string
            if(++count[c] == 1 && c > 0) { // count[c]++; then check if c is a new value
                if(c < cmin) cmin = c; //Record the minimum value
                nc++; //New bucket
            }
        }

        pile[0] = ak = a + count[0];         
        count[0] = 0;
        for(cp = count+cmin; nc > 0; cp++, nc--) {
            while(*cp == 0) cp++; //find next cp s.t. cp != 0;
            if(*cp > 1)
                push(ak, *cp, b+1); // Sort ak by the (b+1)-th byte, with total num cp
            pile[cp - count] = ak += *cp;
            *cp = 0;
        }

        for(ak = ta+n, ai = a+n; ak > ta; )
            *--ak = *--ai;
        for(ak = ta+n; ak-- > ta; )
            *--pile[(*ak)[b]] = *ak;
    } 
    free(stack);
    free(ta);
}

void SelectionSort(string *a, int n, int b) {
    string ak;
    int base = b;
    int i, j, m;
    for (i = 0; i < n; i++) {
        for (j = i, m = i; j < n; j++) {
            b = base;
            while(b < READ_CODE_LEN && (*(a+j))[b] == (*(a+m))[b]){
                b++;
            }
            if ((*(a+j))[b] < (*(a+m))[b])
                m = j;
        }
        ak = a[i];
        a[i] = a[m];
        a[m] = ak;
    }
}