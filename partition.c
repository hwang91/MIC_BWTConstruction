// setBWT.c
// Construct BWT for a large collection of short-reads in parallel
// Version 0.0
// Heng Wang
// National University of Defense Technology
// 2014/09/21

// TODO
/*
    1) Alloc S_Prefix dynamicly
    2) Reduce Convert4BaseToOneUint8 arguments || DONE, 2014/09/26: 10 A.M.
*/

// Changes

// 2014/09/26 
// Heng Wang
// Change type of readsPack and S_Prefix from uint8_t * to uint8_t ** || Fail

#include <stdio.h>
#include <inttypes.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/*
void *
mymtom_malloc(size_t size)
{
        void *p;
        p = malloc(size);
        printf("malloc() size=0x%lu, p=0x%lx\n", (unsigned long)size, (unsigned long)p);
        return p;
}
 
#define malloc(s) mymtom_malloc(s)

#define free(p)  do {                                                   \
                printf("%s:%d:%s:free(0x%lx)\n", __FILE__, __LINE__,    \
                    __func__, (unsigned long)p);                        \
                free(p);                                                \
        } while (0)
*/
#define MAX_LINE_LENGTH 256
#define ALPHABETA_SIZE  4

#define ConvertBase2Pack(c) c == 'A' ? 0 :(c == 'C'? 1 : (c == 'G'? 2 : (c == 'T'? 3 : 0)))


typedef uint8_t * string;

double startTime, lastEventTime, timestamp, partitionTime = 0, sortTime = 0;

int READ_LEN;
unsigned long long READ_NUM;
int READ_CODE_LEN;
//unsigned long long PARTION_COUNT;
int threadsNum;
unsigned long long partitionCount[256] = {0}; //Number of suffixes with specific prefix


inline void String2Pack(char * source, uint8_t * destination, unsigned int len);

void  GetFragment(string destination, uint8_t * source, uint8_t starPosition);

int main(int argc, char ** argv){


    READ_LEN = atoi(argv[2]);
    READ_NUM = atoi(argv[3]);
    READ_CODE_LEN = READ_LEN/4;
    //PARTION_COUNT = READ_LEN * READ_NUM / 256 * 8; //READ_LEN * READ_NUM should be larger than 256


    unsigned long long i = 0, j;

    uint8_t * readsPack = (uint8_t*)malloc(sizeof(uint8_t) * READ_CODE_LEN * READ_NUM); 
    if(readsPack == NULL){
        printf("Fail to alloc memory for sequences. Exiting ... \n");
        exit(1);
    }
    //fprintf(stderr, "\n%s %9.4f MB\n", "[Memory] Memory Allocation:", READ_CODE_LEN * READ_NUM/1048576.0);

    ///////////////////////////////////////////////////////////////////
    //  Read seuquence and convert ACGT to 0123 respetively. 
    ///////////////////////////////////////////////////////////////////
    FILE *fp = fopen(argv[1], "r");
    if(fp == NULL) {
        fprintf(stderr, "%s\n", "Fail to open file");
        exit(0);
    }


    char readLine[MAX_LINE_LENGTH];
    uint8_t readPack[MAX_LINE_LENGTH];
    //fgets(readLine, MAX_LINE_LENGTH, fp);
    while(fgets(readLine, MAX_LINE_LENGTH, fp)){
        if(readLine[0] == '>') continue;
        
        for(j = READ_LEN; j < READ_LEN + 4; j++) readLine[j] = 'A';
        String2Pack(readLine, readPack, READ_LEN+4);
        /*
        for(j = 0; j < READ_LEN; j++){
            partitionCount[readPack[j]*64 + readPack[j+1]*16 + readPack[j+2]*4 + readPack[j+3]]++;
        }
        */
        for(j = 0; j < READ_CODE_LEN; j++){
            readsPack[i * READ_CODE_LEN + j] = readPack[j*4]*64 + readPack[j*4+1]*16 + readPack[j*4+2]*4 + readPack[j*4+3];                 
        }
        i++;
    }
    fclose(fp);



    ///////////////////////////////////////////////////////////////////
    //  Print readsPack in hexadecimal way
    ///////////////////////////////////////////////////////////////////
    /*
    printf("Reads Pack in hexadecimal way:\n"); 
    for(i = 0; i < READ_NUM * READ_LEN/4; i++){
        printf("%4X", readsPack[i]);
        if(i%25 == 24) putchar('\n');
    }
    */
    ///////////////////////////////////////////////////////////////////
    //  List star position of suffixes with particular 4-prefix
    ///////////////////////////////////////////////////////////////////
    int prefix;
    string S_Prefix = 0;
    for(prefix = 0; prefix < 256; prefix++){ // For every possible 4-prefix
        
        char bwtFileName[MAX_LINE_LENGTH] = {0};
        strcat(bwtFileName, "BWT_");
        char tmp[10];
        sprintf(tmp, "%d",prefix);
        strcat(bwtFileName, tmp);
        FILE *fout = fopen(bwtFileName, "w");
        if(fout == NULL) {
            fprintf(stderr, "%s\n", "Fail to open file");
            exit(0);
        }

        S_Prefix = malloc(sizeof(*S_Prefix) * READ_CODE_LEN);
        if(S_Prefix == NULL) {
            fprintf(stderr, "%s\n", "Fail to alloc memory: S_Prefix");
            exit(0);
        }
        memset(S_Prefix, 0, READ_CODE_LEN);

        int readPos;


        ///////////////////////////////////////////////////////////////////
        //  List the suffixes with specific prefix
        ///////////////////////////////////////////////////////////////////        
        // TODO: Merge the 3 cases
        //Last 3 bases considered specially
        
        for(readPos = READ_LEN -1; readPos > READ_LEN -4; readPos--){
            for(j = 0; j < READ_NUM; j++){
                uint8_t value = (uint8_t)(readsPack[j*READ_CODE_LEN + READ_CODE_LEN-1] << (readPos%4 * 2));
                if(value == prefix){
                    //S_Prefix[0] = prefix;
                    S_Prefix[0] = ((uint8_t)readsPack[j * READ_CODE_LEN + READ_CODE_LEN-1] >> (8 - readPos%4 * 2)) & 0X03 ;
                    fwrite(S_Prefix,sizeof(uint8_t), READ_CODE_LEN, fout);
                    partitionCount[prefix]++;
                }

            }
        }
        

        for(readPos = READ_LEN - 4; readPos >= 0; readPos--){
            uint8_t byteIndex = (uint8_t)(readPos/4);
            uint8_t baseIndex = readPos % 4;
            if(!(baseIndex)) {
                for(j = 0; j < READ_NUM; j++){
                    if(readsPack[j * READ_CODE_LEN + byteIndex] == prefix) {
                        GetFragment(S_Prefix, readsPack + j * READ_CODE_LEN, readPos);
                        
                        S_Prefix[0] = (readPos == 0 ? 4 : readsPack[j * READ_CODE_LEN + byteIndex - 1] & 0x03);
                        fwrite(S_Prefix,sizeof(uint8_t), READ_CODE_LEN, fout);
                        partitionCount[prefix]++;                    
                    }
                } 
            } else {
                for(j = 0; j < READ_NUM; j++){
                    if((uint8_t)(readsPack[j * READ_CODE_LEN + byteIndex] << (baseIndex * 2)) + (uint8_t)(readsPack[j * READ_CODE_LEN + \
                                                                        byteIndex + 1] >> (8 - baseIndex * 2)) == prefix){
                        GetFragment(S_Prefix, readsPack + j * READ_CODE_LEN, readPos);

                        S_Prefix[0] = ((uint8_t)readsPack[j * READ_CODE_LEN + byteIndex] >> (8 - baseIndex * 2)) & 0X03 ;
                        fwrite(S_Prefix,sizeof(uint8_t), READ_CODE_LEN, fout);
                        partitionCount[prefix]++;
                    }
                }
            }  
        }

        free(S_Prefix); S_Prefix = NULL;
        fclose(fout);

        //fprintf(stderr, "Prefix %4d Done\n", prefix);//omp_get_thread_num());
    }

    FILE *partition = fopen("partitionCount.ini", "wb");
    if(partition == NULL) {
        fprintf(stderr, "%s\n", "Fail to open file");
        exit(0);
    }
    fwrite(partitionCount,sizeof(unsigned long long), 256, partition);
    fclose(partition);


    ///////////////////////////////////////////////////////////////////
    //  Free memory allocated
    ///////////////////////////////////////////////////////////////////
    free(readsPack);

    return 0;
}

inline void String2Pack(char * source, uint8_t * destination, unsigned int len){
    int i;
    for(i = 0; i < len; i++){
        destination[i] = ConvertBase2Pack(source[i]);
    }
}

inline void GetFragment(string destination, uint8_t * source, uint8_t starPosition){
    uint8_t mod = starPosition%4;
    uint8_t i,j;
    if(mod == 0){
        for(i = 0, j = starPosition >> 2; j < READ_CODE_LEN; i++,j++){
           destination[i] = source[j];
        }
    } else {
        for(i = 0, j = starPosition >> 2; j < READ_CODE_LEN - 1; i++,j++){
            destination[i] = (uint8_t)(source[j]<<(mod * 2)) + (uint8_t)(source[j+1] >>(8 - mod * 2));
        }
        destination[i] = (uint8_t)(source[READ_CODE_LEN-1] << (mod * 2));
    }
}

