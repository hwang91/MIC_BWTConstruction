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

#include "setBWT.h"

double startTime, lastEventTime, timestamp, partitionTime = 0, sortTime = 0;

int READ_LEN;
unsigned long long READ_NUM;
int READ_CODE_LEN;
//unsigned long long PARTION_COUNT;
int threadsNum;
unsigned long long partitionCount[256] = {0}; //Number of suffixes with specific prefix


inline void String2Pack(char * source, uint8_t * destination, unsigned int len);

void rsort(string *, string *, List*, unsigned long long);
void SelectionSort(string *a, unsigned long long n, int b);
uint8_t Convert4BaseToOneUint8(char*, int);
void  GetFragment(string destination, uint8_t * source, uint8_t starPosition);

int main(int argc, char ** argv){


    READ_LEN = atoi(argv[2]);
    READ_NUM = atoi(argv[3]);
    threadsNum = atoi(argv[4]);
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
/*
    FILE *fout = fopen(argv[argc - 1], "w");
    if(fout == NULL) {
        fprintf(stderr, "%s\n", "Fail to open file");
        exit(0);
    }
*/

    startTime = setStartTime();
    lastEventTime = 0;


    char readLine[MAX_LINE_LENGTH];
    uint8_t readPack[MAX_LINE_LENGTH];
    //fgets(readLine, MAX_LINE_LENGTH, fp);
    while(fgets(readLine, MAX_LINE_LENGTH, fp)){
        if(readLine[0] == '>') continue;
        for(j = READ_LEN; j < READ_LEN + 4; j++) readLine[j] = 'A';
        String2Pack(readLine, readPack, READ_LEN+4);
        
        for(j = 0; j < READ_LEN; j++){
            partitionCount[readPack[j]*64 + readPack[j+1]*16 + readPack[j+2]*4 + readPack[j+3]]++;
        }

        for(j = 0; j < READ_CODE_LEN; j++){
            readsPack[i * READ_CODE_LEN + j] = readPack[j*4]*64 + readPack[j*4+1]*16 + readPack[j*4+2]*4 + readPack[j*4+3];                 
        }
        i++;
    }
    fclose(fp);

    //for(i = 0; i < 256; i++)
     //   fprintf(stderr, "Prefix %d : %llu\n", i, partitionCount[i]);

    timestamp = getElapsedTime(startTime);
    fprintf(stderr, "(Elapsed time (Sequence Loading) : %9.4f seconds)\n", timestamp - lastEventTime);
    lastEventTime = timestamp;


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
    string * S_Prefix = 0;
    string * ta = 0;
    uint8_t * prefixMemory = 0;
    List * stack = 0;
    #pragma omp parallel for num_threads(threadsNum) reduction(+:partitionTime, sortTime) \
            private(prefix, S_Prefix, prefixMemory, stack, ta, i, j) shared(readsPack) schedule(dynamic)
    for(prefix = 0; prefix < 256; prefix++){ // For every possible 4-prefix
        //fprintf(stderr, "Prefix %4d on thread %d\n", prefix, omp_get_thread_num());

        //uint8_t * S_Prefix = (uint8_t*)malloc(sizeof(uint8_t)*(READ_CODE_LEN)*PARTION_COUNT);
        S_Prefix = malloc(sizeof(*S_Prefix) * partitionCount[prefix]);
        if(S_Prefix == NULL) {
            fprintf(stderr, "%s\n", "Fail to alloc memory: S_Prefix");
            exit(0);
        }
        prefixMemory = malloc(sizeof(*prefixMemory)*READ_CODE_LEN* partitionCount[prefix]);
        for(i = 0; i < partitionCount[prefix]; i++){
            S_Prefix[i] = prefixMemory + i * READ_CODE_LEN;
            if(S_Prefix[i] == NULL){
                fprintf(stderr, "Fail to alloc memory: S_Prefix[%llu]\n", i);
            }
            memset(S_Prefix[i], 0, READ_CODE_LEN);
        }
        //fprintf(stderr, "Memory allocated %9.4f on thread %d\n", \
        //    (sizeof(*S_Prefix)*partitionCount[prefix] + sizeof(uint8_t)*READ_CODE_LEN*partitionCount[prefix])/1048576.0, omp_get_thread_num());
    
        //memset(S_Prefix, 0, READ_CODE_LEN*PARTION_COUNT);
    

        unsigned long long S_Prefix_index = 0;

        int readPos;

        ///////////////////////////////////////////////////////////////////
        //  List the suffixes with specific prefix
        ///////////////////////////////////////////////////////////////////        
        // TODO: Merge the 3 cases
        //Last 3 bases considered specially

        lastEventTime = getElapsedTime(startTime);
        
        for(readPos = READ_LEN -1; readPos > READ_LEN -4; readPos--){
            for(j = 0; j < READ_NUM; j++){
                uint8_t value = (uint8_t)(readsPack[j*READ_CODE_LEN + READ_CODE_LEN-1] << (readPos%4 * 2));
                if(value == prefix){
                    //S_Prefix[S_Prefix_index][0] = prefix;
                    //S_Prefix[READ_CODE_LEN*S_Prefix_index] = ((uint8_t)readsPack[j * READ_CODE_LEN + READ_CODE_LEN-1] >> (8 - readPos%4 * 2)) & 0X03 ;
                    uint8_t tmp = (uint8_t)(((uint8_t)readsPack[j * READ_CODE_LEN + READ_CODE_LEN-1] >> (8 - readPos%4 * 2)) & 0X03);
                    (*(S_Prefix + S_Prefix_index))[0] = tmp; //((uint8_t)readsPack[j * READ_CODE_LEN + READ_CODE_LEN-1] >> (8 - readPos%4 * 2)) & 0X03 ;
                    S_Prefix_index++;
                    /*
                    if(S_Prefix_index++ > PARTION_COUNT){ // TODO: realloc
                        string * tmp = realloc(S_Prefix, sizeof(*S_Prefix))
                        fprintf(stderr, "[ERROR] Number of %4d-suffixes overflows\n", prefix);
                        exit(0);
                    }
                    */
                }

            }
        }
        

        for(readPos = READ_LEN - 4; readPos >= 0; readPos--){
            uint8_t byteIndex = (uint8_t)(readPos/4);
            uint8_t baseIndex = readPos % 4;
            if(!(baseIndex)) {
                for(j = 0; j < READ_NUM; j++){
                    if(readsPack[j * READ_CODE_LEN + byteIndex] == prefix) {
                        GetFragment(S_Prefix[S_Prefix_index], readsPack + j * READ_CODE_LEN, readPos);
                        
                        //S_Prefix[READ_CODE_LEN*S_Prefix_index] = (readPos == 0 ? 4 : readsPack[j * READ_CODE_LEN + byteIndex - 1] & 0x03);
                        uint8_t tmp = (readPos == 0 ? 4 : readsPack[j * READ_CODE_LEN + byteIndex - 1] & 0x03);
                        (*(S_Prefix + S_Prefix_index))[0] = tmp;//(readPos == 0 ? 4 : readsPack[j * READ_CODE_LEN + byteIndex - 1] & 0x03);
                        S_Prefix_index++;
                        /*
                        if(S_Prefix_index++ > PARTION_COUNT){
                            fprintf(stderr, "[ERROR] Number of %4d-suffixes overflow\n", prefix);
                        }
                        */
                    }
                } 
            } else {
                for(j = 0; j < READ_NUM; j++){
                    if((uint8_t)(readsPack[j * READ_CODE_LEN + byteIndex] << (baseIndex * 2)) + (uint8_t)(readsPack[j * READ_CODE_LEN + \
                                                                        byteIndex + 1] >> (8 - baseIndex * 2)) == prefix){
                        //GetFragment(S_Prefix + (READ_CODE_LEN)*S_Prefix_index, readsPack + j * READ_CODE_LEN, readPos);
                        GetFragment(S_Prefix[S_Prefix_index], readsPack + j * READ_CODE_LEN, readPos);

                        uint8_t tmp = ((uint8_t)readsPack[j * READ_CODE_LEN + byteIndex] >> (8 - baseIndex * 2)) & 0X03 ;
                        (*(S_Prefix + S_Prefix_index))[0] = tmp; //((uint8_t)readsPack[j * READ_CODE_LEN + byteIndex] >> (8 - baseIndex * 2)) & 0X03 ;
                        S_Prefix_index++;
                        /*
                        if(S_Prefix_index++ > PARTION_COUNT){
                            fprintf(stderr, "[ERROR] Number of %4d-suffixes overflow\n", prefix);
                        }
                        */                        
                    }
                }
            }  
        }

        timestamp = getElapsedTime(startTime);

        partitionTime += timestamp - lastEventTime;

        //Print number of suffixes with particular prefix
        //fprintf(stderr,"%d\t%d\n", prefix, S_Prefix_index); 

        ///////////////////////////////////////////////////////////////////
        //  Print all suffixes, in order of 0000 to 3333
        ///////////////////////////////////////////////////////////////////
        /*
        printf("Suffixes with prefix %d\n", prefix);
        for(i = 0; i < S_Prefix_index; i++){
            for(j = 0; j < READ_CODE_LEN; j++)
                printf("%4X", *(*(S_Prefix + i) + j));
            printf("\n");
        }
        */
        
        /*
        string * aa = (string*)malloc(sizeof(string)*S_Prefix_index);
        for(i = 0; i < S_Prefix_index; i++)
            aa[i] = S_Prefix + READ_CODE_LEN * i;
        */
        
        lastEventTime = getElapsedTime(startTime);

        ta = malloc(S_Prefix_index*sizeof(*ta)); // Total size: 4n
        if(ta == NULL) {
            fprintf(stderr, "%s\n", "Fail to alloc memory: ta");
            exit(0);
        }

        stack = malloc(sizeof(*stack)*READ_NUM);
        if(stack == NULL) {
            fprintf(stderr, "%s\n", "Fail to alloc memory: stack");
            exit(0);
        }

        rsort(S_Prefix, ta, stack, S_Prefix_index);
        free(stack); stack = NULL;
        free(ta); ta = NULL;
        timestamp = getElapsedTime(startTime);

        sortTime += timestamp - lastEventTime;       
        
        ///////////////////////////////////////////////////////////////////
        //  Print sorted suffixes
        ///////////////////////////////////////////////////////////////////
        
        printf("Sorted suffixes with prefix %d\n", prefix);
        for(i = 0; i < S_Prefix_index; i++){
            for(j = 0; j < READ_CODE_LEN; j++)
                printf("%4X", *(*(S_Prefix + i) + j));
            printf("\n");
        }
        

        ///////////////////////////////////////////////////////////////////
        //  Print BWT
        ///////////////////////////////////////////////////////////////////
        /*
        uint8_t * bwt = (uint8_t*)malloc(sizeof(uint8_t)*S_Prefix_index);
        for(i = 0; i < S_Prefix_index; i++){
            bwt[i] = *(aa+i)[0];
        }
        
        fwrite(bwt,sizeof(uint8_t), S_Prefix_index, fout);
        free(bwt);
        */
        free(S_Prefix); S_Prefix = NULL;
        free(prefixMemory); prefixMemory = NULL;
        //fprintf(stderr, "Prefix %4d Done\n", prefix);//omp_get_thread_num());
    }
    //fclose(fout);
    fprintf(stderr, "(Elapsed time (Suffixe Partition) : %9.4f seconds)\n", partitionTime);
    fprintf(stderr, "(Elapsed time (Suffixes Sorting ) : %9.4f seconds)\n", sortTime);

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

//Convert 4 bases to an uint8_t
uint8_t Convert4BaseToOneUint8(char * readsPack, int pos){
    if(pos+3 > strlen(readsPack)){
        fprintf(stderr, "%s\n", "Convert4BaseToOneUint8 Error");
        exit(0);
    }
    uint8_t i1 = (readsPack[pos]   == 'A' ? 0 : (readsPack[pos]   == 'C' ? 1 : (readsPack[pos]   == 'G'? 2 : 3)));
    uint8_t i2 = (readsPack[pos+1] == 'A' ? 0 : (readsPack[pos+1] == 'C' ? 1 : (readsPack[pos+1] == 'G'? 2 : 3)));
    uint8_t i3 = (readsPack[pos+2] == 'A' ? 0 : (readsPack[pos+2] == 'C' ? 1 : (readsPack[pos+2] == 'G'? 2 : 3)));
    uint8_t i4 = (readsPack[pos+3] == 'A' ? 0 : (readsPack[pos+3] == 'C' ? 1 : (readsPack[pos+3] == 'G'? 2 : 3)));
    return i1 * 64 + i2 * 16 + i3 * 4 + i4;
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

void rsort(string * a, string * ta, List * stack, unsigned long long n){ //Sort n strings
    unsigned int pushTimes = 0;
    //List *stack = malloc(sizeof(*stack)*READ_NUM);// It's somewhat dangerous. To be fixed.
    List *sp = stack; //Fetch the top element of stack
    string                   *pile[256], *ai, *ak;
    static unsigned int      count[256] = {0};
    unsigned int             b = 1, c, cmin, *cp, nc = 0; // nc: number of unempty buckets

#define push(a, n, i)   sp->sa = a, sp->sn = n, (sp++)->si = i
#define  pop(a, n, i)   a = (--sp)->sa, n = sp->sn, i = sp->si



    push(a, n, b); // sp->sa = a, sp->sn = n, sp->si = 0, sp++; 
    pushTimes++;

    while(sp > stack) { //stack not empty
        pop(a, n, b); // sp--, a = sp->sa, n = sp->sn, b = sp->si;
       
        /*
        ak = a + n - 1;
        uint8_t fisrtBWTCode = (*ak)[0];
        for(ak = a+n; --ak >= a; ) {
            c = (*ak)[0]; //Fetch the b-th byte of ak-th string
            if(c != fisrtBWTCode) 
                break;
        }
        
        //if(ak < a) {fprintf(stderr, "\tSame BWT Code.\n" ); continue;}
        if(ak < a) 
            continue;
        */
        //When the total num is less thah THRESHOLD, sort them by means of comparison.
        //if( n == 0) fprintf(stderr, "\tPop error: n == 0\n");
        
        if(n < THRESHOLD) {
            SelectionSort(a, n, b);
            //fprintf(stderr, "\tSelection sort invoked for %llu suffixes.\n", n);
            continue;
        }
        
        cmin = 255; //Minimum value in the bytes
        for(ak = a+n; --ak >= a; ) {
            c = (*ak)[b]; //Fetch the b-th byte of ak-th string
            if(++count[c] == 1 && c > 0) { // count[c]++; then check if c is a new value. //Whatif c == 0 ?
                if(c < cmin) cmin = c; //Record the minimum value
                nc++; //New bucket
            }
        }
        if(count[0] >1) push(a, count[0], b+1);
        ak = a + count[0];         
        pile[0]  = ak;         
        count[0] = 0;
        for(cp = count+cmin; nc > 0; cp++, nc--) {
            while(*cp == 0) cp++; //find next cp s.t. cp != 0;
            if(*cp > 1){
                push(ak, *cp, b+1); // Sort ak by the (b+1)-th byte, with total num *cp
                pushTimes++;
                //fprintf(stderr, "\tNew Push.\n");
            }
            ak += *cp;
            pile[cp - count] = ak;
            *cp = 0;
        }

        for(ak = ta+n, ai = a+n; ak > ta; )
            *--ak = *--ai;
        for(ak = ta+n; ak-- > ta; )
            *--pile[(*ak)[b]] = *ak;
    }
    //free(stack); stack = NULL;
    //fprintf(stderr, "Push Times: %u\n", pushTimes);
}

void SelectionSort(string *a, unsigned long long n, int b) {
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