#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>

#define READ_CODE_LEN   25
#define THRESHOLD       25

typedef uint8_t * string;

typedef struct List{
    string * sa;
    int sn;
    int si;
} List;


unsigned long long partitionCount[256];
string * S_Prefix = 0;
string prefixMemory = 0;
string * ta = 0;
List * stack = 0;


void PrinfSuffixes(string * s, unsigned long long num);
void rsort(string *, string *, List*, unsigned long long);
void SelectionSort(string *a, unsigned long long n, int b);


int main(int argc, char** argv){

    uint8_t prefix = atoi(argv[argc-1]);

    /////////////////////////////////////////////////////////////
    //  PartitionCount
    /////////////////////////////////////////////////////////////
    FILE *partition = fopen(argv[2], "rb");
    if(partition == NULL) {
        fprintf(stderr, "%s\n", "Fail to open file");
        exit(0);
    }
    fread(partitionCount, sizeof(unsigned long long), 256, partition);
    fclose(partition);

    /////////////////////////////////////////////////////////////
    //  Load Partition
    /////////////////////////////////////////////////////////////
    FILE *fp = fopen(argv[1], "rb");
    if(fp == NULL) {
        fprintf(stderr, "%s\n", "Fail to open file");
        exit(0);
    }

    S_Prefix = malloc(sizeof(*S_Prefix) * partitionCount[prefix]);
    if(S_Prefix == NULL) {
        fprintf(stderr, "%s\n", "Fail to alloc memory: S_Prefix");
        exit(0);
    }
    
    prefixMemory = malloc(sizeof(*prefixMemory)*READ_CODE_LEN* partitionCount[prefix]);

    unsigned int readNum = 0;
    while(fread(prefixMemory + readNum*READ_CODE_LEN*1024, sizeof(uint8_t), READ_CODE_LEN*1024, fp))
        readNum++;
    fclose(fp);



    unsigned long long i, j;
    for(i = 0; i < partitionCount[prefix]; i++){
        S_Prefix[i] = prefixMemory + i * READ_CODE_LEN;
        if(S_Prefix[i] == NULL){
            fprintf(stderr, "Fail to alloc memory: S_Prefix[%llu]\n", i);
        }
    }

        
   ///////////////////////////////////////////////////////////////////
   //  Print Pre-sort suffixes
   ///////////////////////////////////////////////////////////////////        
   //printf("Pre-sort suffixes with prefix %d\n", prefix);
   //PrinfSuffixes(S_Prefix, partitionCount[prefix]);


   ///////////////////////////////////////////////////////////////////
   //  Sort suffixes
   ///////////////////////////////////////////////////////////////////           
    ta = malloc(sizeof(*ta)*partitionCount[prefix]); // Total size: 4n
    if(ta == NULL) {
        fprintf(stderr, "%s\n", "Fail to alloc memory: ta");
        exit(0);
    }
    stack = malloc(sizeof(*stack)*partitionCount[prefix]);
    if(stack == NULL) {
        fprintf(stderr, "%s\n", "Fail to alloc memory: stack");
        exit(0);
    }
    rsort(S_Prefix, ta, stack, partitionCount[prefix]);
    free(stack); stack = NULL;
    free(ta); ta = NULL;
   

    ///////////////////////////////////////////////////////////////////
    //  Print sorted suffixes
    ///////////////////////////////////////////////////////////////////        
    //printf("Sorted suffixes with prefix %d\n", prefix);
    //PrinfSuffixes(S_Prefix, partitionCount[prefix]);


    free(prefixMemory);
    free(S_Prefix);

}

void PrinfSuffixes(string * s, unsigned long long num){
    unsigned long long i, j;
   for(i = 0; i < num; i++){
       for(j = 0; j < READ_CODE_LEN; j++)
           printf("%4X", *(*(s + i) + j));
       printf("\n");
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