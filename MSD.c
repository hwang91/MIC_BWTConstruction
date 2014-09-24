#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <inttypes.h>

#define SIZE        2
#define THRESHOLD   15

typedef uint8_t * string;
typedef struct List{
    string * sa;
    int sn;
    int si;
} List;

void rsort(string *, int);
void simpleSort(string *, int, int);

int main(){
/*
    uint8_t s0[10] = {1,2,3,4,5,6,7,8,9,2};
    uint8_t s1[10] = {3,2,3,4,5,6,7,8,9,0};
    uint8_t s2[10] = {6,2,3,4,5,6,7,8,9,0};
    uint8_t s3[10] = {5,2,3,4,5,6,7,8,9,0};
    uint8_t s4[10] = {1,2,3,4,5,6,7,8,9,0};
    uint8_t s5[10] = {1,2,3,4,5,6,7,8,9,0};

    string * aa = (string*)malloc(sizeof(string)*6);
    
    aa[0] = s0;
    aa[1] = s1;
    aa[2] = s2;
    aa[3] = s3;
    aa[4] = s4;
    aa[5] = s5;
*/
    uint8_t stest[6][10] = {
        126,2,3,4,5,6,7,8,9,2,
        32,2,3,4,5,6,7,8,9,0,
        63,2,3,4,5,6,7,8,9,0,
        61,2,3,4,5,6,7,8,9,0,
        18,2,3,4,5,6,7,8,9,0,
        18,2,3,4,5,6,7,8,9,0
    };

    string * aa = (string*)malloc(sizeof(string)*6);

    int i;
    for(i = 0; i < 6; i++)
        aa[i] = stest[i];


    rsort(aa, 6);

    int j;
    for(i = 0; i < 6; i++){
        for(j = 0; j < 10; j++)
            printf("%4u", *(*(aa + i) + j));
        printf("\n");
    }

    free(aa);
    return 0;
}

void simpleSort(string * a, int n, int b){
    ;
}

void rsort(string * a, int n){ //Sort n strings

    List stack[SIZE], *sp = stack; 
    string          *pile[256], *ai, *ak, *ta;
    static int      count[256] = {0};
    int             b = 0, c, cmin, *cp, nc = 0; // nc: number of unempty buckets

#define push(a, n, i)   sp->sa = a, sp->sn = n, (sp++)->si = i
#define pop(a, n, i)    a = (--sp)->sa, n = sp->sn, i = sp->si
#define stackempty()    (sp <= stack)


    ta = malloc(n*sizeof(string)); // Total size: n
    push(a, n, 0); // sp->sa = a, sp->sn = n, sp->si = 0, sp++; 
    while(sp > stack) {
        pop(a, n, b); // sp--, a = sp->sa, n = sp->sn, b = sp->si;
        //Ignore the bytes before b
    /*
        //When the total num is less thah THRESHOLD, sort them by means of comparison.
        if(n < THRESHOLD) {
            simpleSort(a, n, b);
            continue;
        }
    */
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
    free(ta);
}