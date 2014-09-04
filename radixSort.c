//radixSort.c
// Wang Heng
// Sept 3, 2014

#include <stdio.h>

int main(){
    char s[15][5] = 
    { 
        "12340",
        "23400",
        "34000",
        "40000",
        "00000",
        "41340",
        "13400",
        "34000",
        "40000",
        "00000",
        "33110",
        "31100",
        "11000",
        "10000",
        "00000"
    };

    int rank_before[15], rank_after[15];
    int i;
    for(i = 0; i < 15; i++){
        rank_after[i] = i;
    }

    //LSD Radix Sort
    for(i = 4; i >= 0; i--){
        //Count sort
        int C[5], k, value, pos;
        // set rank_before to rank_after
        for(k = 0; k < 15; k++)
            rank_before[k] = rank_after[k];

        for(k =0; k < 5; k++)
            C[k] = 0;

        for(k = 0; k < 15; k++)
            C[s[k][i] - '0']++;

        /* Print counts 
        for(k = 0; k < 5; k++)
            printf("%d\t", C[k]);

        putchar('\n');
        */

        for(k = 1; k < 5; k++)
            C[k] = C[k] + C[k-1];
        
        for(k = 14; k>= 0; k--){
            value = s[rank_before[k]][i] - '0';
            pos   = C[value];
            rank_after[pos - 1] = rank_before[k];
            C[value]--;
        }

        for(k = 0; k < 15; k++)
            printf("%d\t", rank_after[k]);

        printf("\n");

    }

    return 0;
}
