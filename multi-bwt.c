//bwt.c
// Wang Heng
// Sept 3, 2014

#include <stdio.h>

#define LEN 7
#define ALPHABETA_SIZE 3

int main(){

    //ACGT$, TAGT$, GGAA$}

    char S[LEN] = "abaaba$";
    char s[LEN] = "1211210";
  

    int rank_before[LEN], rank_after[LEN];
    int i;
    for(i = 0; i < LEN; i++){
        rank_after[i] = i;
    }

    //LSD Radix Sort
    for(i = LEN - 1; i >= 0; i--){
        //Count sort
        int C[ALPHABETA_SIZE], k, value, pos;
        
        // set rank_before to rank_after
        // TODO: change to verctor assignment.
        for(k = 0; k < LEN; k++) 
            rank_before[k] = rank_after[k];

        //Gnerally the ALPHABETA_SIZE is very small
        for(k =0; k < ALPHABETA_SIZE; k++)
            C[k] = 0;

        //For multiple reads, the LEN should be changed to NUM_OF_READS
        // Try to come up with a good parallel way
        // An easy way is to seperate the reads to severeal parts and count respectively
        for(k = 0; k < LEN; k++)
            C[s[(k+i)%(LEN)] - '0']++;

        //Gnerally the ALPHABETA_SIZE is very small
        for(k = 1; k < ALPHABETA_SIZE; k++)
            C[k] = C[k] + C[k-1];
        
        for(k = LEN - 1; k>= 0; k--){
            value = s[(rank_before[k] + i) % LEN] - '0';
            pos   = C[value];
            rank_after[pos - 1] = rank_before[k];
            C[value]--;
        }

        for(k = 0; k < LEN; k++)
            printf("%d\t", rank_after[k]);
        printf("\n");

    }

    for(i = 0; i < LEN; i++)
        printf("%c", S[ (rank_after[i] + LEN - 1) % LEN ]);
    printf("\n");


    return 0;
}
