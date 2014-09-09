//bwt.c
// Wang Heng
// Sept 3, 2014

#include <stdio.h>

#define LEN 5
#define ALPHABETA_SIZE 3

int main(){

    //ACGT$, TAGT$, GGAA$}

    //char S1[LEN] = "ACGT$";
    char s1[LEN] = "12340";
    char s2[LEN] = "41340";
    char s3[LEN] = "33110";
  

    int rank_before[LEN * 3], rank_after[LEN * 3];
    int i;
    for(i = 0; i < LEN * 3; i++){
        rank_after[i] = i;
    }

    //LSD Radix Sort
    for(i = LEN - 1; i >= 0; i--){
        //Count sort
        int C[ALPHABETA_SIZE], k, value, pos;
        
        // set rank_before to rank_after
        // TODO: change to verctor assignment.
        for(k = 0; k < LEN * 3 ; k++) 
            rank_before[k] = rank_after[k];

        //Gnerally the ALPHABETA_SIZE is very small
        for(k =0; k < ALPHABETA_SIZE; k++)
            C[k] = 0;

        //For multiple reads, the LEN should be changed to NUM_OF_READS
        // Try to come up with a good parallel way
        // An easy way is to seperate the reads to severeal parts and count respectively
        for(k = 0; k < LEN; k++){
            C[s1[(k+i)%(LEN)] - '0']++;
            C[s2[(k+i)%(LEN)] - '0']++;
            C[s3[(k+i)%(LEN)] - '0']++;
        }

        //Gnerally the ALPHABETA_SIZE is very small
        for(k = 1; k < ALPHABETA_SIZE; k++)
            C[k] = C[k] + C[k-1];
        
        for(k = LEN - 1; k>= 0; k--){
            value = s[(rank_before[k] + i) % LEN] - '0';
            pos   = C[value];
            rank_after[pos - 1] =÷ rank_before[k];
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
