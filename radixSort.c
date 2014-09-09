//radixSort.c
// Wang Heng
// Sept 3, 2014

#include <stdio.h>
#include <stdlib.h>

//#define PRINT_PASS
//#define PRINT_COUNT

#define ALPHABETA_SIZE  4
#define READ_LEN        4
#define READ_NUM        3

char s[READ_NUM][READ_LEN];
int  Count[READ_LEN][ALPHABETA_SIZE + 1] = {0}; 
char ConvertBase2Int(char c); //Convert ACGT to 1234;

int main(int argc, char** argv){

    int i = 0;
    int j = 0;
    FILE *fp;

    if((fp = fopen("./test.fa", "rt")) == NULL){
        printf("Cannot open file strike any key exit!");
        exit(1);
    }

    char c = fgetc(fp);
    while(c != EOF){
        if(c == '\n'){
            i++ ;
            j = 0;
            c = fgetc(fp);
        }
        s[i][j] = ConvertBase2Int(c);
        Count[j][ConvertBase2Int(c) - '0']++;
        j++ ;
        c = fgetc(fp);
    }


    fclose(fp);

    for(i = 0; i < READ_NUM; i++){
        for(j = 0; j < READ_LEN; j++)
            putchar(s[i][j]);
        putchar('\n');
    }

    for(i = READ_LEN - 2; i >= 0; i--){
        for(j = 0; j < ALPHABETA_SIZE + 1; j++)
            Count[i][j] += Count[i + 1][j];
    }

    // Calculate # of 0 
    for(j = 0; j < READ_LEN; j++)
        Count[j][0] = READ_LEN * READ_NUM  - Count[j][1] - Count[j][2] - Count[j][3] - Count[j][4];
    
    //Accumulation
    for(i = 0; i < READ_LEN; i++){
        for(j = 1; j < ALPHABETA_SIZE + 1; j++)
            Count[i][j] += Count[i][j - 1];
    }

#ifdef PRINT_COUNT

    for(i = 0; i < READ_LEN; i++){
        for(j = 0; j < ALPHABETA_SIZE + 1; j++)
            printf("%d\t", Count[i][j]);
        putchar('\n');
    }

#endif

    int rank_before[READ_LEN * READ_NUM], rank_after[READ_LEN * READ_NUM];

    for(i = 0; i < READ_LEN * READ_NUM; i++){
        rank_after[i] = i;
    }

    //LSD Radix Sort
    for(i = READ_LEN - 1; i >= 0; i--){
        //Count sort
        int k, value, pos;
        // set rank_before to rank_after
        for(k = 0; k < READ_LEN * READ_NUM; k++)
            rank_before[k] = rank_after[k];
        

        for(k = READ_NUM * READ_LEN - 1; k >= 0; k--){
            value = ((rank_before[k] % READ_LEN + i) >= READ_LEN )? 0 : s[rank_before[k] / READ_LEN][rank_before[k] % READ_LEN + i] - '0';
            pos   = Count[i][value];
            rank_after[pos - 1] = rank_before[k];
            Count[i][value]--;
        }
        
    #ifdef PRINT_PASS

        for(k = 0; k < READ_LEN * READ_NUM; k++)
            printf("%d\t", rank_after[k]);

        printf("\n");

    #endif

    }

    for(i = 0; i < READ_LEN * READ_NUM; i++){
        if(rank_after[i] % READ_LEN == 0) {
            putchar('0');
            continue;
        }
        printf("%c", s[rank_after[i]/READ_LEN][rank_after[i] % READ_LEN - 1]);
    }
    printf("\n");

    return 0;
}

char ConvertBase2Int(char c){
    char tmp;
    switch(c){
        case 'A': tmp = '1'; break;
        case 'a': tmp = '1'; break;
        case 'C': tmp = '2'; break;
        case 'c': tmp = '2'; break;
        case 'G': tmp = '3'; break;
        case 'g': tmp = '3'; break;
        case 'T': tmp = '4'; break;
        case 't': tmp = '4'; break;
        default : tmp = '0'; 
    }
    return tmp;
}
