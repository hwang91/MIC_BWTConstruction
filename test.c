#include <stdio.h>

char * Rotate(char *s, int i);

int main(){

	char s[5][5] = 
			{	
				"31025",
				"40123",
				"42150",
				"32302",
				"12433"
			};

	int rank_before[5], rank_after[5];
	int i;
	for(i = 0; i < 5; i++){
		rank_after[i] = i;
	}

	for(i = 4; i >= 0; i--){

		int C[6], k, value, pos;
		// set rank_before to rank_after
		for(k = 0; k < 5; k++)
			rank_before[k] = rank_after[k];

		for(k =0; k < 6; k++)
			C[k] = 0;

		for(k = 0; k < 5; k++)
			C[s[k][i] - '0']++;

	    for(k = 1; k < 6; k++)
            C[k] = C[k] + C[k-1];
        
        for(k = 4; k>= 0; k--){
        	value = s[rank_before[k]][i] - '0';
        	pos = C[value];
        	rank_after[pos - 1] = rank_before[k];
        	C[value]--;
        }

       	for(k = 0; k < 5; k++)
			printf("%d\t", rank_after[k]);
		printf("\n");


	}


	return 0;
}
