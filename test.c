#include <stdio.h>

int main(){
	char s[8] = "21300222";
	int count[4] = {0};
	int i;

	for (i = 0; i < 8; i++)
		count[s[i] - '0'] ++;

	for (i = 0; i < 4; i++)
		printf("%d\n",count[i]);
}