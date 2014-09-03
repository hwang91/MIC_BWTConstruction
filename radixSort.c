// radixSort.c
// Wang Heng
// Sept 3, 2014
// Sort a single int array
// No arguments


#include <stdio.h>

void COUNTINGSORT(int *A, int *R, int array_size, int k);
        
int main()
{
        
   int A[8] = {2, 5, 3, 0, 2, 3, 0, 3}, R[8], i;
   for (i=0; i<= 7; i++){
      printf("%d ", A[i]);
   }
   printf("\n");

   COUNTINGSORT(A, R, 8, 5);
   for (i=0; i<= 7; i++)
   {
      printf("%d ", R[i]);
   }
   printf("\n");

    char test[4] = "2013";
   for (i = 0; i < 4; i++)
      printf("%d\n",test[i] - 48);
    printf("\n");   
   
   return 0;
}


void COUNTINGSORT(int *A, int *R, int array_size, int k)
{
        int C[k+1], i, value, pos;
        for(i=0; i<=k; i++)
        {
            C[i] = 0;
        }
        for(i=0; i< array_size; i++)
        {
            C[A[i]] ++;
        }
        for(i=1; i<=k; i++)
        {
            C[i] = C[i] + C[i-1];
        }

        for(i=array_size-1; i>=0; i--)
        {
            value = A[i];
            pos = C[value];
            //  B[pos-1] = value;
            R[pos-1] = i;
            C[value]--;

        }
}
