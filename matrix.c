#include <stdio.h>
#include <stdlib.h>

int main(){
    double A[9], *B;
    for(int i = 0; i<9; i++){
        A[i]=i;
    }
    B = malloc(9*sizeof(double));
    for(int j=0; j<9; j++){
        B[j]=A[j];
    }
    
    for(int j=0; j<9; j++){
        printf("\n %f", B[j]);
    }



}