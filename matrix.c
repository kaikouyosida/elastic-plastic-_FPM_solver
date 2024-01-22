#include<stdio.h>
#include<stdlib.h>

double** matrix(int m, int n){
	double **A;
	if((A=(double **)calloc(m,sizeof(double *))) == NULL){
		printf("Error:Matrix Memory is not enough\n");
		exit(-1);
	}
	if((A[0]=(double *)calloc(m*n,sizeof(double))) == NULL){
		printf("Error:Matrix Memory is not enough\n");
		exit(-1);
	}
	for(int i=1;i<m;i++) A[i] = A[i-1] + n;

	for(int i=0;i<m;i++){
		for(int j=0;j<n;j++) A[i][j] = 0.;
	}
	return A;
}
void free_matrix(double **A){
	free( (double *) A[0] );
	free( (double **) A );
}
double ***threetimes_tensor(int m, int n, int l){
    double ***A;
    if((A = (double ***)calloc(m, sizeof(double **))) == NULL){
        printf("Error:Tensor Memory is not enough\n");
        exit(-1);
    }
    if((A[0] = (double **)calloc(m * n, sizeof(double *))) == NULL){
        printf("Error:Tensor Memory is not enough\n");
        exit(-1);
    }
    if((A[0][0] = (double *)calloc(m * n * l, sizeof(double))) == NULL){
        printf("Error:TensorMemory is not enough\n");
        exit(-1);
    }

    for(int i = 1; i < m; i++) A[i] = A[i-1] + n;
	
    for(int i = 0; i < m; i++){
        for(int j = 0; j < n; j++) A[i][j] = A[0][0] + i * n * l + j * l;
    }
    for(int i = 0; i < m; i++){
        for(int j = 0; j < n; j++){
            for(int k = 0; k < l; k++) A[i][j][k] = 0.0;
        }
    }
    return A;
}
void free_tensor(double ***A){
    free((double *) A[0][0]);
    free((double **) A[0]);
    free((double ***) A);
}