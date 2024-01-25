#include<stdio.h>
#include<stdlib.h>
#include<math.h>

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
void identify3x3Matrix(double (*matrix)[3]){
    matrix[0][0] = 1.0;
    matrix[0][1] = 0.0;
    matrix[0][2] = 0.0;
    matrix[1][0] = 0.0;
    matrix[1][1] = 1.0;
    matrix[1][2] = 0.0;
    matrix[2][0] = 0.0;
    matrix[2][1] = 0.0;
    matrix[2][2] = 1.0;
}
void zeroize3x3Matrix(double (*matrix)[3])
{
    matrix[0][0] = 0.0;
    matrix[0][1] = 0.0;
    matrix[0][2] = 0.0;
    matrix[1][0] = 0.0;
    matrix[1][1] = 0.0;
    matrix[1][2] = 0.0;
    matrix[2][0] = 0.0;
    matrix[2][1] = 0.0;
    matrix[2][2] = 0.0;
}
// 正則行列の逆行列をガウスの消去法で計算 //
void inverse_mat3x3(int dim, double (*A)[3], double (*inverse_A)[3]){
	double **B;       // Aの成分をコピーして、計算に用いる行列
  int max_j = 0;    // ピボット選択で用いる補助変数
	double tmp = 0.;  // ピボット選択で用いる補助変数
	B = matrix(dim,dim);
	for(int i=0;i<dim;i++){
		for(int j=0;j<dim;j++) B[i][j] = A[i][j];
	}

	for(int i=0;i<dim;i++){
		for(int j=0;j<dim;j++){
			if(j==i) inverse_A[i][j] = 1.0;
			else inverse_A[i][j] = 0.0;
		}
	}

  // 部分ピボット選択 //
	for(int l=0;l<dim;l++){
		max_j = l;
		for(int i=l;i<dim;i++){
			if(fabs(B[l][l]) < fabs(B[i][l])){
				max_j = i;
			}
		}
		if(max_j != l){
			for(int j=l;j<dim;j++){
				tmp = B[l][j];
				B[l][j] = B[max_j][j];
				B[max_j][j] = tmp;

				tmp = inverse_A[l][j];
				inverse_A[l][j] = inverse_A[max_j][j];
				inverse_A[max_j][j] = tmp;
			}
		}

		for(int j=0;j<l;j++) inverse_A[l][j] /= B[l][l];
		for(int j=dim-1;j>l-1;j--){
			inverse_A[l][j] /= B[l][l];
			B[l][j] /= B[l][l];
		}
		for(int i=l+1;i<dim;i++){
			for(int j=0;j<l+1;j++) inverse_A[i][j] -= B[i][l]*inverse_A[l][j];
			for(int j=l+1;j<dim;j++){
				B[i][j] -= B[i][l]*B[l][j];
				inverse_A[i][j] -= B[i][l]*inverse_A[l][j];
			}
		}
	}

	for(int l=dim-1;l>-1;l--){
		for(int i=l-1;i>-1;i--){
			for(int j=0;j<dim;j++) inverse_A[i][j] -= B[i][l]*inverse_A[l][j];
		}
	}

	free_matrix(B);
}
// Am行AnBm列の行列とAnBm行Bn列の行列の積をXに保存 //
void multi_mat(int Am, int Bn, double** A, double** B, int AnBm, double** X){
	for(int i=0;i<Am;i++){
		for(int j=0;j<Bn;j++) X[i][j] = 0.;
	}

	for(int i=0;i<Am;i++){
		for(int j=0;j<Bn;j++){
			for(int k=0;k<AnBm;k++) X[i][j] += A[i][k]*B[k][j];
		}
	}
}
// m行n列の行列Xの転置行列をXTに保存 //
void trans_mat(int m, int n, double** X, double** XT){
  for(int i=0;i<n;i++){
		for(int j=0;j<m;j++) XT[i][j] = X[j][i];
	}
}
// 正則行列の逆行列をガウスの消去法で計算 //
void inverse_mat(int mn, double** A, double** inverse_A){
	double **B;       // Aの成分をコピーして、計算に用いる行列
  int max_j = 0;    // ピボット選択で用いる補助変数
	double tmp = 0.;  // ピボット選択で用いる補助変数
	B = matrix(mn,mn);
	for(int i=0;i<mn;i++){
		for(int j=0;j<mn;j++) B[i][j] = A[i][j];
	}

	for(int i=0;i<mn;i++){
		for(int j=0;j<mn;j++){
			if(j==i) inverse_A[i][j] = 1.0;
			else inverse_A[i][j] = 0.0;
		}
	}

  // 部分ピボット選択 //
	for(int l=0;l<mn;l++){
		max_j = l;
		for(int i=l;i<mn;i++){
			if(fabs(B[l][l]) < fabs(B[i][l])){
				max_j = i;
			}
		}
		if(max_j != l){
			for(int j=l;j<mn;j++){
				tmp = B[l][j];
				B[l][j] = B[max_j][j];
				B[max_j][j] = tmp;

				tmp = inverse_A[l][j];
				inverse_A[l][j] = inverse_A[max_j][j];
				inverse_A[max_j][j] = tmp;
			}
		}

		for(int j=0;j<l;j++) inverse_A[l][j] /= B[l][l];
		for(int j=mn-1;j>l-1;j--){
			inverse_A[l][j] /= B[l][l];
			B[l][j] /= B[l][l];
		}
		for(int i=l+1;i<mn;i++){
			for(int j=0;j<l+1;j++) inverse_A[i][j] -= B[i][l]*inverse_A[l][j];
			for(int j=l+1;j<mn;j++){
				B[i][j] -= B[i][l]*B[l][j];
				inverse_A[i][j] -= B[i][l]*inverse_A[l][j];
			}
		}
	}

	for(int l=mn-1;l>-1;l--){
		for(int i=l-1;i>-1;i--){
			for(int j=0;j<mn;j++) inverse_A[i][j] -= B[i][l]*inverse_A[l][j];
		}
	}

	free_matrix(B);
}
double calc_3x3matrix_norm(double (*matrix)[3]){
    double temp = 0.0;
    int i, j;

    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            temp += matrix[i][j] * matrix[i][j];

    return sqrt(temp);
}
void calc_3x3_matrix_square(double(*matrix_out)[3], double (*matrix_in)[3]){
    int i, j, k;

    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
        {
            double temp = 0.0;

            for (k = 0; k < 3; k++)
                temp
                    += matrix_in[i][k] * matrix_in[k][j];

            matrix_out[i][j] = temp;
        }
}
double calc_3x3matrix_determinant(double (*matrix)[3])
{
    return
        matrix[0][0] * matrix[1][1] * matrix[2][2]
        + matrix[0][1] * matrix[1][2] * matrix[2][0]
        + matrix[0][2] * matrix[1][0] * matrix[2][1]
        - matrix[0][0] * matrix[1][2] * matrix[2][1]
        - matrix[0][1] * matrix[1][0] * matrix[2][2]
        - matrix[0][2] * matrix[1][1] * matrix[2][0];
}