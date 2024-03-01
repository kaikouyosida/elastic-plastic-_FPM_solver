#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"type.h"
#include"scalar.h"

extern Option option;

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

    for(int i = 0; i < m; i++) A[i] = A[0] + n * i;
	
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
void calculate3x3MatrixSquare(double matrix_out[3][3],
                              double matrix_in[3][3])
{
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
//法線ベクトルの計算
void calc_Ne(int dim, int subdomain_n1, int subdomain_n2, int face, int *vertex_offset, int *node, double *node_xyz, double *center_xyz, double (*N_matrix)[6])
{
	double unit_n_vec[3]; // 単位法線ベクトル
	double direction[3];  // サブドメインE1からE2の方へ向くベクトル
	double norm_vec = 0.; // ベクトルの大きさ
	int ref_offset = 0;	  // 面を分割するときの基準オフセット(3次元のときのみ使用)
	double edge1[3];
	double edge2[3]; // 面を分割してできた三角形の各辺のベクトル(3次元のみ)
	double N_e[3]; //法線ベクトル

	for (int i = 0; i < dim; i++)
		direction[i] = center_xyz[dim * subdomain_n2 + i] - center_xyz[dim * subdomain_n1 + i];

	if (dim == 2)
	{	
		unit_n_vec[0] = node_xyz[dim * node[vertex_offset[face]] + 1] - node_xyz[dim * node[vertex_offset[face] + 1] + 1];
		unit_n_vec[1] = node_xyz[dim * node[vertex_offset[face] + 1]] - node_xyz[dim * node[vertex_offset[face]]];
		norm_vec = norm(unit_n_vec, dim);
		for (int i = 0; i < dim; i++)
			unit_n_vec[i] /= norm_vec;

		if (dot_product(dim, unit_n_vec, direction) < 0.)
		{
			for (int i = 0; i < dim; i++)
				unit_n_vec[i] *= -1.0;
		}

		for (int i = 0; i < dim; i++)
			N_e[i] = unit_n_vec[i];
	}
	else if (dim == 3)
	{
		for (int k = 0; k < dim; k++)
			N_e[k] = 0.;

		ref_offset = vertex_offset[face];
		for (int i = ref_offset + 2; i < vertex_offset[face + 1]; i++)
		{
			for (int k = 0; k < dim; k++)
				edge1[k] = node_xyz[dim * node[i - 1] + k] - node_xyz[dim * node[ref_offset] + k];
			for (int k = 0; k < dim; k++)
				edge2[k] = node_xyz[dim * node[i] + k] - node_xyz[dim * node[ref_offset] + k];

			cross_product(dim, edge1, edge2, unit_n_vec);
			norm_vec = norm(unit_n_vec, dim);
			for (int k = 0; k < dim; k++)
				unit_n_vec[k] /= norm_vec;

			if (dot_product(dim, unit_n_vec, direction) < 0.)
			{
				for (int k = 0; k < dim; k++)
					unit_n_vec[k] *= -1.0;
			}

			for (int k = 0; k < dim; k++)
				N_e[k] += unit_n_vec[k];
		}
		for (int i = 0; i < dim; i++)
			N_e[i] /= (double)(vertex_offset[face + 1] - ref_offset - 2);

			N_matrix[0][0] = N_e[0];	N_matrix[0][1] = 0.0;		N_matrix[0][2] = 0.0;		N_matrix[0][3] = N_e[1];		N_matrix[0][4] = 0.0;		N_matrix[0][5] = N_e[2];
			N_matrix[1][0] = 0.0;       N_matrix[1][1] = N_e[1];	N_matrix[1][2] = 0.0;		N_matrix[1][3] = N_e[0];		N_matrix[1][4] = N_e[2]; 	N_matrix[1][5] = 0.0;
			N_matrix[2][0] = 0.0;		N_matrix[2][1] = 0.0;		N_matrix[2][2] = N_e[2];	N_matrix[2][3] = 0.0   ;		N_matrix[2][4] = N_e[1]; 	N_matrix[2][5] = N_e[0];	
	}
}
void calc_Ne_diagonal(int dim , int subdomain_n1, int subdomain_n2, int face, int *vertex_offset, int *node, double *node_xyz, double *center_xyz, double (*Ne_d)[9]){
	double Ne[3][6];
	calc_Ne(option.dim, subdomain_n1, subdomain_n2, face, vertex_offset, node, node_xyz, center_xyz, Ne);

	Ne_d[0][0] = Ne[0][0];	Ne_d[0][1] = 0.0;		Ne_d[0][2] = 0.0;		Ne_d[0][3] = Ne[1][1];		Ne_d[0][4] = 0.0;		Ne_d[0][5] = 0.0;      Ne_d[0][6] = Ne[2][2];		Ne_d[0][7] = 0.0;		Ne_d[0][8] = 0.0; 
	Ne_d[1][0] = 0.0;       Ne_d[1][1] = Ne[0][0];	Ne_d[1][2] = 0.0;		Ne_d[1][3] = 0.0   ;		Ne_d[1][4] = Ne[1][1]; 	Ne_d[1][5] = 0.0;      Ne_d[1][6] = 0.0;			Ne_d[1][7] = Ne[2][2]; 	Ne_d[1][8] = 0.0;
	Ne_d[2][0] = 0.0;		Ne_d[2][1] = 0.0;		Ne_d[2][2] = Ne[0][0];	Ne_d[2][3] = 0.0   ;		Ne_d[2][4] = 0.0;	 	Ne_d[2][5] = Ne[1][1]; Ne_d[2][6] = 0.0;	 		Ne_d[2][7] = 0.0;		Ne_d[2][8] = Ne[2][2];	
}