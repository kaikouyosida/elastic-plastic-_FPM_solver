#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"type.h"
#include"scalar.h"
#include"vector.h"

#define NUMBER_OF_NODE_IN_FACE 4
#define NUMBER_OF_NODE_IN_SUBDOMAIN 8

extern Global global;
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

// Am行AnBm列の行列とAnBm行Bn列の行列の積をXに保存 //
void multi_mat(int Am, int Bn, double** A, double** B, int AnBm, double** X){
	for(int i = 0; i < Am; i++){
		for(int j = 0; j < Bn; j++)
			X[i][j] = 0.;
	}

	for(int i = 0; i < Am; i++){
		for(int j = 0; j < Bn; j++){
			for(int k = 0; k < AnBm; k++)
				X[i][j] += A[i][k] * B[k][j];
		}
	}
}

// m行n列の行列Xの転置行列をXTに保存 //
void trans_mat(int m, int n, double** X, double** XT){
  for(int i = 0; i < n; i++){
		for(int j = 0; j < m; j++)
			XT[i][j] = X[j][i];
	}
}

// 正則行列の逆行列をガウスの消去法で計算 //
void inverse_mat(int mn, double** A, double** inverse_A){
	double **aug;
	aug = matrix(mn, 2 * mn);

	for(int i = 0; i < mn; i++){
		for(int j = 0; j < mn; j++){
			aug[i][j] = A[i][j];
			aug[i][j + mn] = (i == j) ? 1 : 0;
		}
	}
	//対角要素を１にする
	for(int i = 0; i < mn; i++){
		double diag = aug[i][i];
        if (diag == 0){
            printf("Inverse matrix is not exist!\n");
            exit(-1);
        }
        for (int j = 0; j < 2 * mn; j++) {
            aug[i][j] /= diag;
        }
		// 他の行の該当列を0にする
        for (int k = 0; k < mn; k++) {
            if (k != i) {
                double factor = aug[k][i];
                for (int j = 0; j < 2 * mn; j++) {
                    aug[k][j] -= factor * aug[i][j];
                }
            }
        }
	}
	// 逆行列を取り出す
    for (int i = 0; i < mn; i++) {
        for (int j = 0; j < mn; j++) {
            inverse_A[i][j] = aug[i][j + mn];
        }
    }
	free_matrix(aug);
}

//３×３マトリクスの二重積
double calc_3x3matrix_norm(double (*matrix)[3]){
    double temp = 0.0;
    int i, j;

    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            temp += matrix[i][j] * matrix[i][j];

    return sqrt(temp);
}

//３×３マトリクスの行列式
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

//３×３マトリクスの積
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
void calc_Ne(int dim, int subdomain_n1, int subdomain_n2, int face, int *vertex_offset, int *node, double *node_xyz, double *center_xyz, double *Ne)
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
				Ne[k] += unit_n_vec[k];
		}
		for (int i = 0; i < dim; i++)
			Ne[i] /= (double)(vertex_offset[face + 1] - ref_offset - 2);
	}
}

void calc_unit_vector(double unit_vector[3], const int face_n, const int subdomain1, const int subdomain2, const double *current_point_XYZ){
	int ref_num = global.subdomain.vertex_offset[face_n];
	int subdomain_node[NUMBER_OF_NODE_IN_SUBDOMAIN];
	int node_id[NUMBER_OF_NODE_IN_FACE];
	double edge1[3];
	double edge2[3]; // 面を分割してできた三角形の各辺のベクトル
	double edge1_cross_edge2[3];
	double direction[3];
	double area_norm;
	int face_node[NUMBER_OF_NODE_IN_FACE];

	//単位ベクトルの方向を決定＋単位法線ベクトルを初期化
	for(int i = 0; i < option.dim; i++){
		direction[i] = current_point_XYZ[option.dim * subdomain2 + i] - current_point_XYZ[option.dim * subdomain1 + i];
		unit_vector[i] = 0.;
	}
	
	//face_nにおけるノード番号のアドレスを取得
	if(option.solver_type == 1){
		generate_subdomain_node(subdomain1, subdomain_node);
		generate_node_id(face_n, subdomain1, subdomain_node, node_id);
	}

	//四辺形を2つの三角形に分割→二つの法線ベクトルの平均をその面の法線ベクトルとして計算
	for(int i = 0; i < 2; i++){

		if(option.solver_type == 1){
			generate_current_edge_vector(edge1, subdomain1, node_id[2 * i + 1], node_id[2 * i], subdomain_node);
			generate_current_edge_vector(edge2, subdomain1, node_id[(2 + i) % 3 + 1], node_id[2 * i], subdomain_node);
		}else{
			for(int j = 0; j < NUMBER_OF_NODE_IN_FACE; j++){
				face_node[j] = global.subdomain.node[global.subdomain.vertex_offset[face_n] + j];
			}
			
			for(int j = 0; j < option.dim; j++){
				edge1[j] = global.subdomain.node_XYZ[option.dim * face_node[2 * i + 1] + j] - global.subdomain.node_XYZ[option.dim * face_node[2 * i] + j];
				edge2[j] = global.subdomain.node_XYZ[option.dim * face_node[(2 + i) % 3 + 1] + j] - global.subdomain.node_XYZ[option.dim * face_node[2 * i] + j];
			}
		}

		cross_product(option.dim, edge1, edge2, edge1_cross_edge2);

		area_norm = norm(edge1_cross_edge2, option.dim);
		for(int j = 0; j < option.dim; j++)
			edge1_cross_edge2[j] /= area_norm;

		if(dot_product(option.dim, direction, edge1_cross_edge2) < 0.)
			for(int j = 0; j < option.dim; j++)
				edge1_cross_edge2[j] *= -1.0;

		for(int j = 0; j < option.dim; j++){
			unit_vector[j] += edge1_cross_edge2[j];
		}
	}
	
	for(int i = 0; i < option.dim; i++)
		unit_vector[i] /= 2.0;
}

void generate_unit_vec_to_mat3x6(const int face_n, const int subdomain1, const int subdomain2, const double *current_point_XYZ, double (*N_matrix)[6]){
	double Ne[3];
	calc_unit_vector(Ne, face_n, subdomain1, subdomain2, current_point_XYZ);
	N_matrix[0][0] = Ne[0];	N_matrix[0][1] = 0.0;		N_matrix[0][2] = 0.0;		N_matrix[0][3] = Ne[1];		N_matrix[0][4] = 0.0;		N_matrix[0][5] = Ne[2];
	N_matrix[1][0] = 0.0;   N_matrix[1][1] = Ne[1];   	N_matrix[1][2] = 0.0;		N_matrix[1][3] = Ne[0];		N_matrix[1][4] = Ne[2]; 	N_matrix[1][5] = 0.0;
	N_matrix[2][0] = 0.0;	N_matrix[2][1] = 0.0;		N_matrix[2][2] = Ne[2];		N_matrix[2][3] = 0.0 ;		N_matrix[2][4] = Ne[1]; 	N_matrix[2][5] = Ne[0];
}

void generate_unit_vec_to_mat3x9(int face_n, int subdomain1, int subdomain2, double *current_point_XYZ, double (*N_matrix)[9]){
	double Ne[3];
	calc_unit_vector(Ne, face_n, subdomain1, subdomain2, current_point_XYZ);
	N_matrix[0][0] = Ne[0]; 	N_matrix[0][1] = 0.0;		N_matrix[0][2] = 0.0;		N_matrix[0][3] = Ne[1];			N_matrix[0][4] = 0.0;		N_matrix[0][5] = 0.0;      N_matrix[0][6] = Ne[2];			N_matrix[0][7] = 0.0;		N_matrix[0][8] = 0.0; 
	N_matrix[1][0] = 0.0;       N_matrix[1][1] = Ne[0];		N_matrix[1][2] = 0.0;		N_matrix[1][3] = 0.0   ;		N_matrix[1][4] = Ne[1]; 	N_matrix[1][5] = 0.0;      N_matrix[1][6] = 0.0;			N_matrix[1][7] = Ne[2]; 	N_matrix[1][8] = 0.0;
	N_matrix[2][0] = 0.0;		N_matrix[2][1] = 0.0;		N_matrix[2][2] = Ne[0];		N_matrix[2][3] = 0.0   ;		N_matrix[2][4] = 0.0;	 	N_matrix[2][5] = Ne[1]; 	N_matrix[2][6] = 0.0;	 		N_matrix[2][7] = 0.0;		N_matrix[2][8] = Ne[2];	
}

void calc_Ne_3x6(int subdomain_n1, int subdomain_n2, int face, int *vertex_offset, int *node, double *node_xyz, double *center_xyz, double (*N_matrix)[6]){
	double Ne[3];
	calc_Ne(option.dim, subdomain_n1, subdomain_n2, face, vertex_offset, node, node_xyz, center_xyz, Ne);

	N_matrix[0][0] = Ne[0];	N_matrix[0][1] = 0.0;		N_matrix[0][2] = 0.0;		N_matrix[0][3] = Ne[1];		N_matrix[0][4] = 0.0;		N_matrix[0][5] = Ne[2];
	N_matrix[1][0] = 0.0;   N_matrix[1][1] = Ne[1];   	N_matrix[1][2] = 0.0;		N_matrix[1][3] = Ne[0];		N_matrix[1][4] = Ne[2]; 	N_matrix[1][5] = 0.0;
	N_matrix[2][0] = 0.0;	N_matrix[2][1] = 0.0;		N_matrix[2][2] = Ne[2];		N_matrix[2][3] = 0.0 ;		N_matrix[2][4] = Ne[1]; 	N_matrix[2][5] = Ne[0];

}

void calc_Ne_3x9(int subdomain_n1, int subdomain_n2, int face, int *vertex_offset, int *node, double *node_xyz, double *center_xyz, double (*N_matrix)[9]){
	double Ne[3];
	calc_Ne(option.dim, subdomain_n1, subdomain_n2, face, vertex_offset, node, node_xyz, center_xyz, Ne);

	N_matrix[0][0] = Ne[0]; 	N_matrix[0][1] = 0.0;		N_matrix[0][2] = 0.0;		N_matrix[0][3] = Ne[1];			N_matrix[0][4] = 0.0;		N_matrix[0][5] = 0.0;      N_matrix[0][6] = Ne[2];			N_matrix[0][7] = 0.0;		N_matrix[0][8] = 0.0; 
	N_matrix[1][0] = 0.0;       N_matrix[1][1] = Ne[0];		N_matrix[1][2] = 0.0;		N_matrix[1][3] = 0.0   ;		N_matrix[1][4] = Ne[1]; 	N_matrix[1][5] = 0.0;      N_matrix[1][6] = 0.0;			N_matrix[1][7] = Ne[2]; 	N_matrix[1][8] = 0.0;
	N_matrix[2][0] = 0.0;		N_matrix[2][1] = 0.0;		N_matrix[2][2] = Ne[0];		N_matrix[2][3] = 0.0   ;		N_matrix[2][4] = 0.0;	 	N_matrix[2][5] = Ne[1]; 	N_matrix[2][6] = 0.0;	 		N_matrix[2][7] = 0.0;		N_matrix[2][8] = Ne[2];	
}


double calculate3x3MatrixDeterminant(double matrix[3][3])
{
    return
        matrix[0][0] * matrix[1][1] * matrix[2][2]
        + matrix[0][1] * matrix[1][2] * matrix[2][0]
        + matrix[0][2] * matrix[1][0] * matrix[2][1]
        - matrix[0][0] * matrix[1][2] * matrix[2][1]
        - matrix[0][1] * matrix[1][0] * matrix[2][2]
        - matrix[0][2] * matrix[1][1] * matrix[2][0];
}

double invert3x3Matrix(double matrix_out[3][3],
                       double matrix[3][3])
{

    double determinant = calculate3x3MatrixDeterminant(matrix);
    double inverse_determinant = 1.0 / determinant;

    matrix_out[0][0] = (matrix[1][1] * matrix[2][2] - matrix[2][1] * matrix[1][2]) * inverse_determinant;
    matrix_out[0][1] = (matrix[2][1] * matrix[0][2] - matrix[0][1] * matrix[2][2]) * inverse_determinant;
    matrix_out[0][2] = (matrix[0][1] * matrix[1][2] - matrix[1][1] * matrix[0][2]) * inverse_determinant;

    matrix_out[1][0] = (matrix[1][2] * matrix[2][0] - matrix[2][2] * matrix[1][0]) * inverse_determinant;
    matrix_out[1][1] = (matrix[2][2] * matrix[0][0] - matrix[0][2] * matrix[2][0]) * inverse_determinant;
    matrix_out[1][2] = (matrix[0][2] * matrix[1][0] - matrix[1][2] * matrix[0][0]) * inverse_determinant;

    matrix_out[2][0] = (matrix[1][0] * matrix[2][1] - matrix[2][0] * matrix[1][1]) * inverse_determinant;
    matrix_out[2][1] = (matrix[2][0] * matrix[0][1] - matrix[0][0] * matrix[2][1]) * inverse_determinant;
    matrix_out[2][2] = (matrix[0][0] * matrix[1][1] - matrix[1][0] * matrix[0][1]) * inverse_determinant;

    return determinant;
}
