//マトリクスの動的確保
double** matrix(int m, int n);

//マトリクスのメモリ開放
void free_matrix(double **A);

//３階の配列の動的確保
double ***threetimes_tensor(int m, int n, int l);

//３階の配列のメモリ開放
void free_tensor(double ***A);

//3×3の単位行列
void identify3x3Matrix(double (*matrix)[3]);

//3×3のゼロ行列
void zeroize3x3Matrix(double (*matrix)[3]);

//3×3の逆行列
void inverse_mat3x3(int dim, double (*A)[3], double (*inverse_A)[3]);

//3×3のマトリクスの二乗
void calculate3x3MatrixSquare(double matrix_out[3][3],
                              double matrix_in[3][3]);

//マトリクスどうしの積
void multi_mat(int Am, int Bn, double** A, double** B, int AnBm, double** X);

//マトリクスの転置
void trans_mat(int m, int n, double** X, double** XT);

//マトリクスの逆行列
void inverse_mat(int mn, double** A, double** inverse_A);

//３×３マトリクスのノルム
double calc_3x3matrix_norm(double (*matrix)[3]);

//3×3マトリクスのデターミナント
double calc_3x3matrix_determinant(double (*matrix)[3]);

//法線ベクトルの計算, それぞれの用途に合った形のマトリクスで用意
void calc_Ne(int dim, int subdomain_n1, int subdomain_n2, int face, int *vertex_offset, int *node, double *node_xyz, double *center_xyz, double (*N_matrix)[6]);
void calc_Ne_diagonal(int dim , int subdomain_n1, int subdomain_n2, int face, int *vertex_offset, int *node, double *node_xyz, double *center_xyz, double (*Ne_d)[9]);

//3×3の逆行列
double invert3x3Matrix(double matrix_out[3][3],
                       double matrix[3][3]);