//INTEL MKLの演算機能を使用して連立一次方程式を求解
//LU分解で密行列を解く
void LU_dens(double *global_k, double *du, double *Q, int DoF_free);

//LU分解でスパース行列を解く
void Paradiso(int n, int nnz, double *a, int *ia, int *ja, double *b, double *x);

//LAPACKE_dgesvdの計算を行い、ランクを計算する
int calculate_rank(double *matrix, int n);