//接線剛性マトリクスのBTDBのBマトリクスを計算
void generate_linear_b_matrix(double (*b_t_matrix)[6], int point_n);

//接線剛性マトリクスの初期応力項中のBマトリクスを計算
double generate_nonlinear_b_matrix(double (*b_t_matrix)[9], int point_n);

//FPMの変位の離散化（u = {G}uE）につかうGマトリクスを計算
void calc_G(int dim, int point_n, double *point_xyz, int *support_offset, int *support, double **G);

//サブドメイン番号point_nの形状関数を計算
void calc_shape(double *xyz, int dim, int point_n, double *point_xyz, int *support_offset, double (*shapeF_t)[3]);

//サブドメイン番号point_nの試行関数を計算
void trial_u(double *xyz, int point_n, double *point_XYZ, double *u_h);