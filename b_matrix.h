//接線剛性マトリクスのBTDBのBマトリクスを計算
void generate_linear_b_matrix(double (*b_t_matrix)[6], const int point_n);

//接線剛性マトリクスの初期応力項中のBマトリクスを計算
double generate_nonlinear_b_matrix(double (*b_t_matrix)[9], const int point_n);

//FPMの変位の離散化（u = {G}uE）につかうGマトリクスを計算
void calc_G(const int dim, const int point_n, const double *point_xyz, const int *support_offset, const int *support, double **G);

//サブドメイン番号point_nの形状関数を計算
void calc_shape(const double *xyz, const int dim, const int point_n, const double *point_xyz, const int *support_offset, double (*shapeF_t)[3]);

//サブドメイン番号point_nの試行関数を計算(pm=0のときは変位, 1のときは変位増分を計算)
void trial_u(const double *xyz, const int point_n, const double *point_XYZ, double *u_h, const int pm);