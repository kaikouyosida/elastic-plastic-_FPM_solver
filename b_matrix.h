void generate_linear_b_matrix(double (*b_t_matrix)[6], int point_n);
void calc_G(int dim, int point_n, double *point_xyz, int *support_offset, int *support, double **G);
void calc_shape(double *xyz, int dim, int point_n, double *point_xyz, int *support_offset, double (*shapeF_t)[3]);