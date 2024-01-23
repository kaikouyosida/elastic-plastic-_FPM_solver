double** matrix(int m, int n);
void free_matrix(double **A);
double ***threetimes_tensor(int m, int n, int l);
void free_tensor(double ***A);
void identify3x3Matrix(double (*matrix)[3]);
void multi_mat(int Am, int Bn, double** A, double** B, int AnBm, double** X);
void trans_mat(int m, int n, double** X, double** XT);
void inverse_mat(int mn, double** A, double** inverse_A);
