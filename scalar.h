//スカラー値の逆数をとる
double calculateInverse(double value);

//スカラー値の二乗をとる
double calculateSquare(double value);

//スカラー値の3乗をとる
double calculateCube(double value);

//value1とvalue2を入れ替える（double用）
void swapReals(double *value1, double *value2);

////value1とvalue2を入れ替える（int用）
void swapIntegers(int *value1, int *value2);

//内積の計算
double dot_product(int N, double *vec1, double *vec2);

//外積の計算
void cross_product(int dim, double *vecA, double *vecB, double *AcrossB);

//サブドメインの体積を計算(サブドメインを8分割→分割した四角すいを三角錐に分割→三角錐の総和をとってサブドメイン全体の体積を計算)
double calc_subdomain_volume(int point_n);
double calc_initial_subdomain_volume(int point_n);

//物理空間座標→正規化座標に変換するためのスカラー値を計算
double calc_mapping_parameter_for_av_area(const double face_node_XYZ[4][3], int s, int t, double *X);

double calc_mapping_parameter(int face_n, int point_n, int s, int t, double *X);

//ノルムの計算
double norm(double *vec, int n);

//マトリクス用のノルム
double norm_for_mat(double **vec,int m ,int n);

//ポイント間の距離を計算
double distance(int dim, int i, int j, double *point_xyz);

//3×3行列のトレースをとる
double calculate3x3MatrixTrace(double matrix[3][3]);

//歪みエネルギーノルムの相対誤差の計算
double generate_strain_energy_rate_parameter(int time_step);
