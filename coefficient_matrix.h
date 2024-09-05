//接線剛性マトリクスを計算
void generate_coefficient_matrix();

//接線剛性マトリクスの領域積分項（BTDB+初期応力項）の計算
void generate_subdomain_coefficient_matrix_for_volume(const int point_n,
                                                    double (*current_deformation_gradients)[3], const double *current_stress, const double *trial_elastic_strains,
                                                    const double *equivalemt_plastic_strains, const double *equivalent_plastic_strain_increments, const double *back_stresses);

//接線剛性マトリクスの境界積分項を計算
void generate_subdomain_coefficient_matrix_for_PenaltyTerm(const int point_n1, const int point_n2, const int face_n,
                                            double (*current_deformation_gradients)[3], const double *current_stress, double *trial_elastic_strains,
                                            const double *equivalemt_plastic_strains, const double *equivalent_plastic_strain_increments, const double *back_stresses, const int flag);

//接線剛性マトリクスの境界積分項を計算（安定化項）
void generate_subdomain_coefficient_matrix_for_StabilizationTerm(const int point_n1, const int point_n2, const int face_n, const int flag);

//剛性マトリクスを要素→全体へ組み立てる
void assemble_coefficient_matrix(double (*element_K)[60], double *Global_K, int point_n1, int point_n2);

//線形弾性のルーチン
//係数マトリクスの計算
void generate_coefficient_linear();

//ペナルティ項を計算
void generate_Linear_coefficient_penalty(int face_n, int point_n1, int point_n2, double (*ke_matrix)[60], int flag);

//ペナルティ項を計算（安定化項）
void generate_Linear_coefficient_stabilization(int face_n, int point_n1, int point_n2, double (*ke_matrix)[60], int flag);