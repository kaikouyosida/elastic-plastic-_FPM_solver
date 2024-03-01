//Dマトリクスの計算（Henkeyモデル）
void generateElasticDMatrix(double (*d_matrix)[6]);

//update Lagrangeの有限ひずみに適用するために修正されたDマトリクス
void modify_d_matrix_with_finite_strain(double (*d_matrix)[6], double *current_stresses, double *trial_elastic_strains, double (*current_deformation_gradient)[3]);

//update Lagrangeの有限ひずみに適用するために修正されたDマトリクス(ペナルティ項用)
void modify_d_matrix_with_finite_strain_for_PenaltyTerm(double (*c_matrix)[9], double (*d_matrix)[6], double *current_stresses, double *trial_elastic_strains, double (*current_deformation_gradient)[3]);