//弾性Dマトリクスの計算（一般化hookの法則)
void generateElasticDMatrix(double (*d_matrix)[6]);

//弾塑性Dマトリクスの計算
void generate_elastic_plastic_d_matrix(double d_matrix[6][6], const double trial_elastic_strains[6], const double equivalent_plastic_strain, const double equivalent_plastic_strain_increment, const double back_stresses[6]);

//update Lagrangeの有限ひずみに適用するために修正されたDマトリクス
void modify_d_matrix_with_finite_strain(double (*d_matrix)[6], const double *current_stresses, const double *trial_elastic_strains, double (*current_deformation_gradient)[3]);

//update Lagrangeの有限ひずみに適用するために修正されたDマトリクス(ペナルティ項用)
void modify_d_matrix_with_finite_strain_for_PenaltyTerm(double (*c_matrix)[9], double (*d_matrix)[6], double *current_stresses, double *trial_elastic_strains, double (*current_deformation_gradient)[3]);