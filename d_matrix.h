void generateElasticDMatrix(double (*d_matrix)[6]);
void modify_d_matrix_with_finite_strain(double (*d_matrix)[6], double *current_stresses, double *trial_elastic_strains, double (*current_deformation_gradient)[3]);
void modify_d_matrix_with_finite_strain_for_PenaltyTerm(double (*d_matrix)[6], double *current_stresses, double *trial_elastic_strains, double (*current_deformation_gradient)[3]);