void generate_coefficient_matrix();
void generate_subdomain_coefficient_matrix(int point_n, double (*ke_matrix)[60], 
                                            double (*current_deformation_gradients)[3], double *current_stress, double *trial_elastic_strains,
                                            double *equivalemt_plastic_strains, double *equivalent_plastic_strain_increments, double *back_stresses);
void assemble_coefficient_matrix_matrix_domain(double (*element_K)[60], double *Global_K, int point_n1, int point_n2);