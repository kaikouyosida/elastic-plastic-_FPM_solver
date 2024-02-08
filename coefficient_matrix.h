void generate_coefficient_matrix();
void generate_subdomain_coefficient_matrix(int point_n, double (*ke_matrix)[60], 
                                            double (*current_deformation_gradients)[3], double *current_stress, double *trial_elastic_strains,
                                            double *equivalemt_plastic_strains, double *equivalent_plastic_strain_increments, double *back_stresses);
void generate_subdomain_coefficient_matrix_for_panaltyterm(int point_n1, int point_n2, int face_n, double (*ke_matrix)[60], 
                                            double (*current_deformation_gradients)[3], double *current_stress, double *trial_elastic_strains,
                                            double *equivalemt_plastic_strains, double *equivalent_plastic_strain_increments, double *back_stresses, double flag);
void generate_subdomain_coefficient_matrix_for_StabilizationTerm(int point_n1, int point_n2, int face_n, double (*ke_matrix)[60],int flag);
void assemble_coefficient_matrix_matrix_domain(double (*element_K)[60], double *Global_K, int point_n1, int point_n2);
double calc_area_change_factor(int subdomain_n1, int subdomain_n2, double (*Ne_d)[9]);
void generate_coefficient_linear();
void generate_Linear_coefficient_penalty(int face_n, int point_n1, int point_n2, double (*ke_matrix)[60], int flag);
void generate_Linear_coefficient_stabilization(int face_n, int point_n1, int point_n2, double (*ke_matrix)[60], int flag);