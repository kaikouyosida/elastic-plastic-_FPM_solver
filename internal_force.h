void update_field_and_internal_forces();
void zero_fill_displacement_increments();
double calc_equivalent_plastic_strain_increment(double trial_relative_equivalent_stress,
                                                double equivalent_plastic_strain,
                                                double yield_stress);
void calc_internal_force_penalty(double **all_stress,int N_qu);
