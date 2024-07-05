//変形勾配, 対数ひずみ, 応力の更新, 内力ベクトルの計算
void update_field_and_internal_forces();

//相当塑性ひずみ増分の計算
double calc_equivalent_plastic_strain_increment(double trial_relative_equivalent_stress,
                                                double equivalent_plastic_strain,
                                                double yield_stress);

//内力ベクトルの体積積分項を計算
void calc_internal_force_volume(double **all_stress);

//内力ベクトルのペナルティ項を計算
void calc_internal_force_penalty(double **all_stress);

//内力ベクトルの安定化項を計算
void calc_internal_force_penalty_stabilization();

//残差ベクトルの計算
double calc_global_force_residual_norm(int iteration_step);

//ポイント変位の増分を更新
void update_point_displaecment_increment(double *current_point_xyz, double *du);

//節点変位の増分を更新
void update_nodal_displacement_increment(double *latest_point_xyz);

//変位、応力等などの更新
void increment_field();