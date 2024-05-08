//変形勾配, 対数ひずみ, 応力の更新, 内力ベクトルの計算
void update_field_and_internal_forces();

//変数の増分をゼロ処理
void zero_fill_displacement_increments();

//相当塑性ひずみ増分の計算
double calc_equivalent_plastic_strain_increment(double trial_relative_equivalent_stress,
                                                double equivalent_plastic_strain,
                                                double yield_stress);

//内力ベクトルのペナルティ項を計算
void calc_internal_force_penalty(double **all_stress,int N_qu);

//内力ベクトルのペナルティ項を計算（安定化項）
void calc_internal_force_penalty_stabilization(int N_qu);

//残差ベクトルの計算
double calc_global_force_residual_norm(int iteration_step);

//節点変位の増分を更新
void update_nodal_displacement_increment();