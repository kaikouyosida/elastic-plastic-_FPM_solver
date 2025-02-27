//変形勾配, 対数ひずみ, 応力の更新, 内力ベクトルの計算
void update_field_and_internal_forces();

//相当塑性ひずみ増分の計算
double calc_equivalent_plastic_strain_increment(double trial_relative_equivalent_stress,
                                                double equivalent_plastic_strain,
                                                double yield_stress);

//内力ベクトルの体積積分項を計算
void calc_internal_force_volume(double **all_stress);

//内力ベクトルのペナルティ項を計算
void calc_internal_force_penalty(double **current_stress);

//内力ベクトルの安定化項を計算
void calc_internal_force_penalty_stabilization();

void update_field_and_internal_infinitesimal();

//残差ベクトルの計算
double calc_global_force_residual_norm(int iteration_step);

//ポイント変位の増分を更新
void update_point_displaecment_increment(double *du);

//節点変位の増分を更新
void update_nodal_displacement_by_current_NT(double *current_point_xyz);
void update_nodal_displacement_increment(double *current_point_XYZ, double *Initial_point_xyz);


//変位、応力等などの更新
void increment_field();

//変位ノルムの計算
double incremental_deformation_norm();

//塑性ひずみの更新
void update_plastic_strains(double plastic_strains[6], const double stresses[6], const double equivarent_stresses, const double equivarent_strain_increment);

//サブドメインの中にポイントが１つづつ含まれているか判定
void whether_points_is_in_the_subdomain();

//反復回数が既定の回数に到達したら変位増分などを初期化して再計算
void cut_back();