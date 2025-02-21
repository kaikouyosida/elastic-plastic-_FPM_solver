#include<stdio.h>
#include<stdlib.h>
#include"type.h"
#include"tensor.h"
#include"matrix.h"
#include"stress.h"
#include"ss_curve.h"

extern Global global;
extern Option option;
extern SS_CURVE ss_curve;

double identity_tensor[3][3] = {{1.0, 0.0, 0.0},
                                {0.0, 1.0, 0.0},
                                {0.0, 0.0, 1.0}};

//弾性Dマトリクスの計算
void generateElasticDMatrix(double (*d_matrix)[6]){
    double young_modulus
        = global.material.E_mod;
    double poisson_ratio
        = global.material.nu_mod;


    double d00
        = young_modulus * (1.0 - poisson_ratio)
        / (1.0 + poisson_ratio) / (1.0 - 2.0 * poisson_ratio);
    double d01
        = young_modulus * poisson_ratio
        / (1.0 + poisson_ratio) / (1.0 - 2.0 * poisson_ratio);
    double d33
        = 0.5 * young_modulus / (1.0 + poisson_ratio);

    d_matrix[0][0] = d00; d_matrix[0][1] = d01; d_matrix[0][2] = d01;
    d_matrix[0][3] = 0.0; d_matrix[0][4] = 0.0; d_matrix[0][5] = 0.0;

    d_matrix[1][0] = d01; d_matrix[1][1] = d00; d_matrix[1][2] = d01;
    d_matrix[1][3] = 0.0; d_matrix[1][4] = 0.0; d_matrix[1][5] = 0.0;

    d_matrix[2][0] = d01; d_matrix[2][1] = d01; d_matrix[2][2] = d00;
    d_matrix[2][3] = 0.0; d_matrix[2][4] = 0.0; d_matrix[2][5] = 0.0;

    d_matrix[3][0] = 0.0; d_matrix[3][1] = 0.0; d_matrix[3][2] = 0.0;
    d_matrix[3][3] = d33; d_matrix[3][4] = 0.0; d_matrix[3][5] = 0.0;

    d_matrix[4][0] = 0.0; d_matrix[4][1] = 0.0; d_matrix[4][2] = 0.0;
    d_matrix[4][3] = 0.0; d_matrix[4][4] = d33; d_matrix[4][5] = 0.0;

    d_matrix[5][0] = 0.0; d_matrix[5][1] = 0.0; d_matrix[5][2] = 0.0;
    d_matrix[5][3] = 0.0; d_matrix[5][4] = 0.0; d_matrix[5][5] = d33;
}

//コンシステント弾性Dマトリクスの計算
void modify_d_matrix_with_finite_strain(double (*d_matrix)[6], const double *current_stresses, const double *trial_elastic_strains, double (*current_deformation_gradient)[3]){
    double consistent_d_tensor[3][3][3][3];
    double d_tensor[3][3][3][3];
    double l_tensor[3][3][3][3];
    double b_tensor[3][3][3][3];
    double current_stress_tensor[3][3];
    double trial_elastic_left_cauchy_green_deformations[3][3];
    double trial_elastic_strain_tensor[3][3];
    double inverse_volume_change;

    //DマトリクスからDテンソルへ変換
    convertSymmetric4thOrderMatrixToTensor(d_tensor, d_matrix);

    //試行弾性左コーシーグリーンテンソルの計算
    trial_elastic_strain_tensor[0][0] = 2.0 * trial_elastic_strains[0];
    trial_elastic_strain_tensor[0][1] = 2.0 * 0.5 * trial_elastic_strains[3];
    trial_elastic_strain_tensor[0][2] = 2.0 * 0.5 * trial_elastic_strains[5];
    trial_elastic_strain_tensor[1][0] = 2.0 * 0.5 * trial_elastic_strains[3];
    trial_elastic_strain_tensor[1][1] = 2.0 * trial_elastic_strains[1];
    trial_elastic_strain_tensor[1][2] = 2.0 * 0.5 * trial_elastic_strains[4];
    trial_elastic_strain_tensor[2][0] = 2.0 * 0.5 * trial_elastic_strains[5];
    trial_elastic_strain_tensor[2][1] = 2.0 * 0.5 * trial_elastic_strains[4];
    trial_elastic_strain_tensor[2][2] = 2.0 * trial_elastic_strains[2];

    calculateTensorExponent(trial_elastic_left_cauchy_green_deformations,
                            trial_elastic_strain_tensor);
    
    //[L] = d(ln([B]^trial))/d[B]^trialの計算
    calculateTensorLogarithmDerivative(l_tensor,
                                       trial_elastic_left_cauchy_green_deformations);

    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            for (int k = 0; k < 3; k++)
                for (int l = 0; l < 3; l++)
                    b_tensor[i][j][k][l]
                        = identity_tensor[i][k] * trial_elastic_left_cauchy_green_deformations[j][l]
                        + identity_tensor[j][k] * trial_elastic_left_cauchy_green_deformations[i][l];
    
    //応力をvoigt表記からテンソル表記へ変換
    current_stress_tensor[0][0] = current_stresses[0];
    current_stress_tensor[0][1] = current_stresses[3];
    current_stress_tensor[0][2] = current_stresses[5];
    current_stress_tensor[1][0] = current_stresses[3];
    current_stress_tensor[1][1] = current_stresses[1];
    current_stress_tensor[1][2] = current_stresses[4];
    current_stress_tensor[2][0] = current_stresses[5];
    current_stress_tensor[2][1] = current_stresses[4];
    current_stress_tensor[2][2] = current_stresses[2];

    inverse_volume_change = 1.0 / calc_3x3matrix_determinant(current_deformation_gradient);

    // 1 / (2 * J) * ([D] : [L] : [B])_ijkl - (sigma_il * delta_jk + sigma_jl * delta_ik) を計算
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            for (int k = 0; k < 3; k++)
                for (int l = 0; l < 3; l++)
                {
                    double d_i_j_k_l = 0.0;

                    for (int ii = 0; ii < 3; ii++)
                        for (int jj = 0; jj < 3; jj++)
                        {
                            double d_ii_jj_kk_ll_times_b_kk_ll_k_l = 0.0;

                            for (int kk = 0; kk < 3; kk++)
                                for (int ll = 0; ll < 3; ll++)
                                    d_ii_jj_kk_ll_times_b_kk_ll_k_l
                                        += l_tensor[ii][jj][kk][ll]
                                        *  b_tensor[kk][ll][k][l];

                            d_i_j_k_l
                                += d_tensor[i][j][ii][jj]
                                *  d_ii_jj_kk_ll_times_b_kk_ll_k_l;
                        }

                    d_i_j_k_l
                        *= 0.5 * inverse_volume_change;

                    d_i_j_k_l
                        -= current_stress_tensor[i][l] * identity_tensor[j][k];
                    d_i_j_k_l
                        -= current_stress_tensor[j][l] * identity_tensor[i][k];

                    consistent_d_tensor[i][j][k][l]
                        = d_i_j_k_l;
                }

    //Dマトリクスをテンソルからマトリクスに変換
    convertSymmetric4thOrderTensorToMatrix(d_matrix,
                                           consistent_d_tensor);
}

//材料非線形性のDマトリクスの計算
void generate_elastic_plastic_d_matrix(double d_matrix[6][6], const double trial_elastic_strains[6], const double equivalent_plastic_strain, const double equivalent_plastic_strain_increment, const double back_stresses[6]){
    //体積弾性係数
    const double bulk_modulus = global.material.E_mod / (3.0 * (1.0 - 2.0 * global.material.nu_mod));
    //せん断弾性係数
    const double shear_modulus = 0.5 * global.material.E_mod / (1.0 + global.material.nu_mod);
    //硬化係数
    const double hardening_modulus = get_hardening_modulus(equivalent_plastic_strain + equivalent_plastic_strain_increment);

    double trial_relative_equivalent_stress;        //試行相対相当応力
    double trial_relative_deviatoric_stresses[6];   //試行偏差応力
    double trial_relative_stresses[6];              //試行応力
    double hydro_static_stress;                     //静水圧
    double elastic_d_matrix[6][6];                  //弾性Dマトリクス

    //l_d = 0.5 * (δ_ik * δ_jl + δ_il * δ_jk) - δ_ij * δ_kl / 3.0
    double l_d[6][6] = {{2.0 / 3.0, -1.0 / 3.0, -1.0 / 3.0, 0, 0, 0}, 
                        {-1.0 / 3.0, 2.0 / 3.0, -1.0 / 3.0, 0, 0, 0}, 
                        {-1.0 / 3.0, -1.0 / 3.0, 2.0 / 3.0, 0, 0, 0},
                        {0, 0, 0, 1.0 / 2.0, 0, 0}, 
                        {0, 0, 0, 0, 1.0 / 2.0, 0}, 
                        {0, 0, 0, 0, 0, 1.0 / 2.0}};

    //δ = δij * δkl
    double delta_matrix[6][6] = {{1.0, 1.0, 1.0, 0.0, 0.0, 0.0}
                                , {1.0, 1.0, 1.0, 0.0, 0.0, 0.0}
                                , {1.0, 1.0, 1.0, 0.0, 0.0, 0.0}
                                , {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}
                                , {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}
                                , {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};

    
    //試行応力の計算
    generateElasticDMatrix(elastic_d_matrix);

    for(int i = 0; i < 6; i++){
        double trial_relative_stresses_i = 0;
        for(int j = 0; j < 6; j++){
            trial_relative_stresses_i += elastic_d_matrix[i][j] * trial_elastic_strains[j];
        }
        trial_relative_stresses[i] = trial_relative_stresses_i;
    }

    //試行相当応力の計算
    trial_relative_equivalent_stress = calc_equivalent_stress(trial_relative_stresses);
    
    //静水圧, 偏差応力の計算
    hydro_static_stress = (trial_relative_stresses[0] + trial_relative_stresses[1] + trial_relative_stresses[2]) / 3.0;
    for(int i = 0; i < 3; i++)
        trial_relative_deviatoric_stresses[i] = trial_relative_stresses[i] - hydro_static_stress;
    for(int i = 3; i < 6; i++)
        trial_relative_deviatoric_stresses[i] = trial_relative_stresses[i];

    //弾塑性Dマトリクスの計算
    for(int i = 0; i < 6; i++){
        for(int j = 0; j < 6; j++){
            d_matrix[i][j] = 2.0 * shear_modulus
                             * (1.0 - 3.0 * shear_modulus * equivalent_plastic_strain_increment / trial_relative_equivalent_stress) * l_d[i][j]

                             + 9.0 * shear_modulus * shear_modulus
                             * (equivalent_plastic_strain_increment / trial_relative_equivalent_stress - 1.0 /(3.0 * shear_modulus + hardening_modulus))
                             / trial_relative_equivalent_stress / trial_relative_equivalent_stress * trial_relative_deviatoric_stresses[i] * trial_relative_deviatoric_stresses[j]
                             
                             + bulk_modulus * delta_matrix[i][j];
        }
    }
}
