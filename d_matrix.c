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
    const double bulk_modulus = global.material.E_mod / (3.0 * (1.0 - 2.0 * global.material.nu_mod));
    const double shear_modulus = 0.5 * global.material.E_mod / (1.0 + global.material.nu_mod);
    const double hardening_modulus = get_hardening_modulus(equivalent_plastic_strain + equivalent_plastic_strain_increment);

    double temp;
    double coefficient;
    double trial_volumetric_strain;
    double trial_hydrostatic_stress;
    double trial_relative_equivalent_stress;
    double trial_relative_deviatoric_stresses[6];
    double trial_relative_stresses[6];

    //試行弾性ひずみの体積成分
    trial_volumetric_strain
        = trial_elastic_strains[0]
        + trial_elastic_strains[1]
        + trial_elastic_strains[2];
    
    //静水圧の計算
    trial_hydrostatic_stress = bulk_modulus * trial_volumetric_strain;

    //試行相対応力の偏差成分を計算
    for(int i = 0; i < 3; i++)
        trial_relative_deviatoric_stresses[i]
            = 2.0 * shear_modulus
            * (trial_elastic_strains[i]
               - trial_volumetric_strain / 3.0)
            - back_stresses[i];
    for(int i = 3; i < 6; i++)
        trial_relative_deviatoric_stresses[i]
            = 2.0 * shear_modulus
            * 0.5 * trial_elastic_strains[i]
            - back_stresses[i];

    //試行相対応力の計算
    for (int i = 0; i < 3; i++)
        trial_relative_stresses[i]
            = trial_relative_deviatoric_stresses[i]
            + trial_hydrostatic_stress;
    for (int i = 3; i < 6; i++)
        trial_relative_stresses[i]
            = trial_relative_deviatoric_stresses[i];
    
    //試行相当相対応力の計算
    trial_relative_equivalent_stress = calc_equivalent_stress(trial_relative_stresses);

    //弾性Dマトリクスの計算
    generateElasticDMatrix(d_matrix);

    //弾塑性Dマトリクスの計算
    temp = -equivalent_plastic_strain_increment
            * 6.0 * shear_modulus * shear_modulus / trial_relative_equivalent_stress;
    
    for(int i = 0; i < 3; i++)
        d_matrix[i][i] += temp;
    for(int i = 3; i < 6; i++)
        d_matrix[i][i] += 0.5 * temp;
    for(int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            d_matrix[i][j] -= temp / 3.0;
    
    coefficient
        = 9.0 * shear_modulus * shear_modulus
        * (equivalent_plastic_strain_increment / trial_relative_equivalent_stress
           - 1.0 / (3.0 * shear_modulus + hardening_modulus))
        / (trial_relative_equivalent_stress * trial_relative_equivalent_stress);
    for (int i = 0; i < 6; i++)
        for (int j = 0; j < 6; j++)
            d_matrix[i][j]
                += coefficient
                *  trial_relative_deviatoric_stresses[i]
                *  trial_relative_deviatoric_stresses[j];
}

//ペナルティ項用のコンシステント弾性Dマトリクスの計算
void modify_d_matrix_with_finite_strain_for_PenaltyTerm(double (*c_matrix)[9], double (*d_matrix)[6], double *current_stresses, double *trial_elastic_strains, double (*current_deformation_gradient)[3]){
    double consistent_d_tensor[3][3][3][3];
    double F_invF[3][3][3][3];
    double d_tensor[3][3][3][3];
    double l_tensor[3][3][3][3];
    double b_tensor[3][3][3][3];
    double current_stress_tensor[3][3];
    double trial_elastic_left_cauchy_green_deformations[3][3];
    double trial_elastic_strain_tensor[3][3];
    double inverse_deformation_grad[3][3];
    double inverse_volume_change;

    //変形勾配テンソルの逆行列
    invert3x3Matrix(inverse_deformation_grad, current_deformation_gradient);

    //DマトリクスからDテンソルへ変換
    convertSymmetric4thOrderMatrixToTensor(d_tensor, d_matrix);

    //試行弾性左コーシーグリーンテンソルの計算(B=exp[2ε^e])
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

    for (int i = 0; i < option.dim; i++)
        for (int j = 0; j < option.dim; j++)
            for (int k = 0; k < option.dim; k++)
                for (int l = 0; l < option.dim; l++)
                    b_tensor[i][j][k][l]
                        = identity_tensor[i][k] * trial_elastic_left_cauchy_green_deformations[j][l]
                        + identity_tensor[j][k] * trial_elastic_left_cauchy_green_deformations[i][l];
    
    inverse_volume_change = 1.0 / calc_3x3matrix_determinant(current_deformation_gradient);

    //応力をvoigt表記からテンソル表記へ変換,Cauchy応力を計算
    current_stress_tensor[0][0] = current_stresses[0];
    current_stress_tensor[0][1] = current_stresses[3];
    current_stress_tensor[0][2] = current_stresses[5];
    current_stress_tensor[1][0] = current_stresses[3];
    current_stress_tensor[1][1] = current_stresses[1];
    current_stress_tensor[1][2] = current_stresses[4];
    current_stress_tensor[2][0] = current_stresses[5];
    current_stress_tensor[2][1] = current_stresses[4];
    current_stress_tensor[2][2] = current_stresses[2];

    //変形勾配テンソルFと逆テンソルのテンソル積を計算
    for (int i = 0; i < option.dim; i++)
        for (int j = 0; j < option.dim; j++)
            for (int k = 0; k < option.dim; k++)
                for (int l = 0; l < option.dim; l++)
                    F_invF[i][j][k][l] = current_deformation_gradient[i][j] * inverse_deformation_grad[k][l];

    //接線係数の計算
    for (int i = 0; i < option.dim; i++)
        for (int j = 0; j < option.dim; j++)
            for (int k = 0; k < option.dim; k++)
                for (int l = 0; l < option.dim; l++)
                {
                    double d_i_j_k_l = 0.;

                    for(int ii = 0; ii < option.dim; ii++){
                        for(int jj = 0; jj < option.dim; jj++){
                            double F_invF_D_L_i_j_ii_jj = 0.;

                            for(int kk = 0; kk < option.dim; kk++){
                                for(int ll = 0; ll < option.dim; ll++){
                                    double F_invF_D_i_j_kk_ll = 0.;

                                    for(int mm = 0; mm < option.dim; mm++){
                                        for(int nn = 0; nn < option.dim; nn++){
                                            F_invF_D_i_j_kk_ll += F_invF[i][mm][j][nn] * d_tensor[mm][nn][kk][ll];
                                        }    
                                    }

                                    F_invF_D_L_i_j_ii_jj += F_invF_D_i_j_kk_ll * l_tensor[kk][ll][ii][jj];
                                }
                            }

                            d_i_j_k_l += F_invF_D_L_i_j_ii_jj * b_tensor[ii][jj][k][l];
                        }
                    }

                    d_i_j_k_l
                    *= 0.5 * inverse_volume_change;

                    double F_current_stress_invF_i_l = 0.;
                    for(int ii = 0; ii < option.dim; ii++){
                        F_current_stress_invF_i_l += current_deformation_gradient[i][ii] * current_stress_tensor[ii][l];
                    }
                    d_i_j_k_l += F_current_stress_invF_i_l * inverse_deformation_grad[j][k];

                    consistent_d_tensor[i][j][k][l] = d_i_j_k_l;
                }
    
    conver4thOrderTensorToMatrix(c_matrix, consistent_d_tensor);
}