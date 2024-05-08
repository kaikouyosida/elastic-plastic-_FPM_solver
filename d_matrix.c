#include<stdio.h>
#include<stdlib.h>
#include"type.h"
#include"tensor.h"
#include"matrix.h"

extern Global global;
extern Option option;

static const double identity_tensor[3][3] = {{1.0, 0.0, 0.0},
                                             {0.0, 1.0, 0.0},
                                             {0.0, 0.0, 1.0}};


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

void modify_d_matrix_with_finite_strain(double (*d_matrix)[6], double *current_stresses, double *trial_elastic_strains, double (*current_deformation_gradient)[3]){
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

    //for(int i = 0; i < 6; i++)
        //printf("%+15.14e    ", trial_elastic_strains[i]);
    
    //printf("\n");

    //for(int i = 0; i < 6; i++)
        //printf("%+15.14e    ", current_stresses[i]);
    //printf("\n");
    //for(int i = 0; i < 3; i++){
        //for(int j = 0; j < 3; j++){
            //printf("%+8.7e  ", current_deformation_gradient[i][j]);
        //}
        //printf("\n");
    //}
    //printf("\n");

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
void modify_d_matrix_with_finite_strain_for_PenaltyTerm(double (*c_matrix)[9], double (*d_matrix)[6], double *current_stresses, double *trial_elastic_strains, double (*current_deformation_gradient)[3]){
    double consistent_d_tensor[3][3][3][3];
    double A_tensor[3][3][3][3];
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

    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            for (int k = 0; k < 3; k++)
                for (int l = 0; l < 3; l++)
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

    //接線係数の計算
    for (int i = 0; i < option.dim; i++)
        for (int j = 0; j < option.dim; j++)
            for (int k = 0; k < option.dim; k++)
                for (int l = 0; l < option.dim; l++)
                {
                    double d_i_j_k_l = 0.0;

                    for (int ii = 0; ii < option.dim; ii++)
                        for (int jj = 0; jj < option.dim; jj++)
                        {
                            double d_ii_jj_kk_ll_times_b_kk_ll_k_l = 0.0;

                            for (int kk = 0; kk < option.dim; kk++)
                                for (int ll = 0; ll < option.dim; ll++)
                                    d_ii_jj_kk_ll_times_b_kk_ll_k_l
                                        += l_tensor[ii][jj][kk][ll]
                                        *  b_tensor[kk][ll][k][l];

                            d_i_j_k_l
                                += d_tensor[i][j][ii][jj]
                                *  d_ii_jj_kk_ll_times_b_kk_ll_k_l;
                        }

                    d_i_j_k_l *= 0.5 * inverse_volume_change;

                    A_tensor[i][j][k][l]
                        = d_i_j_k_l;
                }
    
    for(int i = 0; i < option.dim; i++){
        for(int j = 0; j < option.dim; j++){
            for(int k = 0; k < option.dim; k++){
                for(int l = 0; l < option.dim; l++){
                    double scalar1 = 0.;

                    for(int ii = 0; ii < option.dim; ii++){
                        double scalar2 = 0.;

                        for(int jj = 0; jj < option.dim; jj++){
                            double scalar3 = 0.;

                            for(int kk = 0; kk < option.dim; kk++){
                                scalar3 += A_tensor[ii][jj][k][kk] * current_deformation_gradient[l][kk];
                            }
                            scalar2 += inverse_deformation_grad[j][jj] * scalar3;
                        }
                        scalar1 += current_deformation_gradient[i][ii] * scalar2;
                   }

                   double scalar4 = 0.;
                   for(int ii = 0; ii < option.dim; ii++){
                       scalar4 += current_deformation_gradient[i][ii] * current_stress_tensor[ii][l];
                   }
                   scalar1 -= inverse_deformation_grad[j][k] * scalar4;

                   consistent_d_tensor[i][j][k][l] = scalar1;
                }
            }
        }
    }
    
    conver4thOrderTensorToMatrix(c_matrix, consistent_d_tensor);
}