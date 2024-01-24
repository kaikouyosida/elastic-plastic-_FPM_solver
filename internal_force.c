#include<stdio.h>
#include<stdlib.h>
#include"type.h"
#include"d_matrix.h"
#include"b_matrix.h"
#include"matrix.h"
#include"tensor.h"
#include"stress.h"
extern Global global;
extern Option option;


void zero_fill_displacement_increments(){
    for(int point = 0; point < global.subdomain.N_point; point++){
        for(int i  = 0; i < option.dim; i++){
            global.subdomain.displacement[point][i] = 0.;
        }
    }
}

void update_field_and_internal_forces(){
    FILE *fp_debug;                                 //デバッグ用のファイル
    fp_debug = fopen("debag.dat", "w");

    double subdomain_internal_force[20][3];
    double d_matrix[6][6];          //Dマトリクス
    double inverse_relative_deformation_gradient[3][3];         //相対変形勾配テンソルの逆テンソル
    double relative_deformation_gradient[3][3];                 //相対変形勾配テンソル
    double elastic_strain_tensor[3][3];                         //弾性ひずみテンソル
    double deformation_gradients[3][3];                         //変形勾配テンソル
    double current_deformation_gradients[3][3];                 //現配置の変形勾配テンソル
    double elastic_left_cauchy_green_deformations[3][3];        //弾性左コーシーグリーンテンソル
    double trial_elastic_left_cauchy_green_deformations[3][3];  //試行弾性左コーシーグリーンテンソル
    double trial_relative_stresses[6];                          //試行相対応力
    double b_t_matrix[60][6];                                   //Bマトリクスの転置
    int support[60];                                         //サポートドメイン内のポイント数
    double displacement_increment[3];                           //サポートの変位増分
    double X[27][3];
    double w[27];                                               //ガウス求積に使う正規化座標と重み関数
    double elastic_strains[6];                                  //弾性ひずみ
    double current_elastic_strains[6];                          //現配置の弾性ひずみ
    double trial_elastic_strains[6];                            //試行弾性ひずみ
    double current_stresses[6];                                 //現配置の応力
    double back_stresses[6];                                    //背応力
    double current_back_stresses[6];                            //現配置の背応力

    //internal_forceをゼロ処理
    for(int i = 0; i < global.subdomain.N_point; i++){
        for(int j = 0; j < option.dim; j++){
            global.subdomain.internal_force[i][j] = 0.;
        }
    }

    for(int point = 0; point < global.subdomain.N_point; point++){
        double kinematic_hardening_fractions = global.material.kinematic_hardening_fractions;
        
        double N_support = global.subdomain.support_offset[point + 1] - global.subdomain.support_offset[point];
        for(int i = 0 ; i < N_support; i++)
            support[i] = global.subdomain.support[option.dim * global.subdomain.support_offset[point] + i];

        //サブドメインの内力ベクトルをゼロ処理
        for(int i = 0; i < N_support; i++){
            for(int j = 0; j < option.dim; j++){
                subdomain_internal_force[i][j] = 0.;
            }
        }
        for(int i = 0; i < 3; i++){
            for(int j = 0; j < 3; j++){
                deformation_gradients[i][j] = global.subdomain.deformation_gradients[i][j][point];
                current_deformation_gradients[i][j] = global.subdomain.current_deformation_gradients[i][j][point];           
            }
        }
        for(int i = 0; i < 6; i++){
            elastic_strains[i] = global.subdomain.elastic_strains[point][i];
            current_elastic_strains[i] = global.subdomain.current_elastic_strains[point][i];
            trial_elastic_strains[i] = global.subdomain.trial_elastic_strains[point][i];
            current_stresses[i] = global.subdomain.current_stresses[point][i];
            back_stresses[i] = global.subdomain.back_stresses[point][i];
            current_back_stresses[i] = global.subdomain.current_back_stresses[point][i];
        }
        double *equivalent_plastic_strain = &global.subdomain.equivalent_plastic_strains[point];
        double *equivalent_plastic_strain_increment = &global.subdomain.equivalent_plastic_strain_increments[point];
        double *yield_stress = &global.subdomain.yield_stresses[point];
        double *current_yield_stress = &global.subdomain.current_yield_stresses[point];
        double trial_relative_equivalent_stress;
        double factor;

        generate_linear_b_matrix(b_t_matrix, point);
        
        generateElasticDMatrix(d_matrix);

        //相対変形勾配テンソルを計算。dF=(I-d(du)/d(x+u+du))^-1//
        identify3x3Matrix(inverse_relative_deformation_gradient);
        //inverse_relative_deformation_gradientの対角項を計算
        for(int i = 0; i < N_support; i++){
            for(int j = 0; j < option.dim; j++)
                displacement_increment[j] = global.subdomain.displacement_increment[support[i]][j];

            for(int j = 0; j < option.dim; j++){
                double displacement_increment_j
                    = displacement_increment[j];
                for (int k = 0; k < 3; k++)
                            inverse_relative_deformation_gradient[k][k]
                                -= b_t_matrix[3 * (i + 1) + j][k]
                                *  displacement_increment_j;
            }
        }
        for(int i = 0; i < option.dim; i++)
            displacement_increment[i] = global.subdomain.displacement[point][i];
        for(int j = 0; j < option.dim; j++){
            double displacement_increment_j = displacement_increment[j];
            inverse_relative_deformation_gradient[j][j] 
                -= b_t_matrix[j][j] * displacement_increment_j;
        }

        //inverse_relative_deformation_gradientの非対角項を計算
        for(int i = 0 ; i < N_support; i++){
            for(int j = 0; j < option.dim; j++)
                displacement_increment[j] = global.subdomain.displacement_increment[support[i]][j];

            for (int j = 0; j < 3; j++)
                    {
                        double displacement_increment_j = displacement_increment[j];

                        inverse_relative_deformation_gradient[j][(j + 1) % 3]
                            -= b_t_matrix[3 * (i + 1) + j][3 + j]
                            *  displacement_increment_j;
                        inverse_relative_deformation_gradient[j][(j + 2) % 3]
                            -= b_t_matrix[3 * (i + 1) + j][3 + (j + 2) % 3]
                            *  displacement_increment_j;
                    }
        }
        for(int i = 0; i < option.dim; i++)
                displacement_increment[i] = global.subdomain.displacement[point][i];
        for(int i = 0; i < option.dim; i++){
            double displacement_increment_i = displacement_increment[i];
            for(int j = 0; j < option.dim; j++){
                if(i != j) inverse_relative_deformation_gradient[i][j] -= b_t_matrix[j][j] * displacement_increment_i;
            }
        }

        inverse_mat3x3(option.dim, inverse_relative_deformation_gradient, relative_deformation_gradient);

        //変形勾配テンソル（初期配置に対してのテンソル）
        for (int i = 0; i < 3; i++)
                    for (int j = 0; j < 3; j++)
                    {
                        double deformation_gradient_i_j = 0.0;

                        for (int k = 0; k < 3; k++)
                            deformation_gradient_i_j
                                += relative_deformation_gradient[i][k]
                                *  deformation_gradients[k][j];

                        current_deformation_gradients[i][j]
                            = deformation_gradient_i_j;
                    }
        
        //弾性左コーシーグリーンテンソルの計算([B]^e = exp(2 * {epsilon}^e))
        elastic_strain_tensor[0][0] = 2.0 * elastic_strains[0];
        elastic_strain_tensor[0][1] = 2.0 * 0.5 * elastic_strains[3];
        elastic_strain_tensor[0][2] = 2.0 * 0.5 * elastic_strains[5];
        elastic_strain_tensor[1][0] = 2.0 * 0.5 * elastic_strains[3];
        elastic_strain_tensor[1][1] = 2.0 * elastic_strains[1];
        elastic_strain_tensor[1][2] = 2.0 * 0.5 * elastic_strains[4];
        elastic_strain_tensor[2][0] = 2.0 * 0.5 * elastic_strains[5];
        elastic_strain_tensor[2][1] = 2.0 * 0.5 * elastic_strains[4];
        elastic_strain_tensor[2][2] = 2.0 * elastic_strains[2];
        
        calculateTensorExponent(elastic_left_cauchy_green_deformations, elastic_strain_tensor);

        //試行弾性左コーシーグリーンテンソルの計算
        for (int i = 0; i < 3; i++)
                    for (int j = 0; j < 3; j++)
                    {
                        double left_cauchy_green_deformation_i_j = 0.0;

                        for (int k = 0; k < 3; k++)
                        {
                            double left_cauchy_green_deformation_k_l_times_relative_deformation_gradient_j_l = 0.0;

                            for (int l = 0; l < 3; l++)
                                left_cauchy_green_deformation_k_l_times_relative_deformation_gradient_j_l
                                    += elastic_left_cauchy_green_deformations[k][l]
                                    *  relative_deformation_gradient[j][l];

                            left_cauchy_green_deformation_i_j
                                += relative_deformation_gradient[i][k]
                                *  left_cauchy_green_deformation_k_l_times_relative_deformation_gradient_j_l;
                        }

                        trial_elastic_left_cauchy_green_deformations[i][j]
                            = left_cauchy_green_deformation_i_j;
                    }
        //試行弾性ひずみを計算（ {epsilon}^trial = 1/2 * ln([B]^trial)　
        calculateTensorLogarithm(elastic_strain_tensor, trial_elastic_left_cauchy_green_deformations);
        current_elastic_strains[0] = 0.5 * elastic_strain_tensor[0][0];
        current_elastic_strains[1] = 0.5 * elastic_strain_tensor[1][1];
        current_elastic_strains[2] = 0.5 * elastic_strain_tensor[2][2];
        current_elastic_strains[3] = 0.5 * (elastic_strain_tensor[0][1]
                                            + elastic_strain_tensor[1][0]);
        current_elastic_strains[4] = 0.5 * (elastic_strain_tensor[1][2]
                                            + elastic_strain_tensor[2][1]);
        current_elastic_strains[5] = 0.5 * (elastic_strain_tensor[2][0]
                                            + elastic_strain_tensor[0][2]);

        for (int i = 0; i < 6; i++)
                trial_elastic_strains[i]
                    = current_elastic_strains[i];
        
        //試行応力の計算({sigma}^trial = [D] * {epsilon}^trial)
        for (int i = 0; i < 6; i++)
            {
                double stress_i = 0.0;

                for (int j = 0; j < 6; j++)
                    stress_i += d_matrix[i][j] * current_elastic_strains[j];

                current_stresses[i] = stress_i;
            }

        //試行相対応力の計算({sigma}^trial = {sigma}^trial - {beta})
        for (int i = 0; i < 6; i++)
                trial_relative_stresses[i]
                    = current_stresses[i] - back_stresses[i];
        //試行相対相当応力の計算(sigma^trial_e)
        trial_relative_equivalent_stress
                = calc_equivalent_stress(trial_relative_stresses);
        if (trial_relative_equivalent_stress <= (*yield_stress))
            {   
                *equivalent_plastic_strain_increment = 0.0;
                *current_yield_stress = *yield_stress;
                for (int i = 0; i < 6; i++)
                    current_back_stresses[i] = back_stresses[i];
            }else{
                
            }

    }        
    fclose(fp_debug);
    exit(0);





}