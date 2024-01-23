#include<stdio.h>
#include<stdlib.h>
#include"type.h"
#include"d_matrix.h"
#include"b_matrix.h"
#include"matrix.h"
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
        for(int i = 0; i < 3 * N_support + 3; i++){
            for(int j = 0; j < 6; j++){
                fprintf(fp_debug, "%+3.2e ", b_t_matrix[i][j]);
            }
            fprintf(fp_debug,"\n");
        }
            fprintf(fp_debug,"\n");
        
        
        generateElasticDMatrix(d_matrix);

        identify3x3Matrix(inverse_relative_deformation_gradient);

        for(int i = 0; i < N_support; i++){
            for(int j = 0; j < option.dim; j++)
                displacement_increment[j] = global.subdomain.displacement_increment[point][j];
            for(int j = 0; j < option.dim; j++){
                double displacement_increment_j
                    = displacement_increment[j];
            }
        }



    }
    fclose(fp_debug);
    exit(0);





}