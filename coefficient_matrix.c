#include<stdio.h>
#include<stdlib.h>
#include"type.h"
#include"coefficient_matrix.h"
#include"scalar.h"
#include"b_matrix.h"
#include"d_matrix.h"

extern Global global;
extern Option option;

void generate_coefficient_matrix(){
    double ke_matrix[60][60];
    double current_deformation_gradient[3][3];
    double current_stresses[6];
    double back_stress[6];
    double trial_elastic_strains[6];

    for(int point = 0; point < global.subdomain.N_point; point++){
        for(int i = 0; i < option.dim; i++)
            for(int j = 0; j < option.dim; j++)
                current_deformation_gradient[i][j] = global.subdomain.deformation_gradients[i][j][point];
    
        for(int i = 0; i < 6; i++){
            current_stresses[i] = global.subdomain.current_stresses[point][i];
            trial_elastic_strains[i] = global.subdomain.trial_elastic_strains[point][i];
            back_stress[i] = global.subdomain.back_stresses[point][i];
        }
    }

    for(int point = 0; point < global.subdomain.N_point; point++){
        generate_subdomain_coefficient_matrix(point, ke_matrix, current_deformation_gradient, current_stresses, trial_elastic_strains,
        global.subdomain.equivalent_plastic_strains, global.subdomain.equivalent_plastic_strain_increments, back_stress);
    }
}

void generate_subdomain_coefficient_matrix(int point_n, double (*ke_matrix)[60], 
                                            double (*current_deformation_gradients)[3], double *current_stress, double *trial_elastic_strains,
                                            double *equivalemt_plastic_strains, double *equivalent_plastic_strain_increments, double *back_stresses){
    double b_t_matrix[60][6];
    double d_matrix[6][6];
    double BTD[60][6];
    double jacobian;
    int N_support = global.subdomain.support_offset[point_n + 1] - global.subdomain.support_offset[point_n];

    //ke_matrixをゼロ処理
    for(int i = 0; i < option.dim * (N_support + 1); i++)
        for(int j = 0; j < option.dim * (N_support + 1); j++)
            ke_matrix[i][j] = 0.;

    //ヤコビアンの計算
    jacobian = calc_subdomain_volume(point_n);

    //bマトリクスの計算
    generate_linear_b_matrix(b_t_matrix, point_n);

    //弾性Dマトリクスの計算
    generateElasticDMatrix(d_matrix);

    //有限ひずみのDマトリクスに修正
    modify_d_matrix_with_finite_strain(d_matrix, current_stress, trial_elastic_strains, current_deformation_gradients);

    //弾性Bマトリクスの計算
    generate_linear_b_matrix(b_t_matrix, point_n);

    for(int i = 0; i < option.dim * (N_support + 1); i++){
        for(int j = 0; j < 6; j++){
            double BTD_ij = 0.;
            for(int k = 0; k < 6; k++){
                BTD_ij += b_t_matrix[i][k] * d_matrix[k][j];
            }
            BTD[i][j] = BTD_ij;
        }
    }

    for(int i = 0; i < option.dim * (N_support + 1); i++){
        for(int j = 0; j < option.dim * (N_support + 1); j++){
            double ke_ij = 0.;
            for(int k = 0; k < 6; k++){
                ke_ij += BTD[i][k] * b_t_matrix[j][k];
            }
            ke_matrix[i][j] = ke_ij * jacobian;
        }
    }

    FILE *fp_debug;                                 //デバッグ用のファイル
    fp_debug = fopen("debag.dat", "w");

    for(int i = 0; i < option.dim * (N_support + 1); i++){
        for(int j = 0; j < option.dim * (N_support + 1); j++){
            fprintf(fp_debug, "%+3.2e  ", ke_matrix[i][j]);
        }
        fprintf(fp_debug,"\n");
    }
    fclose(fp_debug);
    exit(0);




}