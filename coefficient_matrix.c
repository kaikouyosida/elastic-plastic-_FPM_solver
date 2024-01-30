#include<stdio.h>
#include<stdlib.h>
#include"type.h"
#include"coefficient_matrix.h"
#include"scalar.h"
#include"b_matrix.h"
#include"d_matrix.h"
#include"s_matrix.h"

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
        assemble_coefficient_matrix_matrix_domain(ke_matrix, global.subdomain.Global_K, point, point);
    }

    FILE *fp_debug;                                 //デバッグ用のファイル
    fp_debug = fopen("debag.dat", "w");

    for(int i = 0; i < option.dim * global.subdomain.N_point; i++){
        for(int j = 0; j < option.dim * global.subdomain.N_point; j++){
            fprintf(fp_debug, "%+3.2e  ", global.subdomain.Global_K[i*global.subdomain.N_point * option.dim + j]);
        }
        fprintf(fp_debug,"\n");
    }
    fclose(fp_debug);
    exit(0);
    
}

void generate_subdomain_coefficient_matrix(int point_n, double (*ke_matrix)[60], 
                                            double (*current_deformation_gradients)[3], double *current_stress, double *trial_elastic_strains,
                                            double *equivalemt_plastic_strains, double *equivalent_plastic_strain_increments, double *back_stresses){
    double b_t_matrix[60][6];
    double b_NL_matrix[60][9];
    double s_matrix[9][9];
    double d_matrix[6][6];
    double BTD[60][6];
    double GTS[60][9];
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

    //非線形Bマトリクスの計算
    generate_nonlinear_b_matrix(b_NL_matrix, point_n);

    //Sマトリクスの計算
    generateSMatrix(s_matrix, current_stress);

    for(int i = 0; i < option.dim * (N_support + 1); i++){
        for(int j = 0; j < 9; j++){
            double GTS_ij = 0.;
            for(int k = 0; k < 9; k++){
                GTS_ij += b_NL_matrix[i][k] * s_matrix[k][j];
            }
            GTS[i][j] = GTS_ij;
        }
    }
    
    for(int i = 0; i < option.dim * (N_support + 1); i++){
        for(int j = 0; j < option.dim * (N_support + 1); j++){
            double ke_ij = 0;
            for(int k = 0; k < option.dim * (N_support + 1); k++){
                ke_ij += GTS[i][k] * b_NL_matrix[j][k];
            }
            ke_matrix[i][j] += ke_ij * jacobian;
        }
    }

}
void assemble_coefficient_matrix_matrix_domain(double (*element_K)[60], double *Global_K, int point_n1, int point_n2){
    int ref_num1 = global.subdomain.support_offset[point_n1];
    int ref_num2 = global.subdomain.support_offset[point_n2];
    int N1_support = global.subdomain.support_offset[point_n1 + 1] - global.subdomain.support_offset[point_n1];
    int N2_support = global.subdomain.support_offset[point_n2 + 1] - global.subdomain.support_offset[point_n2];
    int DoF_free = option.dim * global.subdomain.N_point;
    
    for(int i = 0; i < option.dim; i++){
        for(int j = 0; j < option.dim; j++){
            Global_K[DoF_free * (option.dim * point_n1 + i) + option.dim * point_n2 + j]
                    += element_K[i][j];
        }
    }
    for(int i = 0; i < N2_support; i++){
        for(int j = 0; j < option.dim; j++){
            for(int k = 0; k < option.dim; k++){
                Global_K[DoF_free * (option.dim * point_n1 + j) + option.dim * global.subdomain.support[ref_num2 + i] + k]
                    += element_K[j][option.dim * (i + 1) + k];
            }
        }
    }
    for(int i = 0; i < N1_support; i++){
        for(int j = 0; j < option.dim; j++){
            for(int k = 0; k < option.dim; k++){
                Global_K[DoF_free * (option.dim * global.subdomain.support[ref_num1 + i] + j) + option.dim * point_n2 + k]
                    += element_K[option.dim * (i + 1) + j][k];
            }
        }
    }
    for(int i = 0; i < N1_support; i++){
        for(int j = 0; j < N2_support; j++){
            for(int k = 0; k < option.dim; k++){
                for(int l = 0; l < option.dim; l++){
                    Global_K[DoF_free * (option.dim * global.subdomain.support[ref_num1 + i] + k) + (option.dim * global.subdomain.support[ref_num2 + j] + l)]
                    += element_K[option.dim * (i + 1) + k][option.dim * (j + 1) + l];
                }
            }
        }
    }
    
}