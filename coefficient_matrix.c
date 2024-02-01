#include<stdio.h>
#include<stdlib.h>
#include"type.h"
#include"coefficient_matrix.h"
#include"scalar.h"
#include"matrix.h"
#include"b_matrix.h"
#include"d_matrix.h"
#include"tensor.h"
#include"s_matrix.h"
#include"GetGaussPoints.h"

extern Global global;
extern Option option;

void generate_coefficient_matrix(){
    double ke_matrix[60][60];
    double current_deformation_gradient[3][3];
    double current_stresses[6];
    double back_stress[6];
    double trial_elastic_strains[6];

    //接線剛性マトリクスの領域積分の項を計算
    for(int point = 0; point < global.subdomain.N_point; point++){
        for(int i = 0; i < option.dim; i++)
            for(int j = 0; j < option.dim; j++)
                current_deformation_gradient[i][j] = global.subdomain.deformation_gradients[i][j][point];
    
        for(int i = 0; i < 6; i++){
            current_stresses[i] = global.subdomain.current_stresses[point][i];
            trial_elastic_strains[i] = global.subdomain.trial_elastic_strains[point][i];
            back_stress[i] = global.subdomain.back_stresses[point][i];
        }
        
        generate_subdomain_coefficient_matrix(point, ke_matrix, current_deformation_gradient, current_stresses, trial_elastic_strains,
        global.subdomain.equivalent_plastic_strains, global.subdomain.equivalent_plastic_strain_increments, back_stress);
        assemble_coefficient_matrix_matrix_domain(ke_matrix, global.subdomain.Global_K, point, point);
    }
    //ペナルティ項（安定化項以外）の項を計算
    for(int face = 0; face < global.subdomain.N_int_boundary; face++){
        for(int i = 0; i < option.dim; i++)
            for(int j = 0; j < option.dim; j++)
                current_deformation_gradient[i][j] = global.subdomain.deformation_gradients[i][j][global.subdomain.pair_point_ib[2 * face]];
    
        for(int i = 0; i < 6; i++){
            current_stresses[i] = global.subdomain.current_stresses[global.subdomain.pair_point_ib[2 * face]][i];
            trial_elastic_strains[i] = global.subdomain.trial_elastic_strains[global.subdomain.pair_point_ib[2 * face]][i];
            back_stress[i] = global.subdomain.back_stresses[global.subdomain.pair_point_ib[2 * face]][i];
        }

        generate_subdomain_coefficient_matrix_for_panaltyterm(global.subdomain.pair_point_ib[2*face],global.subdomain.pair_point_ib[2*face], face, ke_matrix,
                                                    current_deformation_gradient, current_stresses, trial_elastic_strains,
                                                    global.subdomain.equivalent_plastic_strains, global.subdomain.equivalent_plastic_strain_increments, back_stress, 0);
        assemble_coefficient_matrix_matrix_domain(ke_matrix, global.subdomain.Global_K, global.subdomain.pair_point_ib[2*face], global.subdomain.pair_point_ib[2*face]);
        generate_subdomain_coefficient_matrix_for_panaltyterm(global.subdomain.pair_point_ib[2*face + 1],global.subdomain.pair_point_ib[2*face], face, ke_matrix,
                                                    current_deformation_gradient, current_stresses, trial_elastic_strains,
                                                    global.subdomain.equivalent_plastic_strains, global.subdomain.equivalent_plastic_strain_increments, back_stress, 1);
        assemble_coefficient_matrix_matrix_domain(ke_matrix, global.subdomain.Global_K, global.subdomain.pair_point_ib[2 * face + 1], global.subdomain.pair_point_ib[2*face]);

        for(int i = 0; i < option.dim; i++)
            for(int j = 0; j < option.dim; j++)
                current_deformation_gradient[i][j] = global.subdomain.deformation_gradients[i][j][global.subdomain.pair_point_ib[2 * face + 1]];
    
        for(int i = 0; i < 6; i++){
            current_stresses[i] = global.subdomain.current_stresses[global.subdomain.pair_point_ib[2 * face + 1]][i];
            trial_elastic_strains[i] = global.subdomain.trial_elastic_strains[global.subdomain.pair_point_ib[2 * face + 1]][i];
            back_stress[i] = global.subdomain.back_stresses[global.subdomain.pair_point_ib[2 * face + 1]][i];
        }

        generate_subdomain_coefficient_matrix_for_panaltyterm(global.subdomain.pair_point_ib[2*face],global.subdomain.pair_point_ib[2*face + 1], face, ke_matrix,
                                                    current_deformation_gradient, current_stresses, trial_elastic_strains,
                                                    global.subdomain.equivalent_plastic_strains, global.subdomain.equivalent_plastic_strain_increments, back_stress, 1);
        assemble_coefficient_matrix_matrix_domain(ke_matrix, global.subdomain.Global_K, global.subdomain.pair_point_ib[2*face], global.subdomain.pair_point_ib[2*face+1]);
        generate_subdomain_coefficient_matrix_for_panaltyterm(global.subdomain.pair_point_ib[2*face + 1],global.subdomain.pair_point_ib[2*face+1], face, ke_matrix,
                                                    current_deformation_gradient, current_stresses, trial_elastic_strains,
                                                    global.subdomain.equivalent_plastic_strains, global.subdomain.equivalent_plastic_strain_increments, back_stress, 0);
        assemble_coefficient_matrix_matrix_domain(ke_matrix, global.subdomain.Global_K, global.subdomain.pair_point_ib[2 * face + 1], global.subdomain.pair_point_ib[2*face+1]);
    }

    for(int face = 0; face < global.subdomain.N_int_boundary; face++){
        generate_subdomain_coefficient_matrix_for_StabilizationTerm(global.subdomain.pair_point_ib[2 *face], global.subdomain.pair_point_ib[2 *face], face, ke_matrix, 0);
        assemble_coefficient_matrix_matrix_domain(ke_matrix, global.subdomain.Global_K, global.subdomain.pair_point_ib[2 *face], global.subdomain.pair_point_ib[2 *face]);
        generate_subdomain_coefficient_matrix_for_StabilizationTerm(global.subdomain.pair_point_ib[2 *face], global.subdomain.pair_point_ib[2 *face + 1], face, ke_matrix, 1);
        assemble_coefficient_matrix_matrix_domain(ke_matrix, global.subdomain.Global_K, global.subdomain.pair_point_ib[2 *face], global.subdomain.pair_point_ib[2 *face + 1]);
        generate_subdomain_coefficient_matrix_for_StabilizationTerm(global.subdomain.pair_point_ib[2 *face+ 1], global.subdomain.pair_point_ib[2 *face], face, ke_matrix, 1);
        assemble_coefficient_matrix_matrix_domain(ke_matrix, global.subdomain.Global_K, global.subdomain.pair_point_ib[2 *face+ 1], global.subdomain.pair_point_ib[2 *face]);
        generate_subdomain_coefficient_matrix_for_StabilizationTerm(global.subdomain.pair_point_ib[2 *face + 1], global.subdomain.pair_point_ib[2 *face + 1], face, ke_matrix, 0);
        assemble_coefficient_matrix_matrix_domain(ke_matrix, global.subdomain.Global_K, global.subdomain.pair_point_ib[2 *face + 1], global.subdomain.pair_point_ib[2 *face + 1]);
    }


    FILE *fp_debug;
    fp_debug = fopen("debag.dat","w");
    for(int i = 0; i < 3*global.subdomain.N_point;i++){
        for(int j = 0; j < 3 * global.subdomain.N_point; j++){
            fprintf(fp_debug,"%+4.3e  ", global.subdomain.Global_K[i*3*global.subdomain.N_point+j]);
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



/*
ここでの変形勾配テンソル等はpoint_n2のサブドメイン内で定義したもの
*/
void generate_subdomain_coefficient_matrix_for_panaltyterm(int point_n1, int point_n2, int face_n, double (*ke_matrix)[60], 
                                            double (*current_deformation_gradients)[3], double *current_stress, double *trial_elastic_strains,
                                            double *equivalemt_plastic_strains, double *equivalent_plastic_strain_increments, double *back_stresses, double flag){

    int N_qu = 1;
    double factor;
    double xyz[3];
    double X[27], w[27];
    double jacobian; 
    int face_node[4];
    double face_node_XYZ[4][3];
    int N1_support = global.subdomain.support_offset[point_n1 + 1] - global.subdomain.support_offset[point_n1];
    int N2_support = global.subdomain.support_offset[point_n2 + 1] - global.subdomain.support_offset[point_n2];
    double **G;
    double NT[60][3];
    double Ne_d[3][9];
    double d_matrix[6][6];
    double c_matrix[9][9];
    double *latest_point_XYZ;
    double *node_XYZ;
    double CG[9][60];
    double neCG[3][60];
    double sign;

    if(flag == 0){
        sign = 1.0;
    }
    else if(flag == 1){
        sign = -1.0;
    }

    //ke_matrixをゼロ処理
    for(int i = 0; i < option.dim * (N1_support + 1); i++)
        for(int j = 0; j < option.dim * (N2_support + 1); j++)
            ke_matrix[i][j] = 0.;

    if((latest_point_XYZ = (double *)calloc(option.dim * global.subdomain.N_point, sizeof(double))) == NULL){
        exit(-1);
    }
    if((node_XYZ = (double *)calloc(option.dim * global.subdomain.N_node, sizeof(double))) == NULL){
        printf("Error:node_XYZ's memory is not enough\n");
        exit(-1);
    }

    for(int i = 0; i < global.subdomain.N_point; i++){
        for(int j = 0; j < option.dim; j++){
            latest_point_XYZ[option.dim * i + j] = global.subdomain.point_XYZ[option.dim * i + j]
                    + global.subdomain.displacement[i][j]
                    + global.subdomain.displacement_increment[i][j];
        }
    }
    for(int i = 0; i < global.subdomain.N_node; i++){
        for(int j = 0; j < option.dim; j++){
            node_XYZ[option.dim * i + j] = global.subdomain.node_XYZ[option.dim * i + j]
                                        + global.subdomain.nodal_displacements[i][j]
                                        + global.subdomain.nodal_displacement_increments[i][j];
        }
    }
    for(int i = 0; i < 4; i++)
            face_node[i] = global.subdomain.node[global.subdomain.vertex_offset[global.subdomain.shared_face[face_n]] + i];
        
        for(int i = 0; i < 4; i++)
            for(int j = 0; j < option.dim; j++)
                face_node_XYZ[i][j] = node_XYZ[option.dim * face_node[i] + j];

    //Gマトリクスの計算
    G = matrix(9, option.dim * (N2_support + 1));
    calc_G(option.dim, point_n2, latest_point_XYZ, global.subdomain.support_offset, global.subdomain.support, G);

    //弾性Dマトリクスの計算
    generateElasticDMatrix(d_matrix);

    //有限ひずみのDマトリクスに修正
    modify_d_matrix_with_finite_strain_for_PenaltyTerm(c_matrix, d_matrix, current_stress, trial_elastic_strains, current_deformation_gradients);

    //ガウス点の座標と重み、ヤコビアンの計算
    Gauss_points_and_weighting_factors(N_qu, X, w);
    jacobian = calc_surface_area(face_n) / 4.0;

    //法線ベクトルの計算
    calc_Ne_diagonal(option.dim, global.subdomain.pair_point_ib[2 * face_n], global.subdomain.pair_point_ib[2 * face_n + 1], face_n, global.subdomain.vertex_offset, global.subdomain.node, node_XYZ, latest_point_XYZ, Ne_d);

    //物質表記→空間表記する際の面積変化率(J/sqrt(nBn))の分子 ,sqrt(nBn) の計算
    factor = calc_area_change_factor(global.subdomain.pair_point_ib[2 * face_n], global.subdomain.pair_point_ib[2 * face_n + 1], Ne_d);

    //ペナルティ項を計算
    for(int s = 0;  s < N_qu; s++){
        for(int t = 0; t < N_qu; t++){
            for(int i = 0; i < option.dim; i++)
                    xyz[i] = 0.25 * (1.0 - X[s]) * (1.0 - X[t]) * face_node_XYZ[0][i]
                            + 0.25 * (1.0 - X[s]) * (1.0 + X[t]) * face_node_XYZ[1][i]
                            + 0.25 * (1.0 + X[s]) * (1.0 + X[t]) * face_node_XYZ[2][i]
                            + 0.25 * (1.0 + X[s]) * (1.0 - X[t]) * face_node_XYZ[3][i];
            
            //形状関数の計算
            calc_shape(xyz, option.dim, point_n1, latest_point_XYZ, global.subdomain.support_offset, NT);

            //被積分関数の計算
            for(int i = 0; i < 9; i++){
                double CG_ij = 0.;
                for(int j = 0; j < option.dim * (N2_support + 1); j++){
                    for(int k = 0; k < 9; k++){
                        CG_ij += c_matrix[i][k] * G[k][j];
                    }
                    CG[i][j] = CG_ij;
                }
            }
            for(int i = 0; i < option.dim; i++){
                for(int j = 0; j < option.dim * (N2_support + 1); j++){
                    double neCG_ij = 0.;
                    for(int k = 0; k < 9; k++){
                        neCG_ij += Ne_d[i][k] * CG[k][j];
                    }
                    neCG[i][j] = neCG_ij;
                }
            }
            for(int i = 0; i < option.dim * (N1_support + 1); i++){
                for(int j = 0; j < option.dim * (N2_support + 1); j++){
                    double ke_ij = 0.;
                    for(int k = 0; k < option.dim; k++){
                        ke_ij += NT[i][k] * neCG[k][j];
                    }
                    ke_matrix[i][j] += 0.5 * ke_ij * factor * jacobian * w[s] * w[t] * sign;
                }
            }

        }
    }

    free_matrix(G);
    free(node_XYZ);
    free(latest_point_XYZ);
    
}

void generate_subdomain_coefficient_matrix_for_StabilizationTerm(int point_n1, int point_n2, int face_n, double (*ke_matrix)[60],int flag){
    int N_qu = 2;
    double N1T[60][3];
    double N2T[60][3];
    double xyz[3];
    double X[27], w[27];
    double Ne_d[3][9];
    double deformation_gradients[3][3];
    double sign;
    int N1_support = global.subdomain.support_offset[point_n1 + 1] - global.subdomain.support_offset[point_n1];
    int N2_support = global.subdomain.support_offset[point_n2 + 1] - global.subdomain.support_offset[point_n2];
    double *latest_point_XYZ;
    double *node_XYZ;
    int face_node[4];
    double face_node_XYZ[4][3];
    double he;                                                  //ポイント間の距離
    double eta = global.material.penalty;                       //ペナルティパラメータ
    double jacobian;
    double factor;
    double inverse_volume_change;

    //ke_matrixをゼロ処理
    for(int i = 0; i < option.dim * (N1_support + 1); i++)
        for(int j = 0; j < option.dim * (N2_support + 1); j++)
            ke_matrix[i][j] = 0.;

    if(flag == 0){
        sign = 1.0;
    }
    else if(flag == 1){
        sign = -1.0;
    }

    if((latest_point_XYZ = (double *)calloc(option.dim * global.subdomain.N_point, sizeof(double))) == NULL){
        exit(-1);
    }
    if((node_XYZ = (double *)calloc(option.dim * global.subdomain.N_node, sizeof(double))) == NULL){
        printf("Error:node_XYZ's memory is not enough\n");
        exit(-1);
    }

    for(int i = 0; i < global.subdomain.N_point; i++){
        for(int j = 0; j < option.dim; j++){
            latest_point_XYZ[option.dim * i + j] = global.subdomain.point_XYZ[option.dim * i + j]
                    + global.subdomain.displacement[i][j]
                    + global.subdomain.displacement_increment[i][j];
        }
    }
    for(int i = 0; i < global.subdomain.N_node; i++){
        for(int j = 0; j < option.dim; j++){
            node_XYZ[option.dim * i + j] = global.subdomain.node_XYZ[option.dim * i + j]
                                        + global.subdomain.nodal_displacements[i][j]
                                        + global.subdomain.nodal_displacement_increments[i][j];
        }
    }
    for(int i = 0; i < 4; i++)
            face_node[i] = global.subdomain.node[global.subdomain.vertex_offset[global.subdomain.shared_face[face_n]] + i];
        
    for(int i = 0; i < 4; i++)
        for(int j = 0; j < option.dim; j++)
            face_node_XYZ[i][j] = node_XYZ[option.dim * face_node[i] + j];

    he = distance(option.dim, global.subdomain.pair_point_ib[2 * face_n], global.subdomain.pair_point_ib[2 * face_n + 1], latest_point_XYZ);

    //ヤコビアンの計算
    jacobian = calc_surface_area(face_n) / 4.0;

    //体積変化率の計算
    for(int i = 0; i < option.dim; i++){
        for(int j = 0; j < option.dim; j++){
            deformation_gradients[i][j] 
                = 0.5 * (global.subdomain.current_deformation_gradients[i][j][global.subdomain.pair_point_ib[2 * face_n]] 
                        + global.subdomain.current_deformation_gradients[i][j][global.subdomain.pair_point_ib[2 * face_n + 1]]);
        }
    }
    inverse_volume_change = 1.0 / calc_3x3matrix_determinant(deformation_gradients);

    //物質表記→空間表記する際の面積変化率(J/sqrt(nBn))の分子 ,sqrt(nBn) の計算
    //法線ベクトルの計算
    calc_Ne_diagonal(option.dim, global.subdomain.pair_point_ib[2 * face_n], global.subdomain.pair_point_ib[2 * face_n + 1], face_n, global.subdomain.vertex_offset, global.subdomain.node, node_XYZ, latest_point_XYZ, Ne_d);
    factor = calc_area_change_factor(global.subdomain.pair_point_ib[2 * face_n], global.subdomain.pair_point_ib[2 * face_n + 1], Ne_d);
    
    for(int s = 0; s < N_qu; s++){
        for(int t = 0; t < N_qu; t++){
            for(int i = 0; i < option.dim; i++)
                    xyz[i] = 0.25 * (1.0 - X[s]) * (1.0 - X[t]) * face_node_XYZ[0][i]
                            + 0.25 * (1.0 - X[s]) * (1.0 + X[t]) * face_node_XYZ[1][i]
                            + 0.25 * (1.0 + X[s]) * (1.0 + X[t]) * face_node_XYZ[2][i]
                            + 0.25 * (1.0 + X[s]) * (1.0 - X[t]) * face_node_XYZ[3][i];
            
            calc_shape(xyz, option.dim, point_n1, latest_point_XYZ, global.subdomain.support_offset, N1T);
            calc_shape(xyz, option.dim, point_n2, latest_point_XYZ, global.subdomain.support_offset, N2T);

            for(int i = 0;  i < option.dim * (N1_support + 1); i++){
                for(int j = 0; j < option.dim * (N2_support + 1); j++){
                    double N1TN_ij = 0.;
                    for(int k = 0; k < option.dim; k++){
                        N1TN_ij += N1T[i][k] * N2T[j][k];
                    }
                    ke_matrix[i][j] += eta / he * N1TN_ij * factor * inverse_volume_change * sign * jacobian * w[s] * w[t];
                }
            }

        }
    }
    free(node_XYZ);
    free(latest_point_XYZ);

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

double calc_area_change_factor(int subdomain_n1, int subdomain_n2, double (*Ne_d)[9]){
    double trial_elastic_left_cauchy_green_deformations1[3][3];
    double trial_elastic_strain_tensor1[3][3];
    double trial_elastic_left_cauchy_green_deformations2[3][3];
    double trial_elastic_strain_tensor2[3][3];
    double trial_elastic_strains1[6];
    double trial_elastic_strains2[6];
    double trial_elastic_left_cauchy_green_deformations[3][3];
    double factor = 0.;
    double Bn[3];

    for(int i = 0; i < 6; i++){
        trial_elastic_strains1[i] = global.subdomain.trial_elastic_strains[subdomain_n1][i];
        trial_elastic_strains2[i] = global.subdomain.trial_elastic_strains[subdomain_n2][i];

    }
    
    //試行弾性左コーシーグリーンテンソルの計算
    trial_elastic_strain_tensor1[0][0] = 2.0 * trial_elastic_strains1[0];
    trial_elastic_strain_tensor1[0][1] = 2.0 * 0.5 * trial_elastic_strains1[3];
    trial_elastic_strain_tensor1[0][2] = 2.0 * 0.5 * trial_elastic_strains1[5];
    trial_elastic_strain_tensor1[1][0] = 2.0 * 0.5 * trial_elastic_strains1[3];
    trial_elastic_strain_tensor1[1][1] = 2.0 * trial_elastic_strains1[1];
    trial_elastic_strain_tensor1[1][2] = 2.0 * 0.5 * trial_elastic_strains1[4];
    trial_elastic_strain_tensor1[2][0] = 2.0 * 0.5 * trial_elastic_strains1[5];
    trial_elastic_strain_tensor1[2][1] = 2.0 * 0.5 * trial_elastic_strains1[4];
    trial_elastic_strain_tensor1[2][2] = 2.0 * trial_elastic_strains1[2];

    calculateTensorExponent(trial_elastic_left_cauchy_green_deformations1,
                            trial_elastic_strain_tensor1);

    //試行弾性左コーシーグリーンテンソルの計算
    trial_elastic_strain_tensor2[0][0] = 2.0 * trial_elastic_strains2[0];
    trial_elastic_strain_tensor2[0][1] = 2.0 * 0.5 * trial_elastic_strains2[3];
    trial_elastic_strain_tensor2[0][2] = 2.0 * 0.5 * trial_elastic_strains2[5];
    trial_elastic_strain_tensor2[1][0] = 2.0 * 0.5 * trial_elastic_strains2[3];
    trial_elastic_strain_tensor2[1][1] = 2.0 * trial_elastic_strains2[1];
    trial_elastic_strain_tensor2[1][2] = 2.0 * 0.5 * trial_elastic_strains2[4];
    trial_elastic_strain_tensor2[2][0] = 2.0 * 0.5 * trial_elastic_strains2[5];
    trial_elastic_strain_tensor2[2][1] = 2.0 * 0.5 * trial_elastic_strains2[4];
    trial_elastic_strain_tensor2[2][2] = 2.0 * trial_elastic_strains2[2];

    calculateTensorExponent(trial_elastic_left_cauchy_green_deformations2,
                            trial_elastic_strain_tensor2);
    
    for(int i = 0; i < option.dim; i++){
        for(int j = 0; j < option.dim; j++){
            trial_elastic_left_cauchy_green_deformations[i][j]
                = 0.5 * (trial_elastic_left_cauchy_green_deformations1[i][j] + trial_elastic_left_cauchy_green_deformations2[i][j]);
        }
    }

    for(int i = 0; i < option.dim; i++){
        double Bn_i =0.;
        for(int j = 0; j < option.dim; j++){
            Bn_i += trial_elastic_left_cauchy_green_deformations[i][j] * Ne_d[0][j];
        }
        Bn[i] = Bn_i;
    }
    for(int i = 0; i < option.dim; i++)
        factor += Ne_d[0][i] * Bn[i];

    return factor;
    

}