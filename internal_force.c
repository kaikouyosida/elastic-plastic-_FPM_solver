#pragma warning(disable: 4100) // 引数が未使用の場合
#pragma warning(disable: 4189) // ローカル変数が未使用の場合
#pragma warning(disable: 4996) //fopenの警告番号
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"type.h"
#include"internal_force.h"
#include"vector.h"
#include"d_matrix.h"
#include"b_matrix.h"
#include"matrix.h"
#include"tensor.h"
#include"stress.h"
#include"ss_curve.h"
#include"scalar.h"
#include"GetGaussPoints.h"
#include"ImposeDirichretCondition.h"

extern Global global;
extern Option option;

FILE *fp_debug;
char FILE_name[128];

double temp;

#define NUMBER_OF_NODE_IN_FACE 4
#define NUMBER_OF_NODE_IN_SUBDOMAIN 8

void update_field_and_internal_forces(){
    long long DoF_Free = option.dim * global.subdomain.N_point; //自由度数
    double d_matrix[6][6];                                      //Dマトリクス
    double inverse_relative_deformation_gradient[3][3];         //相対変形勾配テンソルの逆テンソル
    double relative_deformation_gradient[3][3];                 //相対変形勾配テンソル
    double elastic_strain_tensor[3][3];                         //弾性ひずみテンソル
    double deformation_gradients[3][3];                         //変形勾配テンソル
    double current_deformation_gradients[3][3];                 //現配置の変形勾配テンソル
    double elastic_left_cauchy_green_deformations[3][3];        //弾性左コーシーグリーンテンソル
    double trial_elastic_left_cauchy_green_deformations[3][3];  //試行弾性左コーシーグリーンテンソル
    double trial_relative_stresses[6];                          //試行相対応力
    double b_t_matrix[180][6];                                   //Bマトリクスの転置
    double *current_point_XYZ;                                   //現配置のポイント配置
    double **G;                                                 //GFD法によって導かれるGマトリクス
    int support[180];                                             //サポートドメイン内のポイント数
    double displacement_increment[3];                           //サポートの変位増分
    double trial_relative_equivalent_stress;
    double factor;

    //internal_forceをゼロ処理
    for(int i = 0; i < global.subdomain.N_point; i++){
        for(int j = 0; j < option.dim; j++){
            global.subdomain.global_internal_force[i][j] = 0.;
        }
    }
    
    //各サブドメインで内力ベクトルを作成
    for(int point = 0; point < global.subdomain.N_point; point++){
        int N_support = global.subdomain.support_offset[point + 1] - global.subdomain.support_offset[point];
               
        for(int i = 0 ; i < N_support; i++)
            support[i] = global.subdomain.support[global.subdomain.support_offset[point] + i];
        
        for(int i = 0; i < option.dim; i++){
            for(int j = 0; j < option.dim; j++){
                deformation_gradients[i][j] = global.subdomain.deformation_gradients[i][j][point];
                current_deformation_gradients[i][j] = global.subdomain.current_deformation_gradients[i][j][point];           
            }
        }
        
        double *kinematic_hardening_fractions = &global.material.kinematic_hardening_fractions;
        double *equivalent_plastic_strain = &global.subdomain.equivalent_plastic_strains[point];
        double *equivalent_plastic_strain_increment = &global.subdomain.equivalent_plastic_strain_increments[point];
        double *yield_stress = &global.subdomain.yield_stresses[point];
        double *current_yield_stress = &global.subdomain.current_yield_stresses[point];

        double *elastic_strains = global.subdomain.elastic_strains[point];
        double *current_elastic_strains = global.subdomain.current_elastic_strains[point];
        double *trial_elastic_strains = global.subdomain.trial_elastic_strains[point];
        double *current_stresses = global.subdomain.current_stresses[point];
        double *back_stresses = global.subdomain.back_stresses[point];
        double *current_back_stresses = global.subdomain.current_back_stresses[point];

        //Bマトリクスを作成
        generate_linear_b_matrix(b_t_matrix, point);
        
        //dマトリクスを作成
        generateElasticDMatrix(d_matrix);
        
        if((current_point_XYZ = (double *)calloc(DoF_Free, sizeof(double))) == NULL){
            printf("Error: Latest_point_XYZ's memory is not enough\n");
            exit(-1);
        }
        for(int i = 0; i < global.subdomain.N_point; i++){
            for(int j = 0; j < option.dim; j++){
                current_point_XYZ[option.dim * i + j] = global.subdomain.point_XYZ[option.dim * i + j]
                                                    + global.subdomain.displacement[i][j]
                                                    + global.subdomain.displacement_increment[i][j];
            }
        }

        //変位勾配のマトリクスを作成（動的に確保）
        G = matrix(option.dim * option.dim, option.dim * (N_support + 1));
        calc_G(option.dim, point, current_point_XYZ, global.subdomain.support_offset, global.subdomain.support, G);
        
        
        //現配置の変形勾配テンソルの逆行列を計算
        //相対変形勾配テンソルを計算。dF=(I-d(du)/d(x+u+du))^-1//
        identify3x3Matrix(inverse_relative_deformation_gradient);

        for(int i = 0; i < option.dim; i++){
            for(int j = 0; j < option.dim; j++){
                for(int k = 0; k < N_support; k++){
                    inverse_relative_deformation_gradient[i][j] -= G[option.dim * i + j][option.dim * (k + 1) + i] * global.subdomain.displacement_increment[support[k]][i];
                }
            }
            for(int j = 0; j < option.dim; j++){
                inverse_relative_deformation_gradient[i][j] -= G[option.dim * i + j][i] * global.subdomain.displacement_increment[point][i];
            }
        }
        free_matrix(G);
        free(current_point_XYZ);

        //変形勾配テンソル（初期配置に対してのテンソル）[F] = [dF] * [F]
        invert3x3Matrix(relative_deformation_gradient, inverse_relative_deformation_gradient);
        for(int i = 0; i < option.dim; i++)
            for (int j = 0; j < option.dim; j++){

                double deformation_gradient_i_j = 0.0;

                for (int k = 0; k < option.dim; k++)
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

        //試行弾性左コーシーグリーンテンソルの計算 [B]^trial = [dF] * [B]^e * [dF]^T 
        for (int i = 0; i < option.dim; i++)
            for (int j = 0; j < option.dim; j++)
            {
                double left_cauchy_green_deformation_i_j = 0.0;

                for (int k = 0; k < option.dim; k++)
                {
                    double left_cauchy_green_deformation_k_l_times_relative_deformation_gradient_j_l = 0.0;

                    for (int l = 0; l < option.dim; l++)
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

        //試行対数ひずみへ受け渡す
        for (int i = 0; i < 6; i++)
            trial_elastic_strains[i] = current_elastic_strains[i];        
        
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
                
                double hardening_stress_increment;
                double current_relative_hydrostatic_stress;
                const double three_times_shear_modulus
                                = 3.0 * 0.5 * global.material.E_mod / (1.0 + global.material.nu_mod);

                //相当塑性ひずみ増分の計算
                *equivalent_plastic_strain_increment
                    = calc_equivalent_plastic_strain_increment(trial_relative_equivalent_stress,
                                                                *equivalent_plastic_strain,
                                                                *yield_stress);
                                        
                //硬化応力増分の計算
                hardening_stress_increment
                    = get_hardening_stress((*equivalent_plastic_strain) + (*equivalent_plastic_strain_increment))
                    - get_hardening_stress(*equivalent_plastic_strain);
               
                //試行相対偏差応力の計算
                current_relative_hydrostatic_stress
                    = (1.0 / 3.0)
                    * (trial_relative_stresses[0]
                       + trial_relative_stresses[1]
                       + trial_relative_stresses[2]);
                trial_relative_stresses[0] -= current_relative_hydrostatic_stress;
                trial_relative_stresses[1] -= current_relative_hydrostatic_stress;
                trial_relative_stresses[2] -= current_relative_hydrostatic_stress;


                //最終的な対数ひずみの計算
                factor
                    = (*equivalent_plastic_strain_increment)
                    * 1.5
                    / trial_relative_equivalent_stress;
                for (int i = 0; i < 3; i++)
                    current_elastic_strains[i]
                        -= factor * trial_relative_stresses[i];
                for (int i = 3; i < 6; i++)
                    current_elastic_strains[i]
                        -= 2.0 * factor * trial_relative_stresses[i];
                
                //最終的な降伏応力
                *current_yield_stress
                    = (*yield_stress)
                    + (1.0 - *kinematic_hardening_fractions)
                    * hardening_stress_increment;
                //最終的な背応力の計算
                factor
                    = *kinematic_hardening_fractions
                    * hardening_stress_increment
                    / trial_relative_equivalent_stress;
                for (int i = 0; i < 6; i++)
                    current_back_stresses[i]
                        = back_stresses[i]
                        + factor * trial_relative_stresses[i];   

                //最終的な応力の計算
                factor
                    = (*current_yield_stress)
                    / trial_relative_equivalent_stress;
                for (int i = 0; i < 6; i++)
                    current_stresses[i]
                        = current_back_stresses[i]
                        + factor * trial_relative_stresses[i];

                current_stresses[0] += current_relative_hydrostatic_stress;
                current_stresses[1] += current_relative_hydrostatic_stress;
                current_stresses[2] += current_relative_hydrostatic_stress; 

                //降伏曲面を更新            
                *yield_stress = *current_yield_stress;
            }

            //Kirchhoff応力からCauchy応力の計算
            const double inverse_volume_change
                = 1.0 / calc_3x3matrix_determinant(current_deformation_gradients);

            for (int i = 0; i < 6; i++)
                current_stresses[i] *= inverse_volume_change;

            for(int i = 0; i < option.dim; i++){
                for(int j = 0; j < option.dim; j++){
                    global.subdomain.current_deformation_gradients[i][j][point] = current_deformation_gradients[i][j];
                }
            }

    }
    
    //内力ベクトルの体積積分項を計算
    calc_internal_force_volume(global.subdomain.current_stresses);

    //内力ベクトルのペナルティ項を計算
    calc_internal_force_penalty(global.subdomain.current_stresses);
   
    //内力ベクトルの安定化項を計算
    calc_internal_force_penalty_stabilization();
}

double calc_equivalent_plastic_strain_increment(const double trial_relative_equivalent_stress,
                                                const double equivalent_plastic_strain,
                                                const double yield_stress){
        const double tolerance = 1.0e-8 * yield_stress; /* Piecewise linear SS curve should give the exact solution */
        const int max_iteration_count = 100;

        const double three_times_shear_modulus
        = 3.0 * 0.5 * global.material.E_mod / (1.0 + global.material.nu_mod);

        const double hardening_stress = get_hardening_stress(equivalent_plastic_strain);

        double current_hardening_stress, current_hardening_modulus;
        double residual, residual_gradient;

        double dep = 0.0;

        //Solve 3 * G * dep - (se - sy - (sh (ep + dep) - sh (ep))) = 0の求解
    for (int i = 0; i < max_iteration_count; i++){
        current_hardening_stress
            = get_hardening_stress(equivalent_plastic_strain + dep);

        residual
            = three_times_shear_modulus * dep
            - (trial_relative_equivalent_stress - yield_stress
            - (current_hardening_stress - hardening_stress));
        
        if (fabs(residual) <= tolerance)
            return dep;

        current_hardening_modulus
            = get_hardening_modulus(equivalent_plastic_strain + dep);

        residual_gradient = three_times_shear_modulus + current_hardening_modulus;

        dep -= residual / residual_gradient;
    }

#if 0
    printError("Warning: equivalent plastic strain increment calculation iteration not converged\n");
#endif

    return nan("");
}


void calc_internal_force_volume(double **current_stress){
    double subdomain_internal_force[180];            //サブドメイン単位での内力ベクトル
    double b_t_matrix[180][6];                       //bマトリクス
    double volume;                                  //サブドメインの体積
    
    for(int point = 0; point < global.subdomain.N_point; point++){
        int N_support = global.subdomain.support_offset[point + 1] - global.subdomain.support_offset[point];

        //各サブドメインの体積を計算
        if(option.solver_type == 1){
            volume = calc_subdomain_volume(point);
        }else{
            volume = calc_initial_subdomain_volume(point);
        }
        
        //Bマトリクスを作成
        generate_linear_b_matrix(b_t_matrix, point);

        //内力ベクトル一項目を計算（[B]^T * {sigma} * dV）
        for(int i = 0; i < N_support + 1; i++){
            for(int j = 0; j < option.dim; j++){
                double force_j = 0.;
                for(int k = 0; k < 6; k++)
                    force_j += b_t_matrix[option.dim * i + j][k] * current_stress[point][k];
                
                subdomain_internal_force[option.dim * i + j] = force_j * volume;
            }
        }
        //各サブドメインで内力ベクトルをアセンブル
        assemble_vector(point, global.subdomain.global_internal_force, subdomain_internal_force);
    }
}

void calc_internal_force_penalty(double **all_stress)
{
    FILE *fp_debug_fint;            
    const int N_qu = 1;
    double sign;                                                //項の符号
    double X[27], w[27];                                        //ガウス求積に使う正規化座標と重み関数
    double xyz[3];                                              //求積点の座標
    double Ne[3][6];                                            //内部境界の法線ベクトル
    double NT[180][3];                                           //形状関数の転置
    double subdomain_internal_force[180];                        //サブドメインごとの内力ベクトル
    double face_node_XYZ[NUMBER_OF_NODE_IN_FACE][3];            //面を構成する節点の座標
    int subdomain_node[NUMBER_OF_NODE_IN_SUBDOMAIN];            //サブドメイン中の節点番号
    int node_id[NUMBER_OF_NODE_IN_FACE];                        //サブドメインの節点アドレス
    double mapping_parameter;                                   //正規化座標→物理座標へ変換するための因子
    double *current_point_XYZ;                                  //ポイントの現在の座標
    

    //形状関数を計算するためのpoint現在座標を計算
    if((current_point_XYZ = (double *)calloc(option.dim * global.subdomain.N_point, sizeof(double))) == NULL){
        printf("Error:point_XYZ's memory is not enough\n");
        exit(-1);
    }
    if(option.solver_type == 1){
       for(int i = 0; i < global.subdomain.N_point; i++)
            for(int j = 0; j < option.dim; j++)
                current_point_XYZ[option.dim * i + j] = global.subdomain.point_XYZ[option.dim * i + j]
                                                        + global.subdomain.displacement[i][j]
                                                        + global.subdomain.displacement_increment[i][j]; 
    }else{
        for(int i = 0; i < global.subdomain.N_point; i++){
            for(int j = 0; j < option.dim; j++){
                current_point_XYZ[option.dim * i + j] = global.subdomain.point_XYZ[option.dim * i + j];
            }
        }
    }
    
    //gaussポイントの座標と重みを計算
    Gauss_points_and_weighting_factors(N_qu, X, w);
    
    for(int face = 0; face < global.subdomain.N_int_boundary; face++){
        
        for(int i = 0; i < 2; i++){
            //ポイントiからポイントjへ延びる面shared_face[face]の法線ベクトルを計算 
            generate_unit_vec_to_mat3x6(global.subdomain.shared_face[face], global.subdomain.pair_point_ib[2 * face + i], global.subdomain.pair_point_ib[2 * face + (i + 1) % 2], current_point_XYZ, Ne);
            
            for(int j = 0; j < 2; j++){
                int N_support = global.subdomain.support_offset[global.subdomain.pair_point_ib[2 * face + j] + 1] - global.subdomain.support_offset[global.subdomain.pair_point_ib[2 * face + j]];

                //ペナルティ項の符号
                if((i + j) % 2 == 0){
                    sign = -1.0;
                }else{
                    sign = 1.0;
                }
                
                //両サブドメインにおける形状関数から得た節点の現在座標
                if(option.solver_type == 1){
                    generate_current_node_of_face(face_node_XYZ, global.subdomain.shared_face[face], global.subdomain.pair_point_ib[2 * face + i]);
                }else{
                    for(int a = 0; a < NUMBER_OF_NODE_IN_FACE; a++){
                        for(int b = 0; b < option.dim; b++){
                            face_node_XYZ[a][b] 
                                = global.subdomain.node_XYZ[option.dim * global.subdomain.node[global.subdomain.vertex_offset[global.subdomain.shared_face[face]] + a] + b];
                        }
                    }
                }
                
                
                //各サブドメインで被積分項を計算．サブドメイン単位の内力ベクトルを計算
                for(int s = 0; s < N_qu; s++){
                    for(int t = 0; t < N_qu; t++){
                        
                        //ガウスポイントの重みを計算
                        generate_gauss_point_coordinate(s, t, face_node_XYZ, X, xyz);

                        //サブドメインiの形状関数を計算
                        calc_shape(xyz, option.dim, global.subdomain.pair_point_ib[2 * face + j], current_point_XYZ, global.subdomain.support_offset, NT);
                        
                        //物理空間座標→正規化座標に変換するためのスカラー値を計算
                        mapping_parameter = calc_mapping_parameter(global.subdomain.shared_face[face], global.subdomain.pair_point_ib[2 * face + i], s, t, X);
                      
                        //被積分項を計算 {N^tNeσ}
                        for(int k = 0; k < option.dim * (N_support + 1); k++){
                            double subdomain_internal_force_k = 0.;
                            for(int l = 0; l < 6; l++){
                                double NTne_kl = 0.;
                                for(int m = 0; m < option.dim; m++){
                                    NTne_kl += NT[k][m] * Ne[m][l];
                                }
                                subdomain_internal_force_k += NTne_kl * all_stress[global.subdomain.pair_point_ib[2 * face + i]][l];
                            }

                            subdomain_internal_force[k] =  sign *  0.5 * subdomain_internal_force_k * mapping_parameter * w[s] * w[t];
                        }
                        
                        //サブドメイン単位での内力ベクトルを全体内力ベクトルにアセンブル
                        assemble_vector(global.subdomain.pair_point_ib[2 * face + j], global.subdomain.global_internal_force, subdomain_internal_force);
                          
                    }
                }                
            }           
        }    
    }
    free(current_point_XYZ);

}

void calc_internal_force_penalty_stabilization(){
    const int N_qu = 2;                                         //積分点数
    double sign;                                                //項の符号
    double X[27], w[27];                                        //ガウス求積に使う正規化座標と重み関数
    double subdomain_node[NUMBER_OF_NODE_IN_SUBDOMAIN];         //サブドメインのノード番号
    double ndoe_id[NUMBER_OF_NODE_IN_FACE];                     //ノード番号のアドレス
    double xyz[3];                                              //求積点の座標
    double u_h[3];                                              //試行関数
    double NT[180][3];                                           //形状関数の転置
    double subdomain_internal_force[180];                        //サブドメインごとの内力ベクトル
    double face_node_XYZ[4][3];                                 //面を構成する節点の座標
    double face_node_XYZ1[4][3];
    double face_node_XYZ2[4][3];
    double mapping_parameter;                                   //物理座標→正規化座標へのマッピングに要するパラメータ
    double *current_point_XYZ;                                  //ポイントの現在の座標
    double he;                                                  //ポイント間の距離
    double eta = global.material.penalty;                 //ペナルティパラメータ
    
    //形状関数を計算するためのpoint, nodeにおける現配置の座標を計算
    if((current_point_XYZ = (double *)calloc(option.dim * global.subdomain.N_point, sizeof(double))) == NULL){
        printf("Error:point_XYZ's memory is not enough\n");
        exit(-1);
    }
    if(option.solver_type == 1){
      for(int i = 0; i < global.subdomain.N_point; i++){
            for(int j = 0; j < option.dim; j++){
                current_point_XYZ[option.dim * i + j] = global.subdomain.point_XYZ[option.dim * i + j];
                                                    + global.subdomain.displacement[i][j]
                                                    + global.subdomain.displacement_increment[i][j];
            }
        }  
    }else{
        for(int i = 0; i < global.subdomain.N_point; i++){
            for(int j = 0; j < option.dim; j++){
                current_point_XYZ[option.dim * i + j] = global.subdomain.point_XYZ[option.dim * i + j];
            }
        }
    }

    //gauss積分点の座標を計算
    Gauss_points_and_weighting_factors(N_qu, X, w);
    
    for(int face = 0; face < global.subdomain.N_int_boundary; face++){
        //ポイント間の距離を計算
        he = distance(option.dim, global.subdomain.pair_point_ib[2 * face], global.subdomain.pair_point_ib[2 * face + 1], current_point_XYZ);

        //Γ*の頂点の座標を計算（Γ+とΓ-の平均を計算）
        if(option.solver_type == 1){
            generate_current_node_of_face(face_node_XYZ1, global.subdomain.shared_face[face], global.subdomain.pair_point_ib[2 * face]);
            generate_current_node_of_face(face_node_XYZ2, global.subdomain.shared_face[face], global.subdomain.pair_point_ib[2 * face + 1]);
            for(int i = 0; i < NUMBER_OF_NODE_IN_FACE; i++){
                for(int j = 0; j < option.dim; j++){
                    face_node_XYZ[i][j] = 0.5 * (face_node_XYZ1[i][j] + face_node_XYZ2[i][j]);
                }
            }
        }else{
            for(int i = 0; i < NUMBER_OF_NODE_IN_FACE; i++){
                for(int j = 0; j < option.dim; j++){
                    face_node_XYZ[i][j] 
                        = global.subdomain.node_XYZ[option.dim * global.subdomain.node[global.subdomain.vertex_offset[global.subdomain.shared_face[face]] + i] + j];
                }
            }
        }
        

        for(int i = 0; i < 2; i++){
            for(int j = 0; j < 2; j++){
                int N_support = global.subdomain.support_offset[global.subdomain.pair_point_ib[2 * face + j] + 1] - global.subdomain.support_offset[global.subdomain.pair_point_ib[2 * face + j]];
                
                //符号を計算
                if((i + j) % 2 == 0){
                    sign = 1.0;
                }else{
                    sign = -1.0;
                }

                for(int s = 0; s < N_qu; s++){
                    for(int t = 0; t < N_qu; t++){
                        //物理空間座標→正規化座標に変換するためのスカラー値を計算(Γ*の面積変化率を計算する)
                        mapping_parameter = calc_mapping_parameter_for_av_area(face_node_XYZ, s, t, X);

                        //物理座標におけるガウス点の座標を計算
                        generate_gauss_point_coordinate(s, t, face_node_XYZ, X, xyz);

                        //サブドメインpair_point_ib[2 * face + j]の形状関数を計算
                        calc_shape(xyz, option.dim, global.subdomain.pair_point_ib[2 * face + j], current_point_XYZ, global.subdomain.support_offset, NT);

                        //サブドメインpair_point_ib[2 * face + i]の試行関数を計算
                        trial_u(xyz, global.subdomain.pair_point_ib[2 * face + i], current_point_XYZ, u_h, 1);

                        //要素内力ベクトルの計算
                        for(int k = 0; k < option.dim * (N_support + 1); k++){
                            double subdomain_internal_force_k = 0.;
                            for(int l = 0; l < option.dim; l++){
                                subdomain_internal_force_k += NT[k][l] * u_h[l];
                            }
                            subdomain_internal_force[k] = sign * eta / he * subdomain_internal_force_k * mapping_parameter * w[s] * w[t];
                        }
                        //全体剛性マトリクスにアセンブリ
                        assemble_vector(global.subdomain.pair_point_ib[2 * face + j], global.subdomain.global_internal_force, subdomain_internal_force);
                    }
                }
            }
        }

    }
    
    free(current_point_XYZ);
}

double calc_global_force_residual_norm(int iteration_step){
    double global_f_norm = 0., global_r_norm = 0.;

    for(int i = 0; i < global.subdomain.N_point; i++)
        for(int j = 0; j < option.dim; j++)
            global.subdomain.global_residual_force[i][j]
                 =  global.subdomain.global_external_force[i][j] - global.subdomain.global_internal_force[i][j];
    
    ImposeDirichretResidual(iteration_step);
  
    //外力ベクトルと残差ベクトルのノルムを計算
    global_f_norm = norm_for_mat(global.subdomain.global_external_force, global.subdomain.N_point, option.dim);
    global_r_norm = norm_for_mat(global.subdomain.global_residual_force, global.subdomain.N_point, option.dim);
    
    option.r_abso_norm = global_r_norm;
    
    if(iteration_step == 0)
        temp = global_r_norm;

    printf("%+15.14e\n", temp);
    if(fabs(global_f_norm) <= 1.0e-10){
        return global_r_norm / temp;
    }else{
        return global_r_norm / global_f_norm;
    }
}

//微小変形弾塑性解析における増分
void update_field_and_internal_infinitesimal(){
    int  support[60];                   //サブドメインのサポートドメイン
    double b_t_matrix[180][6];           //Bマトリクス
    double d_matrix[6][6];              //Dマトリクス
    double trial_relative_stresses[6];  //試行相対応力
    
    //全体内力ベクトルを初期化
    for(int i = 0; i < global.subdomain.N_point; i++){
        for(int j = 0; j < option.dim; j++){
            global.subdomain.global_internal_force[i][j] = 0.;
        }
    }
    
    //各サブドメインの内力ベクトルを計算
    for(int point = 0; point < global.subdomain.N_point; point++){
        int N_support = global.subdomain.support_offset[point + 1] - global.subdomain.support_offset[point];

        for(int i = 0 ; i < N_support; i++)
            support[i] = global.subdomain.support[global.subdomain.support_offset[point] + i];

        double *trial_elastic_strain = global.subdomain.trial_elastic_strains[point];
        double *elastic_strains = global.subdomain.elastic_strains[point];
        double *current_elastic_strains = global.subdomain.current_elastic_strains[point];
        double *current_stresses = global.subdomain.current_stresses[point];
        double *back_stresses = global.subdomain.back_stresses[point];
        double *current_back_stresses = global.subdomain.current_back_stresses[point];

        double *yield_stress = &global.subdomain.yield_stresses[point];
        double *equivalent_plastic_strain = &global.subdomain.equivalent_plastic_strains[point];
                
        double *equivalent_plastic_strain_increment = &global.subdomain.equivalent_plastic_strain_increments[point];
        double *current_yield_stress = &global.subdomain.current_yield_stresses[point];

        double trial_relative_equivalent_stress;
        double factor;
    
        //bマトリクスの計算
        generate_linear_b_matrix(b_t_matrix, point);
     
        //dマトリクスの計算
        generateElasticDMatrix(d_matrix);
        
        for(int i = 0; i < 6; i++)
            current_elastic_strains[i] = elastic_strains[i];
        
        //ひずみの更新　{epsilon}^trial = {epsilon} + [B] * {du}
        for(int i = 0; i < 6; i++){
            double current_elastic_strains_i = 0.;
            for(int j = 0; j < option.dim; j++){
                current_elastic_strains_i += b_t_matrix[j][i] * global.subdomain.displacement_increment[point][j];
            }
            for(int j = 0; j < N_support; j++){
                for(int k = 0; k < option.dim; k++){
                    current_elastic_strains_i += b_t_matrix[option.dim * (j + 1) + k][i] * global.subdomain.displacement_increment[support[j]][k];
                }
            }
            current_elastic_strains[i] += current_elastic_strains_i;
        }

        for(int i = 0; i < 6; i++)
            trial_elastic_strain[i] = current_elastic_strains[i];

        //試行応力の計算　{sigma}^trial = [D] * {epsilon}^trial
        for(int i = 0; i < 6; i++){
            double stress_i = 0.;
            for(int j = 0; j < 6; j++)
                stress_i += d_matrix[i][j] * current_elastic_strains[j];
            current_stresses[i] = stress_i;
        }


        //試行相対応力の計算
        for(int i = 0; i < 6; i++)
            trial_relative_stresses[i] = current_stresses[i] - back_stresses[i];
        
        //von-mises応力の計算
        trial_relative_equivalent_stress = calc_equivalent_stress(trial_relative_stresses);

        if (trial_relative_equivalent_stress <= (*yield_stress)){
            *equivalent_plastic_strain_increment = 0.0;
            *current_yield_stress = *yield_stress;
            for (int i = 0; i < 6; i++)
                current_back_stresses[i] = back_stresses[i];
        }else{
            double hardening_stress_increment;
            double current_relative_hydrostatic_stress;

            //相当塑性ひずみ増分の計算;
            *equivalent_plastic_strain_increment
                = calc_equivalent_plastic_strain_increment(trial_relative_equivalent_stress, *equivalent_plastic_strain, *yield_stress);
            
            //硬化応力の更新
            hardening_stress_increment = get_hardening_stress((*equivalent_plastic_strain) + (*equivalent_plastic_strain_increment))
                                        - get_hardening_stress(*equivalent_plastic_strain);
            
            //静水圧と偏差応力を計算
            current_relative_hydrostatic_stress
                = (1.0 / 3.0)
                * (trial_relative_stresses[0]
                +  trial_relative_stresses[1]
                +  trial_relative_stresses[2]);
            trial_relative_stresses[0] -= current_relative_hydrostatic_stress;
            trial_relative_stresses[1] -= current_relative_hydrostatic_stress;
            trial_relative_stresses[2] -= current_relative_hydrostatic_stress;

            //最終的な弾性ひずみの計算（下段の×２はvoigt表記）
            factor
                = (*equivalent_plastic_strain_increment)
                * 1.5
                / trial_relative_equivalent_stress;
            for (int i = 0; i < 3; i++)
                current_elastic_strains[i]
                    -= factor * trial_relative_stresses[i];
            for (int i = 3; i < 6; i++)
                current_elastic_strains[i]
                    -= 2.0 * factor * trial_relative_stresses[i];

             //最終的な降伏応力を計算
            *current_yield_stress
                = (*yield_stress)
                + (1.0 - global.material.kinematic_hardening_fractions)
                * hardening_stress_increment;

            //最終的な背応力を計算
            factor
                = global.material.kinematic_hardening_fractions
                * hardening_stress_increment
                / trial_relative_equivalent_stress;
            for (int i = 0; i < 6; i++)
                current_back_stresses[i]
                    = back_stresses[i]
                    + factor * trial_relative_stresses[i];

            //最終的な応力を計算
            factor
                = (*current_yield_stress)
                / trial_relative_equivalent_stress;
            for (int i = 0; i < 6; i++)
                current_stresses[i]
                    = current_back_stresses[i]
                    + factor * trial_relative_stresses[i];
            current_stresses[0] += current_relative_hydrostatic_stress;
            current_stresses[1] += current_relative_hydrostatic_stress;
            current_stresses[2] += current_relative_hydrostatic_stress;
            
            *yield_stress = *current_yield_stress;
        }
    }
    
   
    //体積積分項の内力ベクトルの計算
    calc_internal_force_volume(global.subdomain.current_stresses);
   
    //penalty項の内力ベクトルの計算
    calc_internal_force_penalty(global.subdomain.current_stresses);
    
    //安定化項の内力ベクトルの計算
    calc_internal_force_penalty_stabilization();
}

void update_point_displaecment_increment(double *du){
    int count = 0;
    for(int i = 0; i < global.subdomain.N_point; i++){
        for(int j = 0; j < option.dim; j++){
            int flag = 0; 
            for(int k = 0; k < global.bc.N_D_DoF; k++)
                if(option.dim * i + j == global.bc.fixed_dof[k]) flag = 1;
            if(flag == 0){
                global.subdomain.displacement_increment[i][j] += du[count];
                count++;
            }
        }
    }
}

void update_nodal_displacement_by_inital_NT(double *Initial_point_xyz){
    double node_xyz[3];
    double u_h[3];

    //u_hに0をfill-in
    for(int i = 0; i < option.dim; i++)
        u_h[i] = 0.;

    for(int point = 0; point < global.subdomain.N_point; point++){
        for(int i = 0; i < NUMBER_OF_NODE_IN_SUBDOMAIN; i++){
            
            for(int j = 0; j < option.dim; j++)
                node_xyz[j] = global.subdomain.node_XYZ[option.dim * global.subdomain.subdomain_node[NUMBER_OF_NODE_IN_SUBDOMAIN * point + i] + j]
                            + global.subdomain.nodal_displacement_sd[point][i][j];
        
            trial_u(node_xyz, point, Initial_point_xyz, u_h, 0);
            
            for(int j = 0; j < option.dim; j++)
                global.subdomain.nodal_displacement_sd[point][i][j] = u_h[j];
        }
    }
}

void update_nodal_displacement_increment(double *current_point_XYZ, double *Initial_point_xyz){
    double node_xyz[3];
    double u_h[3];
    double delta_u[3];

    //0をfill-in
    for(int i = 0; i < option.dim; i++){
        u_h[i] = 0.; delta_u[i] = 0.;
    }

    for(int point = 0; point < global.subdomain.N_point; point++){
        if(global.subdomain.equivalent_plastic_strain_increments[point] == 0.){
            for(int i = 0; i < NUMBER_OF_NODE_IN_SUBDOMAIN; i++){
                
                for(int j = 0; j < option.dim; j++)
                    node_xyz[j] = global.subdomain.node_XYZ[option.dim * global.subdomain.subdomain_node[NUMBER_OF_NODE_IN_SUBDOMAIN * point + i] + j]
                                + global.subdomain.nodal_displacement_sd[point][i][j];
            
                trial_u(node_xyz, point, Initial_point_xyz, u_h, 0);
                
                for(int j = 0; j < option.dim; j++)
                    global.subdomain.nodal_displacement_sd[point][i][j] = u_h[j];
            }
        }else{
            for(int i = 0; i < NUMBER_OF_NODE_IN_SUBDOMAIN; i++){
                for(int j = 0; j < option.dim; j++)
                    node_xyz[i] = global.subdomain.node_XYZ[option.dim * global.subdomain.subdomain_node[NUMBER_OF_NODE_IN_SUBDOMAIN * point + i] + j]
                                + global.subdomain.nodal_displacement_sd[point][i][j];
                
                trial_u(node_xyz, point, current_point_XYZ, u_h, 0);

                for(int j = 0; j < option.dim; j++)
                    global.subdomain.nodal_displacement_sd[point][i][j] = u_h[j];                
            }
        }
    }

}

void increment_field(){
    double stresses[6];
    //変位を更新し、変位増分をゼロ処理
    for(int point = 0;  point < global.subdomain.N_point; point++){
        for(int i = 0; i < option.dim; i++){
            global.subdomain.displacement[point][i] += global.subdomain.displacement_increment[point][i];
            global.subdomain.displacement_increment[point][i] = 0.;
        }
    }
    
    for(int point = 0; point < global.subdomain.N_point; point++){

        if(option.solver_type == 1){
           for(int i = 0; i < option.dim; i++){
                for(int j = 0; j < option.dim; j++){
                    global.subdomain.deformation_gradients[i][j][point] = global.subdomain.current_deformation_gradients[i][j][point];
                }
            }
        }

        for(int i = 0; i < 6; i++){
            global.subdomain.elastic_strains[point][i] = global.subdomain.current_elastic_strains[point][i];
            global.subdomain.stresses[point][i] = global.subdomain.current_stresses[point][i];
            global.subdomain.back_stresses[point][i] = global.subdomain.current_back_stresses[point][i];

            stresses[i] = global.subdomain.stresses[point][i];  //相当応力の計算のためにポイントごとの応力ベクトルを用意
        }

        global.subdomain.equivalent_stresses[point] = calc_equivalent_stress(stresses);

        update_plastic_strains(global.subdomain.current_plastic_strains[point], global.subdomain.stresses[point], global.subdomain.equivalent_stresses[point], global.subdomain.equivalent_plastic_strain_increments[point]);
        global.subdomain.equivalent_plastic_strains[point] += global.subdomain.equivalent_plastic_strain_increments[point];
        global.subdomain.equivalent_plastic_strain_increments[point] = 0.;

        global.subdomain.yield_stresses[point] = global.subdomain.current_yield_stresses[point];        
    }

    for(int node = 0; node < global.subdomain.N_node; node++){
        for(int i = 0; i < option.dim; i++){
            global.subdomain.nodal_displacements[node][i] += global.subdomain.nodal_displacement_increments[node][i];
            global.subdomain.nodal_displacement_increments[node][i] = 0.;
        }
    }
    
}

double incremental_deformation_norm(){
    double deformation_variation_norm = 0.;
    double deformation_norm = 0.;
    double **displacement_increment_variation;
    displacement_increment_variation = matrix(global.subdomain.N_point, option.dim);

    for(int point = 0; point < global.subdomain.N_point; point++){
        for(int i = 0; i < option.dim; i++){
            displacement_increment_variation[point][i] = global.subdomain.previous_displacement_increment[point][i] - global.subdomain.displacement_increment[point][i];
        }
    }

    deformation_variation_norm = norm_for_mat(displacement_increment_variation, global.subdomain.N_point, option.dim);
    deformation_norm = norm_for_mat(global.subdomain.displacement_increment, global.subdomain.N_point, option.dim);

    free_matrix(displacement_increment_variation);

    return  deformation_variation_norm / deformation_norm;    
}

void update_plastic_strains(double plastic_strains[6], const double stresses[6], const double equivarent_stresses, const double equivarent_strain_increment){
    double current_relative_hydrostatic_stress;
    double deviatoric_stresses[6];
    
    current_relative_hydrostatic_stress 
    = (1.0 / 3.0) * (stresses[0] + stresses[1] + stresses[2]);
    deviatoric_stresses[0] = stresses[0] - current_relative_hydrostatic_stress;
    deviatoric_stresses[1] = stresses[1] - current_relative_hydrostatic_stress;
    deviatoric_stresses[2] = stresses[2] - current_relative_hydrostatic_stress;
    for(int i = 3; i < 6; i++){
        deviatoric_stresses[i] = stresses[i];
    }

    double factor
            = equivarent_strain_increment * 1.5 / equivarent_stresses;
    for(int i = 0; i < 3; i++)
        plastic_strains[i] 
        += factor * deviatoric_stresses[i];
    for(int i = 3; i < 6; i++)
        plastic_strains[i] 
        += 2.0 * factor * deviatoric_stresses[i];
}

void cut_back(){
    for(int point = 0;  point < global.subdomain.N_point; point++){
        for(int i = 0; i < option.dim; i++){
            global.subdomain.displacement_increment[point][i] = 0.;
            global.subdomain.equivalent_plastic_strain_increments[point] = 0.;
        }
    }
    for(int node = 0; node < global.subdomain.N_node; node++){
        for(int i = 0; i < option.dim; i++){
            global.subdomain.nodal_displacement_increments[node][i] = 0.;
        }
    }    
}

void whether_points_is_in_the_subdomain(){
    double face_node_XYZ[4][3];
    int node_id[4];
    int subdomain_node[8];
    double edge_1[3]; 
    double edge_2[3];
    double cross[3];
    double point1[3];
    double point2[3];
    double current_point_XYZ1[3];
    double current_point_XYZ2[3];
    int node_id_ref[4];
    int node_id_nei[4];
    double subdomain_node_ref[8];
    double subdomain_node_nei[8];
    double face_node_XYZ_ref[4][3];
    double face_node_XYZ_nei[4][3];

    for(int i = 0; i < global.subdomain.N_int_boundary; i++){
        int reference_point = global.subdomain.pair_point_ib[2*i];
        int neighbor_point = global.subdomain.pair_point_ib[2* i + 1];

        for(int j = 0; j < option.dim; j++){
            current_point_XYZ1[j] = global.subdomain.point_XYZ[option.dim * reference_point + j] + global.subdomain.displacement[reference_point][j];
            current_point_XYZ2[j] = global.subdomain.point_XYZ[option.dim * neighbor_point + j] + global.subdomain.displacement[neighbor_point][j];
        }

        for(int j = 0; j < 8; j++){
            subdomain_node[j] = global.subdomain.subdomain_node[8 * reference_point + j];
        }

        generate_node_id(global.subdomain.shared_face[i], reference_point, subdomain_node, node_id);
        
        generate_current_face_node(face_node_XYZ, node_id, subdomain_node, reference_point);

        for(int j = 0; j < option.dim; j++){
            edge_1[j] = face_node_XYZ[3][j] - face_node_XYZ[1][j];
            edge_2[j] = face_node_XYZ[0][j] - face_node_XYZ[1][j];
            point1[j] = current_point_XYZ1[j] - face_node_XYZ[1][j];
            point2[j] = current_point_XYZ2[j] - face_node_XYZ[1][j];
        }
        cross_product(option.dim, edge_1, edge_2, cross);
        double a = dot_product(option.dim, cross, point1);
        double b = dot_product(option.dim, cross, point2);

        if((a > 0 && b > 0) || (a < 0 && b < 0)){
            printf("overray was occured! ===> reference => %5d neighbor => %5d\n", reference_point, neighbor_point);

        for(int j = 0; j < option.dim; j++){
            printf("%+15.14e    ", global.subdomain.point_XYZ[option.dim * reference_point + j]
                                    + global.subdomain.displacement[reference_point][j]);
        }
        printf("\n");
            // for(int j = 0; j < 8; j++){
            //     subdomain_node_ref[j] = global.subdomain.subdomain_node[8 * reference_point + j];
            //     subdomain_node_nei[j] = global.subdomain.subdomain_node[8 * neighbor_point + j];
            // }

            // generate_node_id(i, reference_point, subdomain_node_ref, node_id_ref);
            // generate_node_id(i, neighbor_point, subdomain_node_nei, node_id_nei);

            // generate_current_face_node(face_node_XYZ_ref, node_id_ref, subdomain_node_ref, reference_point);
            // generate_current_face_node(face_node_XYZ_nei, node_id_nei, subdomain_node_nei, neighbor_point);
            // for(int j = 0; j < option.dim; j++){
            //     printf("%+15.14e    ", global.subdomain.point_XYZ[option.dim * reference_point + j]);
            // }
            // printf("\n");
            
            // for(int j = 0; j < 4; j++){
            //     for(int k = 0; k < option.dim; k++){
            //         printf("%+15.14e    ", face_node_XYZ_ref[j][k]);
            //     }
            //     printf("\n");
            // }
            
            // for(int j = 0; j < option.dim; j++){
            //     printf("%+15.14e    ", global.subdomain.point_XYZ[option.dim * neighbor_point + j]);
            // }
            
            // for(int j = 0; j < 4; j++){
            //     for(int k = 0; k < option.dim; k++){
            //         printf("%+15.14e    ", face_node_XYZ_nei[j][k]);
            //     }
            //     printf("\n");
            // }

        }

        
    }
}