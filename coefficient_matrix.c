#pragma warning(disable: 4100) // 引数が未使用の場合
#pragma warning(disable: 4189) // ローカル変数が未使用の場合
#pragma warning(disable: 4996) //fopenの警告番号

#include<stdio.h>
#include<stdlib.h>
#include<omp.h>
#include"type.h"
#include"coefficient_matrix.h"
#include"scalar.h"
#include"vector.h"
#include"matrix.h"
#include"b_matrix.h"
#include"d_matrix.h"
#include"tensor.h"
#include"s_matrix.h"
#include"GetGaussPoints.h"
#include"external_force.h"

extern Global global;
extern Option option;

#define NUMBER_OF_NODE_IN_SUBDOMAIN 8
#define NUMBER_OF_NODE_IN_FACE 4

void generate_coefficient_matrix(){
    double current_deformation_gradient[3][3];                  //現配置での変形勾配テンソル
    long long DoF_Free = option.dim * global.subdomain.N_point; //マトリクスの自由度数
    long long size = DoF_Free * DoF_Free;                       //接線係数マトリクスの要素数

    //全体剛性マトリクスを初期化
    for(long long i = 0; i < DoF_Free; i++)
        for(long long j = 0; j < DoF_Free; j++){
            long long size_K = DoF_Free * i + j;
            global.subdomain.Global_K[size_K] = 0.;
        }
    printf("Coefficient1\n");
    //接線剛性マトリクスの領域積分の項を計算
    for(int point = 0; point < global.subdomain.N_point; point++){
        if(option.solver_type == 1){
            for(int i = 0; i < option.dim; i++)
                for(int j = 0; j < option.dim; j++)
                    current_deformation_gradient[i][j] = global.subdomain.current_deformation_gradients[i][j][point];
        }
    
        generate_subdomain_coefficient_matrix_for_volume(point, current_deformation_gradient, global.subdomain.current_stresses[point], global.subdomain.trial_elastic_strains[point],
        global.subdomain.equivalent_plastic_strains, global.subdomain.equivalent_plastic_strain_increments,global.subdomain.back_stresses[point]);
    }
    printf("Coefficient2\n");
    //ペナルティ項の第2, 3項を計算
    for(int face = 0; face < global.subdomain.N_int_boundary; face++){
        for(int i = 0; i < 2; i++){
            if(option.solver_type == 1){
                for(int k = 0; k < option.dim; k++)
                    for(int l = 0; l < option.dim;l++)
                        current_deformation_gradient[k][l] = global.subdomain.current_deformation_gradients[k][l][global.subdomain.pair_point_ib[2 * face + i]];
            }
            
            for(int j = 0; j < 2; j++)
                generate_subdomain_coefficient_matrix_for_PenaltyTerm(global.subdomain.pair_point_ib[2 * face + j],global.subdomain.pair_point_ib[2 * face + i], face,
                                                    current_deformation_gradient, global.subdomain.current_stresses[global.subdomain.pair_point_ib[2 * face + i]], global.subdomain.trial_elastic_strains[global.subdomain.pair_point_ib[2 * face + i]],
                                                    global.subdomain.equivalent_plastic_strains, global.subdomain.equivalent_plastic_strain_increments, global.subdomain.back_stresses[global.subdomain.pair_point_ib[2 * face + i]], (j+1)%2);
            
        }
    }
    printf("Coefficient3\n");
    //ペナルティ項（安定化項）の項を計算
    for(int face = 0; face < global.subdomain.N_int_boundary; face++)
        for(int i = 0; i < 2; i++)
            for(int j = 0; j < 2; j++)
                generate_subdomain_coefficient_matrix_for_StabilizationTerm(global.subdomain.pair_point_ib[2 * face + i], global.subdomain.pair_point_ib[2 * face + j], face, 2 * i + j);
  
}

void generate_subdomain_coefficient_matrix_for_volume(const int point_n,
                                            double (*current_deformation_gradients)[3], const double *current_stress, const double *trial_elastic_strains,
                                            const double *equivalemt_plastic_strains, const double *equivalent_plastic_strain_increments, const double *back_stresses){
    double b_t_matrix[180][6];                   //bマトリクス
    double b_t_NL_matrix[180][9];                //Bマトリクス（初期応力項用
    double *current_point_XYZ;
    double **G;
    double s_matrix[9][9];                      //応力マトリクス（初期応力項用）
    double d_matrix[6][6];                      //材料定数マトリクス
    double concictent_d_matrix[3][3][3][3];     //コンシステント接線剛性マトリクス
    double mapping_parameter;                            //体積 
    double ke_matrix[180][180];
    double BTD[180][6];
    double GTS[180][9];
    FILE *fp_debug;
    char FILE_name[128];

    int N_support = global.subdomain.support_offset[point_n + 1] - global.subdomain.support_offset[point_n];

    //形状を更新する場合と更新しない場合で条件分岐
    if(option.solver_type == 1){    
        mapping_parameter = calc_subdomain_volume(point_n);
    }else{
        mapping_parameter = calc_initial_subdomain_volume(point_n);
    }
    
    //ke_matrixをゼロ処理
    for(int i = 0; i < option.dim * (N_support + 1); i++)
        for(int j = 0; j < option.dim * (N_support + 1); j++)
            ke_matrix[i][j] = 0.;

    //bマトリクスの計算
    generate_linear_b_matrix(b_t_matrix, point_n);
    
    if(equivalent_plastic_strain_increments[point_n] > 0.){
        generate_elastic_plastic_d_matrix(d_matrix, trial_elastic_strains, equivalemt_plastic_strains[point_n], equivalent_plastic_strain_increments[point_n], back_stresses);
    }else{
        generateElasticDMatrix(d_matrix);
    }
    
    if(option.solver_type == 1){
        //有限ひずみのDマトリクスに修正
        modify_d_matrix_with_finite_strain(d_matrix, current_stress, trial_elastic_strains, current_deformation_gradients);
    }
    
    for(int i = 0; i < option.dim * (N_support + 1); i++){
        for(int j = 0; j < option.dim * (N_support + 1); j++){
            double ke_ij = 0.;
            for(int k = 0; k < 6; k++){
                double BTD_ik = 0.;
                for(int l = 0; l < 6; l++)
                    BTD_ik += b_t_matrix[i][l] * d_matrix[l][k];
                ke_ij += BTD_ik * b_t_matrix[j][k];
            }
            ke_matrix[i][j] += ke_ij * mapping_parameter;
        }
    }

    if(option.solver_type == 1){
        //非線形Bマトリクスの計算
        if((current_point_XYZ = (double *)calloc(option.dim * global.subdomain.N_point, sizeof(double))) == NULL){
            printf("Error:current_point_XYZ's memory is not enough\n");
            exit(-1);
        }
        for(int i = 0; i < global.subdomain.N_point; i++){
            for(int j = 0; j < option.dim; j++){
                current_point_XYZ[option.dim * i + j] = global.subdomain.point_XYZ[option.dim * i + j]
                                                    + global.subdomain.displacement[i][j]
                                                    + global.subdomain.displacement_increment[i][j];
            }
        }

        G = matrix(option.dim * option.dim, option.dim * (N_support + 1));
        calc_G(option.dim, point_n, current_point_XYZ, global.subdomain.support_offset, global.subdomain.support, G);

        //Sマトリクスの計算
        generateSMatrix(s_matrix, current_stress);
        
        #if 1
        for(int i = 0; i < option.dim * (N_support + 1); i++){
            for(int j = 0; j < option.dim * (N_support + 1); j++){
                double ke_ij = 0.;
                for(int k = 0; k < 9; k++){
                    double BTS_ik = 0.;
                    for(int l = 0; l < 9; l++)
                        BTS_ik += G[l][i] * s_matrix[l][k];
                    ke_ij += BTS_ik * G[k][j];
                }
                ke_matrix[i][j] += ke_ij * mapping_parameter;
            }
        }
        #endif

        free_matrix(G);
        free(current_point_XYZ);
    }
    
    assemble_coefficient_matrix(ke_matrix, global.subdomain.Global_K, point_n, point_n);
}

/*
引数内での変形勾配テンソル等はpoint_n2のサブドメイン内で定義したもの
*/
void generate_subdomain_coefficient_matrix_for_PenaltyTerm(const int point_n1, const int point_n2, const int face_n,
                                            double (*current_deformation_gradients)[3], const double *current_stress, double *trial_elastic_strains,
                                            const double *equivalemt_plastic_strains, const double *equivalent_plastic_strain_increments, const double *back_stresses, const int flag)
{
    const int N_qu = 1;
    double sign;                                //項の符号
    int end_point;                              //法線ベクトルの方向を決定するためのポイント番号
    double xyz[3];                              //ガウスポイントの座標
    double X[27], w[27];                        //正規化座標、重み関数
    double mapping_parameter;                   //正規化座標→物理座標への変換パラメーター
    int face_node[NUMBER_OF_NODE_IN_FACE];      //面内頂点番号
    double face_node_XYZ[4][3];                 //↑の座標
    double b_t_matrix[180][6];                   //ｂマトリクス
    double b_t_NL_matrix[180][9];                //非線形ｂマトリクス
    double NT[180][3];                           //形状関数                        
    double Ne[3][6];                            //法線ベクトルの成分を格納したマトリクス
    double Ne_d[3][9];                          //Neを異なる形式で格納したマトリクス
    double d_matrix[6][6];                      //ｄマトリクス
    double s_matrix[9][9];                      //Sマトリクス
    double *current_point_XYZ;                  //ポイントの現配置座標
    double **G;
    double ke_matrix[180][180];                   //要素剛性マトリクス

    int N1_support = global.subdomain.support_offset[point_n1 + 1] - global.subdomain.support_offset[point_n1];
    int N2_support = global.subdomain.support_offset[point_n2 + 1] - global.subdomain.support_offset[point_n2];

    //ke_matrixをゼロ処理
    for(int i = 0; i < option.dim * (N1_support + 1); i++)
        for(int j = 0; j < option.dim * (N2_support + 1); j++)
            ke_matrix[i][j] = 0.;
        
    //項の符号を判定
    if(point_n1 == point_n2){
        sign = -1.0;
    }else{
        sign = 1.0;
    }
    //ポイントの座標を計算
    if((current_point_XYZ = (double *)calloc(option.dim * global.subdomain.N_point, sizeof(double))) == NULL){
        printf("Error:current_point_XYZ's memory is not enough\n");
        exit(-1);
    }
    if(option.solver_type == 1){
       for(int i = 0; i < global.subdomain.N_point; i++){
            for(int j = 0; j < option.dim; j++){
                current_point_XYZ[option.dim * i + j] = global.subdomain.point_XYZ[option.dim * i + j]
                                                    + global.subdomain.displacement[i][j]
                                                    + global.subdomain.displacement_increment[i][j];
            }
        } 
    }else{
        for(int i = 0; i < global.subdomain.N_point; i++)
            for(int j = 0; j < option.dim; j++)
                current_point_XYZ[option.dim * i + j] = global.subdomain.point_XYZ[option.dim * i + j];
    }

    //ガウスポイントを計算
    Gauss_points_and_weighting_factors(N_qu, X, w);

    //法線ベクトルを計算
    if(point_n2 == global.subdomain.pair_point_ib[2 * face_n]){
        end_point = global.subdomain.pair_point_ib[2*face_n+1];
    }else if(point_n2 == global.subdomain.pair_point_ib[2*face_n+1]){
        end_point = global.subdomain.pair_point_ib[2*face_n];
    }
    generate_unit_vec_to_mat3x6(global.subdomain.shared_face[face_n], point_n2, end_point, current_point_XYZ, Ne);
    
    //Bマトリクスの計算
    generate_linear_b_matrix(b_t_matrix, point_n2);

    //材料定数マトリクスの計算
    if(equivalent_plastic_strain_increments[point_n2] > 0.){
        generate_elastic_plastic_d_matrix(d_matrix, trial_elastic_strains, equivalemt_plastic_strains[point_n2], equivalent_plastic_strain_increments[point_n2], back_stresses);
    }else{
        generateElasticDMatrix(d_matrix);
    }
   
    if(option.solver_type == 1){
        //接線係数をコンシステント接線係数に修正
        modify_d_matrix_with_finite_strain(d_matrix, current_stress, trial_elastic_strains, current_deformation_gradients);
     
        //節点の現配置座標を計算
        generate_current_node_of_face(face_node_XYZ, global.subdomain.shared_face[face_n], point_n2);

    }else{
        //面内に含まれる頂点の座標を計算
        for(int i = 0; i < NUMBER_OF_NODE_IN_FACE; i++){
            int v_offset = global.subdomain.vertex_offset[global.subdomain.shared_face[face_n]];
            face_node[i] = global.subdomain.node[v_offset + i];
        }

        for(int i = 0; i < NUMBER_OF_NODE_IN_FACE; i++)
            for(int j = 0; j < option.dim; j++)
                face_node_XYZ[i][j] = global.subdomain.node_XYZ[option.dim * face_node[i] + j];
    }
    
    //積分
    for(int s = 0; s < N_qu; s++){
        for(int t = 0; t < N_qu; t++){
            //正規化座標→物理座標のマッピングパラメーターを計算
            mapping_parameter = calc_mapping_parameter(global.subdomain.shared_face[face_n], point_n2, s, t, X);
            
            //ガウス点のz座標を計算
            generate_gauss_point_coordinate(s, t, face_node_XYZ, X, xyz);

            //形状関数の計算
            calc_shape(xyz, option.dim, point_n1, current_point_XYZ, global.subdomain.support_offset, NT);
            
            //要素剛性マトリクスの計算
            for(int i = 0; i < option.dim * (N1_support + 1); i++){
                for(int j = 0; j < option.dim * (N2_support + 1); j++){
                    double Ke_ij = 0.;
                    for(int k = 0; k < 6; k++){
                        double NT_Ne_D_ik = 0.;
                        
                        for(int l = 0; l < 6; l++){
                            double NT_Ne_il = 0.;

                            for(int m = 0; m < option.dim; m++){
                                NT_Ne_il += NT[i][m] * Ne[m][l];
                            }
                            NT_Ne_D_ik += NT_Ne_il * d_matrix[l][k];
                        }
                        Ke_ij += NT_Ne_D_ik * b_t_matrix[j][k];
                    }
                    ke_matrix[i][j] += 0.5 * Ke_ij * mapping_parameter * sign * w[s] * w[t];
                }
            }
        }
    }

    if(option.solver_type == 1){
        //法線ベクトルを計算
        generate_unit_vec_to_mat3x9(global.subdomain.shared_face[face_n], point_n2, end_point, current_point_XYZ, Ne_d);

        //Bマトリクスの計算
        //generate_nonlinear_b_matrix(b_t_NL_matrix, point_n2);
        G = matrix(option.dim * option.dim, option.dim * (N2_support + 1));
        calc_G(option.dim, point_n2, current_point_XYZ, global.subdomain.support_offset, global.subdomain.support, G);

        //Sマトリクスの計算
        generateSMatrix(s_matrix, current_stress);
        
        for(int s = 0; s < N_qu; s++){
            for(int t = 0; t < N_qu; t++){
                //正規化座標→物理座標のマッピングパラメーターを計算
                mapping_parameter = calc_mapping_parameter(global.subdomain.shared_face[face_n], point_n2, s, t, X);

                //ガウス点の座標を計算
                generate_gauss_point_coordinate(s, t, face_node_XYZ, X, xyz);

                //形状関数の計算
                calc_shape(xyz, option.dim, point_n1, current_point_XYZ, global.subdomain.support_offset, NT);

                //要素剛性マトリクスの計算
                for(int i = 0; i < option.dim * (N1_support + 1); i++){
                    for(int j = 0; j < option.dim * (N2_support + 1); j++){
                        double Ke_ij = 0.;

                        for(int k = 0; k < 9; k++){
                            double Nt_Ne_d_S_ik = 0.;

                            for(int l = 0; l < 9; l++){
                                double Nt_Ne_d_il = 0.;

                                for(int m = 0; m < option.dim; m++){
                                    Nt_Ne_d_il += NT[i][m] * Ne_d[m][l];
                                }
                                Nt_Ne_d_S_ik += Nt_Ne_d_il * s_matrix[l][k]; 
                            }
                            Ke_ij += Nt_Ne_d_S_ik * G[k][j];
                        }
                        ke_matrix[i][j] += 0.5 * Ke_ij * mapping_parameter * sign * w[s] * w[t];
                    }
                }
            }
        }
        free_matrix(G);
    }

    //全体剛性マトリクスにアセンブル
    assemble_coefficient_matrix(ke_matrix, global.subdomain.Global_K, point_n1, point_n2);

    free(current_point_XYZ);
}

void generate_subdomain_coefficient_matrix_for_StabilizationTerm(const int point_n1, const int point_n2, const int face_n, const int flag){
int N_qu = 2;                                                               //ガウス点数2x2tンを使用
    long long DoF_Free = option.dim * global.subdomain.N_point;
    double N1T[180][3];
    double N2T[180][3];
    double xyz[3];
    double X[27], w[27];
    double sign;
    double *current_point_XYZ;
    double jump_u[3];
    int end_point;
    double Ns[9];
    double A[9][9];
    double Ls[9][180];
    double **G;
    int face_node[NUMBER_OF_NODE_IN_FACE];
    double face_node_XYZ[4][3];
    double face_node_XYZ1[4][3];
    double face_node_XYZ2[4][3];
    double he;                                                  //ポイント間の距離
    double eta = global.material.penalty;                       //ペナルティパラメータ
    double mapping_parameter;
    double ke_matrix[180][180];
    int N1_support = global.subdomain.support_offset[point_n1 + 1] - global.subdomain.support_offset[point_n1];
    int N2_support = global.subdomain.support_offset[point_n2 + 1] - global.subdomain.support_offset[point_n2];

    //ke_matrixをゼロ処理
    for(int i = 0; i < option.dim * (N1_support + 1); i++)
        for(int j = 0; j < option.dim * (N2_support + 1); j++)
            ke_matrix[i][j] = 0.;

    //項の符号を与える
    if(flag == 0 || flag == 3){
        sign = 1.0;
    }
    else if(flag == 1 || flag == 2){
        sign = -1.0;
    }


    //ポイントの座標を計算
    if((current_point_XYZ = (double *)calloc(option.dim * global.subdomain.N_point, sizeof(double))) == NULL){
        printf("Error: current_point_XYZ memory is not enough\n");
        exit(-1);
    }
    
    if(option.solver_type == 1){
       for(int i = 0; i < global.subdomain.N_point; i++)
            for(int j = 0; j < option.dim; j++)
                current_point_XYZ[option.dim * i + j] = global.subdomain.point_XYZ[option.dim * i + j]
                                                    + global.subdomain.displacement[i][j]
                                                    + global.subdomain.displacement_increment[i][j];
    }else{
        for(int i = 0; i < global.subdomain.N_point; i++)
            for(int j = 0; j < option.dim; j++)
                current_point_XYZ[option.dim * i + j] = global.subdomain.point_XYZ[option.dim * i + j];
    }

    //ガウス積分点と重み係数の設定
    Gauss_points_and_weighting_factors(N_qu, X, w);

    //Γ*の頂点の座標を計算（Γ+とΓ-の平均を計算
    if(option.solver_type == 1){
        generate_current_node_of_face(face_node_XYZ1, global.subdomain.shared_face[face_n], global.subdomain.pair_point_ib[2 * face_n]);
        generate_current_node_of_face(face_node_XYZ2, global.subdomain.shared_face[face_n], global.subdomain.pair_point_ib[2 * face_n + 1]);
        for(int i = 0; i < NUMBER_OF_NODE_IN_FACE; i++){
            for(int j = 0; j < option.dim; j++){
                face_node_XYZ[i][j] = 0.5 * (face_node_XYZ1[i][j] + face_node_XYZ2[i][j]);
            }
        }
    }else{
        //面内に含まれる頂点の座標を計算
        for(int i = 0; i < NUMBER_OF_NODE_IN_FACE; i++){
            int v_offset = global.subdomain.vertex_offset[global.subdomain.shared_face[face_n]];
            face_node[i] = global.subdomain.node[v_offset + i];
        }

        for(int i = 0; i < NUMBER_OF_NODE_IN_FACE; i++)
            for(int j = 0; j < option.dim; j++)
                face_node_XYZ[i][j] = global.subdomain.node_XYZ[option.dim * face_node[i] + j];
    }
    
  
    he = distance(option.dim, global.subdomain.pair_point_ib[2 * face_n], global.subdomain.pair_point_ib[2 * face_n + 1], current_point_XYZ);
 
    for(int s = 0; s < N_qu; s++){
        for(int t = 0; t < N_qu; t++){
            //物理空間座標→正規化座標に変換するためのスカラー値を計算
            mapping_parameter = calc_mapping_parameter_for_av_area(face_node_XYZ, s, t, X);
            
            //物理座標におけるガウス点の座標を計算
            generate_gauss_point_coordinate(s, t, face_node_XYZ, X, xyz);
        
            //サブドメイン番号point_n1とpoint_n1の形状関数を計算
            calc_shape(xyz, option.dim, point_n1, current_point_XYZ, global.subdomain.support_offset, N1T);
           
            calc_shape(xyz, option.dim, point_n2, current_point_XYZ, global.subdomain.support_offset, N2T);

            //要素剛性マトリクスの計算
            for(int i = 0;  i < option.dim * (N1_support + 1); i++){
                for(int j = 0; j < option.dim * (N2_support + 1); j++){
                    double N1TN2_ij = 0.;
                    for(int k = 0; k < option.dim; k++){
                        N1TN2_ij += N1T[i][k] * N2T[j][k];
                    }
                    ke_matrix[i][j] += eta / he * N1TN2_ij * mapping_parameter * sign  * w[s] * w[t];
                }
            }
        }
    }
  
    if(option.solver_type == 1){
        if (flag > 1){
            sign = -1.0;
        }else if(flag < 2){
            sign = 1.0;
        }
        //法線ベクトルを計算
        generate_unit_vec_to_mat1x9(global.subdomain.shared_face[face_n], face_node_XYZ, 
        global.subdomain.pair_point_ib[2*face_n], global.subdomain.pair_point_ib[2*face_n+1], current_point_XYZ, Ns);

        //Lsを計算
        cross_minus_1x9(A);
        G = matrix(option.dim * option.dim, option.dim * (N2_support + 1));
        calc_G(option.dim, point_n2, current_point_XYZ, global.subdomain.support_offset, global.subdomain.support, G);
        for(int i = 0; i < 9; i++){
            for(int j = 0; j < option.dim * (N2_support + 1); j++){
                double LS_ij = 0.;
                for(int k = 0; k < 9; k++){
                    LS_ij += A[i][k] * G[k][j];
                }
                Ls[i][j] = LS_ij;
            }
        }

        for(int s = 0; s < N_qu; s++){
            for(int t = 0; t < N_qu; t++){
                //物理空間座標→正規化座標に変換するためのスカラー値を計算
                mapping_parameter = calc_mapping_parameter_for_av_area(face_node_XYZ, s, t, X);

                //物理座標におけるガウス点の座標を計算
                generate_gauss_point_coordinate(s, t, face_node_XYZ, X, xyz);
            
                //サブドメイン番号point_n1とpoint_n1の形状関数を計算
                calc_shape(xyz, option.dim, point_n1, current_point_XYZ, global.subdomain.support_offset, N1T);
           
                jump_trial_u(xyz, global.subdomain.pair_point_ib[2*face_n], global.subdomain.pair_point_ib[2*face_n+1], current_point_XYZ, jump_u, 0);

                for(int i = 0; i < option.dim * (N1_support + 1); i++){
                    for(int j = 0; j < option.dim * (N2_support + 1); j++){
                        double ke_ij = 0.;
                        for(int k = 0; k < 9; k++){
                            double N_u_Ns_ik = 0.;
                            for(int l = 0; l < option.dim; l++){
                                N_u_Ns_ik += N1T[i][l] * jump_u[l] * Ns[k];
                            }
                            ke_ij += N_u_Ns_ik * Ls[k][j];
                        }
                        ke_matrix[i][j] += eta / he * 0.5 * ke_ij * mapping_parameter * sign * w[s] * w[t];
                    }
                }
            }
        }
        free_matrix(G);

    }

    //全体剛性マトリクスにアセンブリ
    assemble_coefficient_matrix(ke_matrix, global.subdomain.Global_K, point_n1, point_n2);

    free(current_point_XYZ);

}

void assemble_coefficient_matrix(double (*element_K)[180], double *Global_K, int point_n1, int point_n2){
    int ref_num1 = global.subdomain.support_offset[point_n1];
    int ref_num2 = global.subdomain.support_offset[point_n2];
    int N1_support = global.subdomain.support_offset[point_n1 + 1] - global.subdomain.support_offset[point_n1];
    int N2_support = global.subdomain.support_offset[point_n2 + 1] - global.subdomain.support_offset[point_n2];
    long long DoF_free = option.dim * global.subdomain.N_point;
    
    for(int i = 0; i < option.dim; i++){
        for(int j = 0; j < option.dim; j++){
            long long size_K =  DoF_free * (option.dim * point_n1 + i) + option.dim * point_n2 + j;
            Global_K[size_K] += element_K[i][j];
        }
    }
    for(int i = 0; i < N2_support; i++){
        for(int j = 0; j < option.dim; j++){
            for(int k = 0; k < option.dim; k++){
                long long size_K = DoF_free * (option.dim * point_n1 + j) + option.dim * global.subdomain.support[ref_num2 + i] + k;
                Global_K[size_K] += element_K[j][option.dim * (i + 1) + k];
            }
        }
    }
    for(int i = 0; i < N1_support; i++){
        for(int j = 0; j < option.dim; j++){
            for(int k = 0; k < option.dim; k++){
                long long size_K = DoF_free * (option.dim * global.subdomain.support[ref_num1 + i] + j) + option.dim * point_n2 + k;
                Global_K[size_K] += element_K[option.dim * (i + 1) + j][k];
            }
        }
    }
    for(int i = 0; i < N1_support; i++){
        for(int j = 0; j < N2_support; j++){
            for(int k = 0; k < option.dim; k++){
                for(int l = 0; l < option.dim; l++){
                    long long size_K = DoF_free * (option.dim * global.subdomain.support[ref_num1 + i] + k) + (option.dim * global.subdomain.support[ref_num2 + j] + l);
                    Global_K[size_K] += element_K[option.dim * (i + 1) + k][option.dim * (j + 1) + l];
                }
            }
        }
    }
    
}
