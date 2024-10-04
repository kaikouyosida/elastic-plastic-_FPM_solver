#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"scalar.h"
#include"tensor.h"
#include"type.h"
#include"vector.h"
#include"matrix.h"

#define NUMBER_OF_NODE_IN_FACE 4
#define NUMBER_OF_NODE_IN_SUBDOMAIN 8

extern Global global;
extern Option option;

double calculateInverse(double value){
    return 1.0 / value;
}

///二乗の計算
double calculateSquare(double value){
    return value * value;
}

//3乗の計算
double calculateCube(double value){
    return value * value * value;
}

//double値の入れ替え
void swapReals(double *value1, double *value2){
    double temp = *value1;
    *value1  = *value2;
    *value2  = temp;
}

//int値の計算
void swapIntegers(int *value1, int *value2)
{
    int temp = *value1;
    *value1  = *value2;
    *value2  = temp;
}

//3×3行列のトレースをとる
double calculate3x3MatrixTrace(double matrix[3][3])
{
    return matrix[0][0] + matrix[1][1] + matrix[2][2];
}

//サブドメインの体積を計算する（六面体用）
double calc_subdomain_volume(int point_n){
    double volume = 0.;                     //サブドメインの体積
    double center[3];                       //サブドメインの重心
    double edge1[3],edge2[3],edge3[3];                            //辺のベクトル
    double edge1_cross_edge2[3];                                  //辺のベクトルの外積
    int subdomain_node[NUMBER_OF_NODE_IN_SUBDOMAIN];              //サブドメインpoint_nがもつ頂点番号
    int node_id[NUMBER_OF_NODE_IN_FACE];                          //節点番号のアドレス(subdomain_node用)

    int N_face = global.subdomain.face_offset[point_n + 1] - global.subdomain.face_offset[point_n];     //サブドメインpoint_nにおける面の数

    //6面体（節点8個）の重心を算術平均で計算
    generate_subdomain_node(point_n, subdomain_node);

    for(int i = 0; i < option.dim; i++){
        double center_i = 0.;
        for(int j = 0; j < NUMBER_OF_NODE_IN_SUBDOMAIN; j++){
            center_i += (global.subdomain.node_XYZ[option.dim * subdomain_node[j] + i]
                     + global.subdomain.nodal_displacement_sd[point_n][j][i]);
        }
        center[i] = center_i / 8.0;
    }

    //各面の頂点と重心の5点で四角すいを形成。その後各面で総和をとって六面体の面積を計算する
    for(int i = 0; i < N_face; i++){

        //内部境界面のノード番号のアドレスを格納
        generate_node_id(global.subdomain.face[global.subdomain.face_offset[point_n] + i], point_n, subdomain_node, node_id);
        
        generate_current_edge_vector(edge1, point_n, node_id[1], node_id[0], subdomain_node);
        generate_current_edge_vector(edge2, point_n, node_id[3], node_id[0], subdomain_node);
        generate_current_points_vector(edge3, center, point_n, node_id[0], subdomain_node);
        

        cross_product(option.dim, edge1, edge2, edge1_cross_edge2);
        volume += fabs(dot_product(option.dim,edge1_cross_edge2, edge3));
        
        generate_current_edge_vector(edge1, point_n, node_id[3], node_id[2], subdomain_node);
        generate_current_edge_vector(edge2, point_n, node_id[1], node_id[2], subdomain_node);
        generate_current_points_vector(edge3, center, point_n, node_id[2], subdomain_node);
        

        cross_product(option.dim, edge1, edge2, edge1_cross_edge2);
        volume += fabs(dot_product(option.dim,edge1_cross_edge2, edge3));
    }
    
    volume /= 6.0;
    
    return volume;
}

double calc_initial_subdomain_volume(int point_n){
    double volume = 0;
    int subdomain_node[NUMBER_OF_NODE_IN_SUBDOMAIN];
    double center[3];
    int face_node[NUMBER_OF_NODE_IN_FACE];
    double edge1[3], edge2[3];
    double node_to_center_edge[3];
    double edge1crossedge2[3];

    int N_face = global.subdomain.face_offset[point_n + 1] - global.subdomain.face_offset[point_n];
    
    //サブドメインに含まれる頂点番号を格納
    generate_subdomain_node(point_n, subdomain_node);

    //サブドメイン重心を計算
    for(int i = 0; i < option.dim; i++){
        double center_i = 0.;
        for(int j = 0; j < NUMBER_OF_NODE_IN_SUBDOMAIN; j++){
            center_i += global.subdomain.node_XYZ[option.dim * subdomain_node[j] + i];
        }
        center[i] = center_i / 8.0;
    }
    
    //辺ベクトルの計算
    for(int i = 0; i < N_face; i++){
        for(int j = 0; j < NUMBER_OF_NODE_IN_FACE; j++)
            face_node[j] 
                = global.subdomain.node[global.subdomain.vertex_offset[global.subdomain.face[global.subdomain.face_offset[point_n] + i]] + j];
        
        for(int j = 0; j < option.dim; j++){
            edge1[j] = global.subdomain.node_XYZ[option.dim * face_node[0] + j] - global.subdomain.node_XYZ[option.dim * face_node[1] + j];
            edge2[j] = global.subdomain.node_XYZ[option.dim * face_node[2] + j] - global.subdomain.node_XYZ[option.dim * face_node[1] + j];
            node_to_center_edge[j] = center[j] - global.subdomain.node_XYZ[option.dim * face_node[1] + j];
        }
        cross_product(option.dim, edge1, edge2, edge1crossedge2);
        volume += fabs(dot_product(option.dim, edge1crossedge2, node_to_center_edge));

        for(int j = 0; j < option.dim; j++){
            edge1[j] = global.subdomain.node_XYZ[option.dim * face_node[2] + j] - global.subdomain.node_XYZ[option.dim * face_node[3] + j];
            edge2[j] = global.subdomain.node_XYZ[option.dim * face_node[0] + j] - global.subdomain.node_XYZ[option.dim * face_node[3] + j];
            node_to_center_edge[j] = center[j] - global.subdomain.node_XYZ[option.dim * face_node[3] + j];
        }
        cross_product(option.dim, edge1, edge2, edge1crossedge2);
        volume += fabs(dot_product(option.dim, edge1crossedge2, node_to_center_edge));
    }
    volume /= 6.0;
    
    return volume;
}

//内部境界Γ*を物理座標→正規化座標にマッピングするためのパラメータを計算
double calc_mapping_parameter_for_av_area(const double face_node_XYZ[4][3], int s, int t, double *X){
    double result = 0.;
    double area_vector[3];
    double dx_ds[3];
    double dx_dt[3];

    //x_h = 1/4*(1-ξ*ξ_i)*(1-η*η_i)*x_iのξ, ηに対する微係数を計算.
    for(int i = 0; i < option.dim ; i++){
        dx_ds[i] = -0.25 * (1.0 - X[t]) * face_node_XYZ[0][i] - 0.25 * (1.0 + X[t]) * face_node_XYZ[1][i] + 0.25 * (1.0 + X[t]) * face_node_XYZ[2][i] + 0.25 * (1.0 - X[t]) * face_node_XYZ[3][i];
        dx_dt[i] = -0.25 * (1.0 - X[s]) * face_node_XYZ[0][i] + 0.25 * (1.0 - X[s]) * face_node_XYZ[1][i] + 0.25 * (1.0 + X[s]) * face_node_XYZ[2][i] - 0.25 * (1.0 + X[s]) * face_node_XYZ[3][i];
    }

    //面積変化率を計算 da/dA=norm((dx_h/dξ)×(dx_h/dη))
    cross_product(option.dim, dx_ds, dx_dt, area_vector);
    result = norm(area_vector, option.dim);
    return result;
}

//内部境界Γ+とΓ-を物理座標→正規化座標にマッピングするためのパラメータを計算
double calc_mapping_parameter(int face_n, int point_n, int s, int t, double *X){
    double area_vector[3];                                    //面積ベクトルの計算
    double dx_ds[3];                                          //微係数dx_h/dξ
    double dx_dt[3];                                          //微係数dx_h/dη
    double x1[3], x2[3], x3[3], x4[3];                        //面内の頂点座標
    int node_id[NUMBER_OF_NODE_IN_FACE];                      //頂点番号のアドレス
    int subdomain_node[NUMBER_OF_NODE_IN_SUBDOMAIN];          //サブドメイン中の節点番号

    if(option.solver_type == 1){
        //ノード番号のアドレスを取得
        generate_subdomain_node(point_n, subdomain_node);
        generate_node_id(face_n, point_n, subdomain_node, node_id);
        
        for(int i = 0; i < option.dim; i++){
            x1[i] = global.subdomain.node_XYZ[option.dim * subdomain_node[node_id[0]] + i] + global.subdomain.nodal_displacement_sd[point_n][node_id[0]][i];// + global.subdomain.nodal_displacement_increment_sd[point_n][node_id[0]][i];
            x2[i] = global.subdomain.node_XYZ[option.dim * subdomain_node[node_id[1]] + i] + global.subdomain.nodal_displacement_sd[point_n][node_id[1]][i];// + global.subdomain.nodal_displacement_increment_sd[point_n][node_id[1]][i];
            x3[i] = global.subdomain.node_XYZ[option.dim * subdomain_node[node_id[2]] + i] + global.subdomain.nodal_displacement_sd[point_n][node_id[2]][i];// + global.subdomain.nodal_displacement_increment_sd[point_n][node_id[2]][i];
            x4[i] = global.subdomain.node_XYZ[option.dim * subdomain_node[node_id[3]] + i] + global.subdomain.nodal_displacement_sd[point_n][node_id[3]][i];// + global.subdomain.nodal_displacement_increment_sd[point_n][node_id[3]][i];
        }
    }else{
        for(int i = 0; i < option.dim; i++){
            x1[i] = global.subdomain.node_XYZ[option.dim * global.subdomain.node[global.subdomain.vertex_offset[face_n]] + i];
            x2[i] = global.subdomain.node_XYZ[option.dim * global.subdomain.node[global.subdomain.vertex_offset[face_n] + 1] + i];
            x3[i] = global.subdomain.node_XYZ[option.dim * global.subdomain.node[global.subdomain.vertex_offset[face_n] + 2] + i];
            x4[i] = global.subdomain.node_XYZ[option.dim * global.subdomain.node[global.subdomain.vertex_offset[face_n] + 3] + i];
        }
    }
    
    //x_h = 1/4*(1-ξ*ξ_i)*(1-η*η_i)*x_iのξ, ηに対する微係数を計算.
    for(int i = 0; i < option.dim ; i++){
        dx_ds[i] = -0.25 * (1.0 - X[t]) * x1[i] - 0.25 * (1.0 + X[t]) * x2[i] + 0.25 * (1.0 + X[t]) * x3[i] + 0.25 * (1.0 - X[t]) * x4[i];
        dx_dt[i] = -0.25 * (1.0 - X[s]) * x1[i] + 0.25 * (1.0 - X[s]) * x2[i] + 0.25 * (1.0 + X[s]) * x3[i] - 0.25 * (1.0 + X[s]) * x4[i];
    }
     
    //面積変化率を計算 da/dA=norm((dx_h/dξ)×(dx_h/dη))
    cross_product(option.dim, dx_ds, dx_dt, area_vector);
    
    double result = norm(area_vector, option.dim);
    return result;
}

//面積変化率（dΓ*/dΓ0)の計算
double generate_area_change_parameter(int subdomain_n1, int subdomain_n2, int face_n, int *vertex_offset, double face_node_XYZ[4][3], double *center_xyz){
    double trial_elastic_left_cauchy_green_deformations1[3][3];     //trial strainから計算したサブドメイン１の左Cauchyグリーンテンソル
    double trial_elastic_strain_tensor1[3][3];                      //サブドメイン1の試行ひずみテンソル
    double trial_elastic_left_cauchy_green_deformations2[3][3];     //trial strainから計算したサブドメイン2の左Cauchyグリーンテンソル
    double trial_elastic_strain_tensor2[3][3];                      //サブドメイン2の試行ひずみテンソル
    double trial_elastic_strains1[6];                               //サブドメイン1のひずみテンソルのvoigt表記
    double trial_elastic_strains2[6];                               //サブドメイン2のひずみテンソルのvoigt表記
    double trial_elastic_left_cauchy_green_deformations[3][3];      //平均化された左Cauchyグリーンテンソル
    double deformation_gradients[3][3];                             //平均化された変形勾配テンソル
    double Ne[3];                                                   //法線ベクトル
    double inverse_deformation_gradient;                            //体積変化率
    double area_change_parameter = 0.;                              //Γ*の面積変化率

    for(int i = 0; i < 6; i++){
        trial_elastic_strains1[i] = global.subdomain.trial_elastic_strains[subdomain_n1][i];
        trial_elastic_strains2[i] = global.subdomain.trial_elastic_strains[subdomain_n2][i];
    }
    
    //subdomain1の試行弾性左コーシーグリーンテンソルの計算
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

    //subdomain2の試行弾性左コーシーグリーンテンソルの計算
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
    
    calc_unit_vector(Ne, global.subdomain.shared_face[face_n], subdomain_n1, subdomain_n2, center_xyz);

    //内部境界Γ*での左Cauchyグリーンテンソルを計算(B+とB-の平均値)
    for(int i = 0; i < option.dim; i++){
        for(int j = 0; j < option.dim; j++){
            trial_elastic_left_cauchy_green_deformations[i][j]
                = 0.5 * (trial_elastic_left_cauchy_green_deformations1[i][j] + trial_elastic_left_cauchy_green_deformations2[i][j]);
        }
    }

    //体積変化率の計算.内部境界を共有するサブドメインどうしの変形勾配の平均がとられ、そのデターミナントを計算.
    for(int i = 0; i < option.dim; i++){
        for(int j = 0; j < option.dim; j++){
            deformation_gradients[i][j] 
                = 0.5 * (global.subdomain.current_deformation_gradients[i][j][global.subdomain.pair_point_ib[2 * face_n]] 
                        + global.subdomain.current_deformation_gradients[i][j][global.subdomain.pair_point_ib[2 * face_n + 1]]);
        }
    }

    inverse_deformation_gradient = 1.0 / calc_3x3matrix_determinant(deformation_gradients);
    
    for(int i = 0; i < option.dim; i++){
        double nB_i = 0.;
        for(int j = 0; j < option.dim; j++){
            nB_i += Ne[j] * trial_elastic_left_cauchy_green_deformations[j][i];
        }
        area_change_parameter += nB_i * Ne[i];
    }

    area_change_parameter = inverse_deformation_gradient * sqrt(area_change_parameter);

    return area_change_parameter;
}

//内積を計算する
double dot_product(int N, double *vec1, double *vec2)
{
	double X = 0.;
	for (int i = 0; i < N; i++)
		X += vec1[i] * vec2[i];
	return X;
}

//外積を計算する
void cross_product(int dim, double *vecA, double *vecB, double *AcrossB)
{
	if (dim == 2)
		AcrossB[0] = vecA[0] * vecB[1] - vecA[1] * vecB[0];
	else if (dim == 3)
	{
		AcrossB[0] = vecA[1] * vecB[2] - vecA[2] * vecB[1];
		AcrossB[1] = vecA[2] * vecB[0] - vecA[0] * vecB[2];
		AcrossB[2] = vecA[0] * vecB[1] - vecA[1] * vecB[0];
	}
}

//ノルムを計算する
double norm(double *vec, int n)
{
	double s = 0.;
	s = dot_product(n, vec, vec);
	return sqrt(s);
}

//ノルムを計算する（二次配列用）
double norm_for_mat(double **vec, int m ,int n){
    double s = 0;
    for(int i = 0; i < m; i++){
        for(int j = 0; j < n; j++){
            s += vec[i][j] * vec[i][j];
        }
    }
    return sqrt(s);
}

//ポイントの距離を計算する
double distance(int dim, int i, int j, double *point_xyz)
{
	double d = 0.;
	for (int k = 0; k < dim; k++)
		d += point_xyz[dim * i + k] * point_xyz[dim * i + k] + point_xyz[dim * j + k] * point_xyz[dim * j + k];
	for (int k = 0; k < dim; k++)
		d -= 2.0 * point_xyz[dim * i + k] * point_xyz[dim * j + k];
	return sqrt(d);
}

//ひずみエネルギ密度のerrorを計算
double generate_strain_energy_rate_parameter(int time_step){
    double strain_sol[6];
    double stress_sol[6];
    double volume;
    double volume_global = 1.0;
    double strain_enrgy_rate = 0.;
    double strain_enrgy_rate_sol = 0.;
    double lamda;

    lamda = (double)(time_step + 1) / (double)option.N_timestep;

    //参照解の応力とひずみを計算
    for(int i = 0; i < 6; i++){
        strain_sol[i] = 0.; stress_sol[i] = 0.;
    }
    strain_sol[2] = 0.5 * log(9.0 / 4.0) * lamda; 
    stress_sol[0] = 1.55948118503140 * 100.0 * lamda;
    stress_sol[1] = 1.55948118503140 * 100.0 * lamda;
    stress_sol[2] = 3.63878943173994 * 100.0 * lamda;

    for(int point = 0; point < global.subdomain.N_point; point++){
        volume = calc_subdomain_volume(point);

        for(int i = 0; i < 6; i++){
            strain_enrgy_rate += (strain_sol[i] - global.subdomain.elastic_strains[point][i]) *(stress_sol[i] - global.subdomain.stresses[point][i]) * volume;
        }
    }

    for(int i = 0; i < 6; i++){
        strain_enrgy_rate_sol += strain_sol[i] * stress_sol[i];
    }
    strain_enrgy_rate_sol *= volume_global;

    return sqrt(strain_enrgy_rate) / sqrt(strain_enrgy_rate_sol);
}