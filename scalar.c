#include<stdio.h>
#include<math.h>
#include"scalar.h"
#include"tensor.h"
#include"type.h"
#include"vector.h"

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
    int N_face = global.subdomain.face_offset[point_n + 1] - global.subdomain.face_offset[point_n];
    int subdomain_node[60];              //サブドメインpoint_nがもつ頂点番号
    double edge1[3],edge2[3],edge3[3];   //辺のベクトル
    double edge1_cross_edge2[3];         //辺のベクトルの外積

    //6面体（節点8個）の重心を算術平均で計算
    for(int i = 0; i < 8; i++)  
        subdomain_node[i] = global.subdomain.subdomain_node[8 * point_n + i];
   
    for(int i = 0; i < option.dim; i++){
        double center_i = 0.;
        for(int j = 0; j < 8; j++){
            center_i += (global.subdomain.node_XYZ[option.dim * subdomain_node[j] + i]
                     + global.subdomain.nodal_displacements[subdomain_node[j]][i]
                     + global.subdomain.nodal_displacement_increments[subdomain_node[j]][i]);
        }
        center[i] = center_i / 8.0;
    }

    //各面の頂点と重心の5点で四角すいを形成。その後各面で総和をとって六面体の面積を計算する
    for(int i = 0; i < N_face; i++){
        int ref_num = global.subdomain.vertex_offset[global.subdomain.face[global.subdomain.face_offset[point_n] + i]];


        generate_current_node_vector(global.subdomain.node[ref_num], global.subdomain.node[ref_num + 1], edge1);
        generate_current_node_vector(global.subdomain.node[ref_num], global.subdomain.node[ref_num + 2], edge2);
        generate_current_node_to_point_vector(global.subdomain.node[ref_num], center, edge3);
        
        cross_product(option.dim, edge1, edge2, edge1_cross_edge2);
        volume += fabs(dot_product(option.dim,edge1_cross_edge2, edge3));
        
        generate_current_node_vector(global.subdomain.node[ref_num + 2], global.subdomain.node[ref_num + 3], edge1);
        generate_current_node_vector(global.subdomain.node[ref_num + 2], global.subdomain.node[ref_num + 1], edge2);
        generate_current_node_to_point_vector(global.subdomain.node[ref_num], center, edge3);
        cross_product(option.dim, edge1, edge2, edge1_cross_edge2);

        volume += fabs(dot_product(option.dim,edge1_cross_edge2, edge3));
    }
    volume /= 6.0;
    
    return volume;
}

//物理空間座標→正規化座標に変換するためのスカラー値を計算
double calc_area_change(int face_n, int s, int t, double *X){
    int ref_num = global.subdomain.vertex_offset[face_n];
    double scalar = 0.;
    double area_vector[3];
    double dx_ds[3];
    double dx_dt[3];
    double x1[3];
    double x2[3];
    double x3[3];
    double x4[3];
    int node1 = global.subdomain.node[ref_num];
    int node2 = global.subdomain.node[ref_num + 1];
    int node3 = global.subdomain.node[ref_num + 2];
    int node4 = global.subdomain.node[ref_num + 3];

    for(int i = 0; i < option.dim; i++){
        x1[i] = global.subdomain.node_XYZ[option.dim * node1 + i] + global.subdomain.nodal_displacements[node1][i] + global.subdomain.nodal_displacement_increments[node1][i];
        x2[i] = global.subdomain.node_XYZ[option.dim * node2 + i] + global.subdomain.nodal_displacements[node2][i] + global.subdomain.nodal_displacement_increments[node2][i];
        x3[i] = global.subdomain.node_XYZ[option.dim * node3 + i] + global.subdomain.nodal_displacements[node3][i] + global.subdomain.nodal_displacement_increments[node3][i];
        x4[i] = global.subdomain.node_XYZ[option.dim * node4 + i] + global.subdomain.nodal_displacements[node4][i] + global.subdomain.nodal_displacement_increments[node4][i];
    }

    //x_h = 1/4*(1-ξ*ξ_i)*(1-η*η_i)*x_iのξ, ηに対する微係数を計算.
    for(int i = 0; i < option.dim ; i++){
        dx_ds[i] = -0.25 * (1.0 - X[t]) * x1[i] - 0.25 * (1.0 + X[t]) * x2[i] + 0.25 * (1.0 + X[t]) * x3[i] + 0.25 * (1.0 - X[t]) * x4[i];
        dx_dt[i] = -0.25 * (1.0 - X[s]) * x1[i] + 0.25 * (1.0 - X[s]) * x2[i] + 0.25 * (1.0 + X[s]) * x3[i] - 0.25 * (1.0 + X[s]) * x4[i];
    }
    
    cross_product(option.dim, dx_ds, dx_dt, area_vector);
    scalar = norm(area_vector, option.dim);
    return scalar;
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