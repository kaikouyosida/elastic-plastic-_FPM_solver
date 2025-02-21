#pragma warning(disable: 4100) // 引数が未使用の場合
#pragma warning(disable: 4189) // ローカル変数が未使用の場合
#pragma warning(disable: 4996) //fopenの警告番号

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"type.h"
#include"scalar.h"
#include"vector.h"

extern Global global;
extern Option option;

#define NUMBER_OF_NODE_IN_FACE 4
#define NUMBER_OF_NODE_IN_SUBDOMAIN 8

void sortVector(double *vector, const int num){
    for (int i = 0; i < num - 1; i++)
        for (int j = num - 1; j > i; j--)
            if (vector[j - 1] > vector[j])
                swapReals(&vector[j - 1], &vector[j]);
}

void reverseVector(double *vector, const int num){
    for (int i = 0; i < num / 2; i++)
        swapReals(&vector[i], &vector[num - 1 - i]);
}

void assemble_vector(int point_n, double **global_vecter, const double *element_vector){
    int support[60];        //サポート番号
    int N_support = global.subdomain.support_offset[point_n + 1] - global.subdomain.support_offset[point_n]; //サポート点の数
   
    for(int i = 0; i < N_support; i++)
        support[i] = global.subdomain.support[global.subdomain.support_offset[point_n] + i];
    
    for(int i = 0; i < N_support; i++){
        for(int j = 0; j < option.dim; j++){
            global_vecter[support[i]][j] += element_vector[option.dim * (i + 1) + j];
        }
    }
    for(int i = 0; i < option.dim; i++)
        global_vecter[point_n][i] += element_vector[i];
}

//頂点2から頂点1に伸びるベクトルの生成
void generate_current_node_vector(int node_1, int node_2, double *vector){
    for(int i = 0; i < option.dim; i++)
        vector[i] = global.subdomain.node_XYZ[option.dim * node_1 + i] - global.subdomain.node_XYZ[option.dim * node_2 + i]
                  + global.subdomain.nodal_displacements[node_1][i] - global.subdomain.nodal_displacements[node_2][i];
}

//頂点から任意の座標に伸びるベクトルの生成
void generate_current_node_to_point_vector(int node,double *point,double *vector){
    for(int i = 0; i < option.dim; i++)
        vector[i] = point[i] - global.subdomain.node_XYZ[option.dim * node + i]
                - global.subdomain.nodal_displacements[node][i];
}

//頂点2から頂点1に伸びるベクトルの生成
void generate_current_edge_vector(double *vector ,int point_n, int node_id_1, int node_id_2, int *subdomain_node){
    for(int i = 0; i < option.dim; i++)
        vector[i] = global.subdomain.node_XYZ[option.dim * subdomain_node[node_id_1] + i] - global.subdomain.node_XYZ[option.dim * subdomain_node[node_id_2] + i]
                  + global.subdomain.nodal_displacement_sd[point_n][node_id_1][i] - global.subdomain.nodal_displacement_sd[point_n][node_id_2][i];
}

//頂点から任意の点の座標に伸びるベクトルの生成
void generate_current_points_vector(double *vector, double *point_xyz, int point_n, int node_id,  int *subdomain_node){
    for(int i = 0; i < option.dim; i++)
        vector[i] = point_xyz[i]
                    - (global.subdomain.node_XYZ[option.dim * subdomain_node[node_id] + i]
                    + global.subdomain.nodal_displacement_sd[point_n][node_id][i]);
}

//point側のサブドメインにおける形状関数から得た節点の現在座標
void generate_current_face_node(double face_node_XYZ[4][3], int node_id[4] ,int subdomain_node[8], int point){
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < option.dim; j++){
            face_node_XYZ[i][j] = global.subdomain.node_XYZ[option.dim * subdomain_node[node_id[i]] + j]
                                + global.subdomain.nodal_displacement_sd[point][node_id[i]][j];
        }
    }
}

//point側のサブドメインにおける形状関数から得た節点の現在座標
void generate_current_node_of_face(double face_node_XYZ[4][3], int face_n, int point){
    int subdomain_node[8];
    int node_id[4];

    //各サブドメインが持つ節点番号を格納する
    generate_subdomain_node(point, subdomain_node);

    //内部境界面のノード番号のアドレスを格納
    generate_node_id(face_n, point, subdomain_node, node_id);
    for(int i = 0; i < NUMBER_OF_NODE_IN_FACE; i++){
        for(int j = 0; j < option.dim; j++){
            face_node_XYZ[i][j] = global.subdomain.node_XYZ[option.dim * subdomain_node[node_id[i]] + j]
                                + global.subdomain.nodal_displacement_sd[point][node_id[i]][j];
        }
    }
}

//サブドメインがもつ節点の番号を格納する
void generate_subdomain_node(int point_n, int subdomain_node[8]){
    for(int i = 0; i < NUMBER_OF_NODE_IN_SUBDOMAIN; i++){
        subdomain_node[i] = global.subdomain.subdomain_node[NUMBER_OF_NODE_IN_SUBDOMAIN * point_n + i];
    }
}

//内部境界面のノード番号のアドレスを格納
void generate_node_id(int face_n, int point_n, int subdomain_node[8], int node_id[4]){
    int ref_number = global.subdomain.vertex_offset[face_n];
    for(int i = 0; i < NUMBER_OF_NODE_IN_FACE; i++)
        for(int j = 0; j < NUMBER_OF_NODE_IN_SUBDOMAIN; j++)
            if(subdomain_node[j] == global.subdomain.node[ref_number + i])
                node_id[i] = j;  
}

//物理座標におけるガウス点の座標を計算
void generate_gauss_point_coordinate(const int s, const int t, const double face_node_XYZ[4][3], const double *X, double xyz[3]){
    for(int i = 0; i < option.dim; i++)
        xyz[i] = 0.25 * (1.0 - X[s]) * (1.0 - X[t]) * face_node_XYZ[0][i]
                + 0.25 * (1.0 - X[s]) * (1.0 + X[t]) * face_node_XYZ[1][i]
                + 0.25 * (1.0 + X[s]) * (1.0 + X[t]) * face_node_XYZ[2][i]
                + 0.25 * (1.0 + X[s]) * (1.0 - X[t]) * face_node_XYZ[3][i];
}