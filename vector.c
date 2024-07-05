#include<stdio.h>
#include<math.h>

#include"type.h"
#include"scalar.h"
#include"vector.h"

extern Global global;
extern Option option;

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

void assemble_vector(int point_n, double **global_vecter, double *element_vector){
    int support[60];        //サポート番号
    int N_support = global.subdomain.support_offset[point_n + 1] - global.subdomain.support_offset[point_n]; 

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
                + global.subdomain.nodal_displacements[node_1][i] - global.subdomain.nodal_displacements[node_2][i]
                + global.subdomain.nodal_displacement_increments[node_1][i] - global.subdomain.nodal_displacement_increments[node_2][i];
}

//頂点から任意の座標に伸びるベクトルの生成
void generate_current_node_to_point_vector(int node,double *point,double *vector){
    for(int i = 0; i < option.dim; i++)
        vector[i] = point[i] - global.subdomain.node_XYZ[option.dim * node + i]
                - global.subdomain.nodal_displacements[node][i]
                - global.subdomain.nodal_displacement_increments[node][i];
}