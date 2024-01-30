#include<stdio.h>
#include<stdlib.h>
#include"type.h"
#include"scalar.h"
#include"b_matrix.h"
#include"d_matrix.h"

extern Global global;
extern Option option;

void generate_coefficient_matrix(){
    double ke_matrix[60][60];

    for(int point = 0; point < global.subdomain.N_point; point++){

    }
}

void generate_subdomain_coefficient_matrix(int point_n, double (*ke_matrix)[60]){
    double b_t_matrix[60][6];
    double d_matrix[6][6];
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



}