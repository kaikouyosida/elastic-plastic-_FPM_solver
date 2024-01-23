#include<stdio.h>
#include<stdlib.h>
#include"type.h"
extern Global global;
extern Option option;

void generateElasticDMatrix(double (*d_matrix)[6]){
    double young_modulus
        = global.material.E_mod;
    double poisson_ratio
        = global.material.nu_mod;


    double d00
        = young_modulus * (1.0 - poisson_ratio)
        / (1.0 + poisson_ratio) / (1.0 - 2.0 * poisson_ratio);
    double d01
        = young_modulus * poisson_ratio
        / (1.0 + poisson_ratio) / (1.0 - 2.0 * poisson_ratio);
    double d33
        = 0.5 * young_modulus / (1.0 + poisson_ratio);

    d_matrix[0][0] = d00; d_matrix[0][1] = d01; d_matrix[0][2] = d01;
    d_matrix[0][3] = 0.0; d_matrix[0][4] = 0.0; d_matrix[0][5] = 0.0;

    d_matrix[1][0] = d01; d_matrix[1][1] = d00; d_matrix[1][2] = d01;
    d_matrix[1][3] = 0.0; d_matrix[1][4] = 0.0; d_matrix[1][5] = 0.0;

    d_matrix[2][0] = d01; d_matrix[2][1] = d01; d_matrix[2][2] = d00;
    d_matrix[2][3] = 0.0; d_matrix[2][4] = 0.0; d_matrix[2][5] = 0.0;

    d_matrix[3][0] = 0.0; d_matrix[3][1] = 0.0; d_matrix[3][2] = 0.0;
    d_matrix[3][3] = d33; d_matrix[3][4] = 0.0; d_matrix[3][5] = 0.0;

    d_matrix[4][0] = 0.0; d_matrix[4][1] = 0.0; d_matrix[4][2] = 0.0;
    d_matrix[4][3] = 0.0; d_matrix[4][4] = d33; d_matrix[4][5] = 0.0;

    d_matrix[5][0] = 0.0; d_matrix[5][1] = 0.0; d_matrix[5][2] = 0.0;
    d_matrix[5][3] = 0.0; d_matrix[5][4] = 0.0; d_matrix[5][5] = d33;
}