#include<stdio.h>
#include<stdlib.h>

#include"type.h"
#include"matrix.h"
#include"b_matrix.h"

extern Global global;
extern Option option;

void calc_G(const int dim, const int point_n, const double *point_xyz, const int *support_offset, const int *support, double **G){
    double **A;
    double **AT;
    double **ATA;
    double **inv_ATA;
    double **invATA_AT;
    double **G1;
    double L2 = 0.;                                                        // 2点間の距離 (サポートドメインが1つのとき)
    int N_support = support_offset[point_n + 1] - support_offset[point_n]; // サポートドメインの数
    int support_point = 0;                                                 // サポートドメインの番号

    A = matrix(N_support, dim);
    for (int i = 0; i < N_support; i++)
    {
        support_point = support[support_offset[point_n] + i];
        for (int j = 0; j < dim; j++)
            A[i][j] = point_xyz[dim * support_point + j] - point_xyz[dim * point_n + j];
    }

    G1 = matrix(dim, N_support + 1);
    // サポートドメインが2つ以上 //
    if (1 < N_support){

        AT = matrix(dim, N_support);
        ATA = matrix(dim, dim);
        inv_ATA = matrix(dim, dim);
        invATA_AT = matrix(dim, N_support);

        trans_mat(N_support, dim, A, AT);
        multi_mat(dim, dim, AT, A, N_support, ATA);
        inverse_mat(dim, ATA, inv_ATA);
        multi_mat(dim, N_support, inv_ATA, AT, dim, invATA_AT);

        for (int i = 0; i < dim; i++)
        {
            for (int j = 0; j < N_support; j++)
                G1[i][0] -= invATA_AT[i][j];
            for (int j = 0; j < N_support; j++)
                G1[i][j + 1] = invATA_AT[i][j];
        }

        free_matrix(invATA_AT);
        free_matrix(inv_ATA);
        free_matrix(ATA);
        free_matrix(AT);
    }
    // サポートドメインが1つの場合 //
    else
    {   printf("support is too small!\n");
        support_point = support[support_offset[point_n]];
        for (int i = 0; i < dim; i++)
            L2 += point_xyz[dim * point_n + i] * point_xyz[dim * point_n + i] + point_xyz[dim * support_point + i] * point_xyz[dim * support_point + i];
        for (int i = 0; i < dim; i++)
            L2 -= 2.0 * point_xyz[dim * point_n + i] * point_xyz[dim * support_point + i];

        for (int i = 0; i < dim; i++)
        {
            G1[i][0] = -1.0 * A[0][i] / L2;
            G1[i][1] = A[0][i] / L2;
        }
    }

    for (int i = 0; i < dim; i++){
        for (int j = 0; j < dim; j++){
            for (int k = 0; k < N_support + 1; k++)
                G[dim * i + j][dim * k + i] = G1[j][k];
        }
    }

    free_matrix(G1);
    free_matrix(A);
}

void generate_linear_b_matrix(double (*b_t_matrix)[6], const int point_n){
    double **G;
    int N_support = global.subdomain.support_offset[point_n + 1] - global.subdomain.support_offset[point_n];
    
    G = matrix(option.dim * option.dim, option.dim * (N_support + 1));

    if(option.solver_type == 1){
        double *current_point_XYZ;

        if((current_point_XYZ = (double *)calloc(option.dim * global.subdomain.N_point, sizeof(double))) == NULL){
            printf("Error: current_point_XYZ's memory is not enough\n");
            exit(-1);
        }
        for(int i = 0; i < global.subdomain.N_point; i++){
            for(int j = 0; j < option.dim; j++){
                current_point_XYZ[option.dim * i + j] = global.subdomain.point_XYZ[option.dim * i + j]
                                                        + global.subdomain.displacement[i][j]
                                                        + global.subdomain.displacement_increment[i][j];
            }
        }

        calc_G(option.dim, point_n, current_point_XYZ, global.subdomain.support_offset, global.subdomain.support, G);

        free(current_point_XYZ);
    }else{
        calc_G(option.dim, point_n, global.subdomain.point_XYZ, global.subdomain.support_offset, global.subdomain.support, G);
    }
        const double dn_dx_0_0 = G[0][0];
        const double dn_dx_1_0 = G[1][0];
        const double dn_dx_2_0 = G[2][0];
        
        b_t_matrix[0][0] = dn_dx_0_0;
        b_t_matrix[0][1] = 0.0;
        b_t_matrix[0][2] = 0.0;
        b_t_matrix[0][3] = dn_dx_1_0;
        b_t_matrix[0][4] = 0.0;
        b_t_matrix[0][5] = dn_dx_2_0;

        b_t_matrix[1][0] = 0.0;
        b_t_matrix[1][1] = dn_dx_1_0;
        b_t_matrix[1][2] = 0.0;
        b_t_matrix[1][3] = dn_dx_0_0;
        b_t_matrix[1][4] = dn_dx_2_0;
        b_t_matrix[1][5] = 0.0;

        b_t_matrix[2][0] = 0.0;
        b_t_matrix[2][1] = 0.0;
        b_t_matrix[2][2] = dn_dx_2_0;
        b_t_matrix[2][3] = 0.0;
        b_t_matrix[2][4] = dn_dx_1_0;
        b_t_matrix[2][5] = dn_dx_0_0;

        for(int i = 1; i < N_support + 1; i++){
            const double dn_dx_0_i = G[0][3 * i];
            const double dn_dx_1_i = G[1][3 * i];
            const double dn_dx_2_i = G[2][3 * i];

            b_t_matrix[option.dim * i][0] = dn_dx_0_i;
            b_t_matrix[option.dim * i][1] = 0.0;
            b_t_matrix[option.dim * i][2] = 0.0;
            b_t_matrix[option.dim * i][3] = dn_dx_1_i;
            b_t_matrix[option.dim * i][4] = 0.0;
            b_t_matrix[option.dim * i][5] = dn_dx_2_i;

            b_t_matrix[option.dim * i + 1][0] = 0.0;
            b_t_matrix[option.dim * i + 1][1] = dn_dx_1_i;
            b_t_matrix[option.dim * i + 1][2] = 0.0;
            b_t_matrix[option.dim * i + 1][3] = dn_dx_0_i;
            b_t_matrix[option.dim * i + 1][4] = dn_dx_2_i;
            b_t_matrix[option.dim * i + 1][5] = 0.0;

            b_t_matrix[option.dim * i + 2][0] = 0.0;
            b_t_matrix[option.dim * i + 2][1] = 0.0;
            b_t_matrix[option.dim * i + 2][2] = dn_dx_2_i;
            b_t_matrix[option.dim * i + 2][3] = 0.0;
            b_t_matrix[option.dim * i + 2][4] = dn_dx_1_i;
            b_t_matrix[option.dim * i + 2][5] = dn_dx_0_i;
        }

    free_matrix(G);
}

//Bマトリクスの計算
double generate_nonlinear_b_matrix(double (*b_t_matrix)[9], const int point_n){
    double **G;
    int N_support = global.subdomain.support_offset[point_n + 1] - global.subdomain.support_offset[point_n];
    double *latest_point_XYZ;
    G = matrix(option.dim * option.dim, option.dim * (N_support + 1));

    if((latest_point_XYZ = (double *)calloc(option.dim * global.subdomain.N_point, sizeof(double))) == NULL){
        exit(-1);
    }
    for(int i = 0; i < global.subdomain.N_point; i++){
        for(int j = 0; j < option.dim; j++){
            latest_point_XYZ[option.dim * i + j] = global.subdomain.point_XYZ[option.dim * i + j]
                                                    + global.subdomain.displacement[i][j]
                                                    + global.subdomain.displacement_increment[i][j];
        }
    }

    calc_G(option.dim, point_n, latest_point_XYZ, global.subdomain.support_offset, global.subdomain.support, G);

    for (int i = 0; i < N_support + 1; i++)
        for (int j = 0; j < option.dim; j++)
        {
            double dn_dx_j_i = G[j][option.dim * i];

            b_t_matrix[3 * i    ][3 * j]     = dn_dx_j_i;
            b_t_matrix[3 * i    ][3 * j + 1] = 0.0;
            b_t_matrix[3 * i    ][3 * j + 2] = 0.0;

            b_t_matrix[3 * i + 1][3 * j]     = 0.0;
            b_t_matrix[3 * i + 1][3 * j + 1] = dn_dx_j_i;
            b_t_matrix[3 * i + 1][3 * j + 2] = 0.0;

            b_t_matrix[3 * i + 2][3 * j]     = 0.0;
            b_t_matrix[3 * i + 2][3 * j + 1] = 0.0;
            b_t_matrix[3 * i + 2][3 * j + 2] = dn_dx_j_i;
        }
}

// 形状関数を計算//
void calc_shape(const double *xyz, const int dim, const int point_n, const double *point_xyz, const int *support_offset, double (*shapeF_t)[3]){
    const int N_support = support_offset[point_n + 1] - support_offset[point_n]; // サポートドメインの数
    double h[3];            //ポイント間距離
    double **G;             //変位勾配マトリクス
    double **shapeF;        //形状関数

    G = matrix(option.dim * option.dim, option.dim * (N_support + 1));
    shapeF = matrix(option.dim, option.dim * (N_support + 1));
    calc_G(dim, point_n, point_xyz, global.subdomain.support_offset, global.subdomain.support, G);

    for (int i = 0; i < dim; i++)
        h[i] = xyz[i] - point_xyz[dim * point_n + i];

    for (int i = 0; i < dim; i++){
        for (int j = 0; j < N_support + 1; j++){
            for(int k = 0; k < option.dim; k++){
                shapeF[i][dim * j + i] += h[k] * G[dim * i + k][dim * j + i];
            }
        }
    }
    for (int i = 0; i < dim; i++)
        shapeF[i][i] += 1.0;
    
    for(int i = 0; i < dim * (N_support + 1); i++){
        for(int j = 0; j < dim; j++){
            shapeF_t[i][j] = shapeF[j][i];
        }
    }
    free_matrix(shapeF);
    free_matrix(G);
}

//試行関数の計算 (pm = 0 or 1で変位と変位増分の試行関数を計算)
void trial_u(const double *xyz, const int point_n, const double *point_XYZ, double *u_h, const int pm){
    int N_support = global.subdomain.support_offset[point_n + 1] - global.subdomain.support_offset[point_n];
    double *u;
    double shapeF_t[60][3];

    if((u = (double *)calloc(option.dim * global.subdomain.N_point, sizeof(double))) == NULL){
        printf("Error:u's memory is not enough\n");
        exit(-1);
    }

    if(pm == 0){
        for(int i = 0; i < global.subdomain.N_point; i++){
            for(int j = 0; j < option.dim; j++){
                u[option.dim * i + j] = global.subdomain.displacement[i][j] + global.subdomain.displacement_increment[i][j];
            }
        }
    }else if(pm == 1){
        for(int i = 0; i < global.subdomain.N_point; i++){
            for(int j = 0; j < option.dim; j++){
                u[option.dim * i + j] = global.subdomain.displacement_increment[i][j];
            }
        }
    }else{
        printf("Undefined number.");
    }

    calc_shape(xyz, option.dim, point_n, point_XYZ, global.subdomain.support_offset, shapeF_t);

    for(int i = 0; i < option.dim; i++)
        u_h[i] = shapeF_t[i][i] * u[option.dim * point_n + i];

    for(int i = 0; i < option.dim; i++){
        double u_i = 0.;
        for(int j = 0; j < N_support; j++){
            u_i += shapeF_t[option.dim * (j + 1) + i][i] * u[option.dim * global.subdomain.support[global.subdomain.support_offset[point_n] + j] + i];
        }
        u_h[i] += u_i;
    } 

    free(u);
}