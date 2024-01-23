#include<stdio.h>
#include<stdlib.h>

#include"type.h"
#include"matrix.h"

extern Global global;
extern Option option;

void calc_G(int dim, int point_n, double *point_xyz, int *support_offset, int *support, double **G)
{
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
    if (1 < N_support)
    {
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
    {
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

    for (int i = 0; i < dim; i++)
    {
        for (int j = 0; j < dim; j++)
        {
            for (int k = 0; k < N_support + 1; k++)
                G[dim * i + j][dim * k + i] = G1[j][k];
        }
    }

    free_matrix(G1);
    free_matrix(A);
}

void generate_linear_b_matrix(double (*b_t_matrix)[6], int point_n){
    double **G;
    double N_support = global.subdomain.support_offset[point_n + 1] - global.subdomain.support_offset[point_n];
    double *latest_point_XYZ;

    if((latest_point_XYZ = (double *)calloc(option.dim * global.subdomain.N_point, sizeof(double))) == NULL){
        exit(-1);
    }
    G = matrix(option.dim * option.dim, option.dim * N_support + option.dim);


    for(int i = 0; i < global.subdomain.N_point; i++){
        for(int j = 0; j < option.dim; j++){
            latest_point_XYZ[option.dim * i + j] = global.subdomain.point_XYZ[option.dim * i + j]
                    + global.subdomain.displacement[i][j]
                    + global.subdomain.displacement_increment[i][j];
        }
    }

    calc_G(option.dim, point_n, latest_point_XYZ, global.subdomain.support_offset, global.subdomain.support, G);
    
        double dn_dx_0_i = G[0][0];
        double dn_dx_1_i = G[1][0];
        double dn_dx_2_i = G[2][0];
        
        b_t_matrix[0][0] = dn_dx_0_i;
        b_t_matrix[0][1] = 0.0;
        b_t_matrix[0][2] = 0.0;
        b_t_matrix[0][3] = dn_dx_1_i;
        b_t_matrix[0][4] = 0.0;
        b_t_matrix[0][5] = dn_dx_2_i;

        b_t_matrix[1][0] = 0.0;
        b_t_matrix[1][1] = dn_dx_1_i;
        b_t_matrix[1][2] = 0.0;
        b_t_matrix[1][3] = dn_dx_0_i;
        b_t_matrix[1][4] = dn_dx_2_i;
        b_t_matrix[1][5] = 0.0;

        b_t_matrix[2][0] = 0.0;
        b_t_matrix[2][1] = 0.0;
        b_t_matrix[2][2] = dn_dx_2_i;
        b_t_matrix[2][3] = 0.0;
        b_t_matrix[2][4] = dn_dx_1_i;
        b_t_matrix[2][5] = dn_dx_0_i;

        for(int i = 1; i < N_support + 1; i++){
            double dn_dx_0_i = G[0][3 * i];
            double dn_dx_1_i = G[1][3 * i + 1];
            double dn_dx_2_i = G[2][3 * i + 2];

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
    free(latest_point_XYZ);
}