#include<stdio.h>
#include<stdlib.h>
#include"type.h"
#include"b_matrix.h"
#include"Output.h"

extern Global global;
extern Option option;

void Output_Linear_strain_data(double *du){
    FILE *fp_strain;
    FILE *fp_u;
    double b_t_matrix[60][6];
    double strain[6];

    fp_strain = fopen("linear_strain.dat", "w");
    fprintf(fp_strain, "subdomain   /   Lagrange strain   xx            yy          zz          xy          yz          zx\n");
    printf("status3\n");
    for(int point = 0; point < global.subdomain.N_point; point++){
        int N_support = global.subdomain.support_offset[point + 1] - global.subdomain.support_offset[point];
        generate_linear_b_matrix(b_t_matrix, point);
        printf("status\n");

        fprintf(fp_strain, "%5d         ", point);
        
        for(int i = 0; i < 6; i++){
            double strain_i = 0.;
            for(int j = 0; j < N_support; j++){
                for(int k = 0; k < option.dim; k++){
                    strain_i += b_t_matrix[option.dim * (j + 1) + k][i] * du[option.dim * (global.subdomain.support[global.subdomain.support_offset[point] + j]) + k];
                }
            }
            for(int j = 0; j < option.dim; j++)
                strain_i += b_t_matrix[j][i] * du[option.dim * point + j];
            
            strain[i] = strain_i;
        }

        for(int i = 0; i < 6; i++){
            
            fprintf(fp_strain, "%+15.14e   ", strain[i]);
        }
        fprintf(fp_strain, "\n");
    }
    fclose(fp_strain);

    fp_u = fopen("Output_displacement.dat", "w");
    fprintf(fp_u, "point        /displacement           x           y           z\n");
    for(int point = 0; point < global.subdomain.N_point; point++){
        fprintf(fp_u, "%5d      ", point);
        fprintf(fp_u, "%+15.14e %+15.14e %+15.14e\n", du[option.dim * point], du[option.dim * point + 1], du[option.dim * point + 2]);
    }
    fclose(fp_u);
}