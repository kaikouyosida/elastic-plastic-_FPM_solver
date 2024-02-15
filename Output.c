#include<stdio.h>
#include<stdlib.h>
#include"type.h"
#include"b_matrix.h"
#include"Output.h"

extern Global global;
extern Option option;

void Output_Linear_data(double *du){
    Output_Linear_strain_data(du);
}

void Output_Linear_strain_data(double *du){
    FILE *fp_strain;
    FILE *fp_debug;
    double b_t_matrix[60][6];
    double strain[6];
    fp_debug = fopen("debug.dat", "w");
    fp_strain = fopen("Data_Files_Output/linear_strain.dat", "w");
    fprintf(fp_strain, "subdomain   /   Lagrange strain   xx            yy          zz          xy          yz          zx\n");
    for(int point = 0; point < global.subdomain.N_point; point++){
        int N_support = global.subdomain.support_offset[point + 1] - global.subdomain.support_offset[point];
        generate_linear_b_matrix(b_t_matrix, point);

        for(int i = 0; i < (N_support+1) * option.dim; i++){
            for(int j = 0; j < 6; j++){
                fprintf(fp_debug, "%+4.3e   ", b_t_matrix[i][j]);
            }
            fprintf(fp_debug, "\n");
        }
        fprintf(fp_debug, "\n");


        fprintf(fp_strain, "%5d         ", point);
        
        for(int i = 0; i < 6; i++){
            double strain_i = 0.;
            for(int j = 0; j < N_support; j++){
                for(int k = 0; k < option.dim; k++){
                    strain_i += b_t_matrix[option.dim * (j + 1) + k][i] * du[option.dim * (global.subdomain.support[global.subdomain.support_offset[point] + j]) + k];
                }
            }
            for(int j = 0; j < option.dim; j++)
                strain_i += b_t_matrix[j][i] * global.subdomain.displacement[point][j];
            
            strain[i] = strain_i;
        }

        for(int i = 0; i < 6; i++){
            
            fprintf(fp_strain, "%+15.14e   ", strain[i]);
        }
        fprintf(fp_strain, "\n");
    }
    fclose(fp_debug);
    fclose(fp_strain);
}