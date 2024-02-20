#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"type.h"
#include"b_matrix.h"
#include"d_matrix.h"
#include"Output.h"

extern Global global;
extern Option option;

void Output_Linear_strain_data(double *du){
    FILE *fp_strain;
    FILE *fp_stress;
    FILE *fp_u;
    double d_matrix[6][6];
    double b_t_matrix[60][6];
    double DB[6][60];
    double strain[6];
    double stress[6];
    double *strain_all;
    double *stress_all;

    if((strain_all = (double *)calloc(6 * global.subdomain.N_point, sizeof(double))) == NULL){
        printf("Error:Strain's memory is not enough\n");
        exit(-1);
    }
    if((stress_all = (double *)calloc(6 * global.subdomain.N_point, sizeof(double))) == NULL){
        printf("Error:Stress's memory is not enough\n");
        exit(-1);
    }

    fp_strain = fopen("Data_Files_Output/Output_strain.dat", "w");
    fprintf(fp_strain, "subdomain   /   Lagrange strain   xx            yy          zz          xy          yz          zx\n");
    
    for(int point = 0; point < global.subdomain.N_point; point++){
        int N_support = global.subdomain.support_offset[point + 1] - global.subdomain.support_offset[point];
        generate_linear_b_matrix(b_t_matrix, point);
       
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

        for(int i = 0; i < 6; i++)
            strain_all[6 * point + i] = strain[i];

        for(int i = 0; i < 6; i++){
            fprintf(fp_strain, "%+15.14e   ", strain[i]);
        }
        fprintf(fp_strain, "\n");
    }
    fclose(fp_strain);

    fp_u = fopen("Data_Files_Output/Output_displacement.dat", "w");
    fprintf(fp_u, "point        /displacement           x           y           z\n");
    for(int point = 0; point < global.subdomain.N_point; point++){
        fprintf(fp_u, "%5d      ", point);
        fprintf(fp_u, "%+15.14e %+15.14e %+15.14e\n", du[option.dim * point], du[option.dim * point + 1], du[option.dim * point + 2]);
    }
    fclose(fp_u);

    fp_stress = fopen("Data_Files_Output/Output_stress.dat", "w");
    fprintf(fp_stress, "subdomain   /   Cauchy stress   xx            yy          zz          xy          yz          zx\n");
    for(int point = 0; point < global.subdomain.N_point; point++){
        int N_support = global.subdomain.support_offset[point + 1] - global.subdomain.support_offset[point];
        generate_linear_b_matrix(b_t_matrix, point);
        generateElasticDMatrix(d_matrix);
        fprintf(fp_stress, "%5d      ", point);
        for(int i = 0; i < 6; i++){
            for(int j = 0; j < option.dim * (N_support + 1); j++){
                double DB_ij = 0.;
                for(int k = 0; k < 6; k++){
                    DB_ij += d_matrix[i][k] * b_t_matrix[j][k];
                }
                DB[i][j] = DB_ij;
            }
        }
        for(int i = 0; i < 6; i++){
            double stress_i = 0.;
            for(int j = 0; j < N_support; j++){
                for(int k = 0; k < option.dim; k++){
                    stress_i += DB[i][option.dim * (j + 1) + k] * du[option.dim * global.subdomain.support[global.subdomain.support_offset[point] + j] + k];
                }
            }
            for(int j = 0; j < option.dim; j++){
                stress_i += DB[i][j] * du[option.dim * point + j];
            }
            stress[i] = stress_i;
        }

        for(int i = 0; i < 6; i++)
            stress_all[6 * point + i] = stress[i];
        
        for(int i = 0; i < 6; i++)
            fprintf(fp_stress, "%+15.14e    ", stress[i]);
        fprintf(fp_stress, "\n");
    }
    fclose(fp_stress);

    #if 1
    FILE *fp_debug;
    FILE *fp_debug_u;
    fp_debug = fopen("Data_Files_Output/debug_stress.dat", "w");
    for(int point = 0; point < global.subdomain.N_point; point++){
        if(fabs(global.subdomain.point_XYZ[3*point+2] - 0.0) < 1.0e-5 && fabs(global.subdomain.point_XYZ[3*point]  - 0.0) < 1.0e-5){
            fprintf(fp_debug, "%5d      ", point);
            fprintf(fp_debug, "%+15.14e %+15.14e %+15.14e  ", global.subdomain.point_XYZ[option.dim * point], global.subdomain.point_XYZ[option.dim * point + 1], global.subdomain.point_XYZ[option.dim * point + 2]);
            fprintf(fp_debug, "%+15.14e %+15.14e %+15.14e %+15.14e %+15.14e %+15.14e\n", stress_all[6 * point], stress_all[6 * point + 1], stress_all[6 * point + 2], stress_all[6 * point + 3],  stress_all[6 * point + 4],  stress_all[6 * point + 5]);
        }
    }
    fclose(fp_debug);

    fp_debug_u = fopen("Data_Files_Output/debug_displacement.dat", "w");
    for(int point = 0; point < global.subdomain.N_point; point++){
        if(fabs(global.subdomain.point_XYZ[3*point+2] - 0.0) < 1.0e-5 && fabs(global.subdomain.point_XYZ[3*point + 1]  - 0) < 1.0e-5){
            fprintf(fp_debug_u, "%5d      ", point);
            fprintf(fp_debug_u, "%+15.14e %+15.14e %+15.14e  ", global.subdomain.point_XYZ[option.dim * point], global.subdomain.point_XYZ[option.dim * point + 1], global.subdomain.point_XYZ[option.dim * point + 2]);
            fprintf(fp_debug_u, "%+15.14e %+15.14e %+15.14e\n", du[3 * point], du[3 * point + 1], du[3 * point + 2]);
        }
    }
    fclose(fp_debug_u);
    #endif 

    free(stress_all);
    free(strain_all);

}