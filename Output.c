#pragma warning(disable: 4100) // 引数が未使用の場合
#pragma warning(disable: 4189) // ローカル変数が未使用の場合
#pragma warning(disable: 4996) //fopenの警告番号

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"type.h"
#include"b_matrix.h"
#include"d_matrix.h"
#include"matrix.h"
#include"Output.h"
#include"scalar.h"

extern Global global;
extern Option option;

void Output_Linear_strain_data(double *du){
    FILE *fp_strain;
    FILE *fp_stress;
    FILE *fp_u;
    double d_matrix[6][6];
    double b_t_matrix[60][6];
    int support[60];
    double deformation_gradient[9];
    double DB[6][60];
    double strain[6];
    double stress[6];
    double **G;
    double *strain_all;
    double *stress_all;
    double strain_energy;

    if((strain_all = (double *)calloc(6 * global.subdomain.N_point, sizeof(double))) == NULL){
        printf("Error:Strain's memory is not enough\n");
        exit(-1);
    }
    if((stress_all = (double *)calloc(6 * global.subdomain.N_point, sizeof(double))) == NULL){
        printf("Error:Stress's memory is not enough\n");
        exit(-1);
    }
    printf("stati");
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

        //for(int i = 0; i < 6; i++){
            //for(int j = 0; j < 6; j++){
                //printf("%+5.4e  ", d_matrix[i][j]);
            //}
            //printf("\n");
        //}
        //exit(-1);

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

    #if 0
    FILE *fp_debug;
    FILE *fp_debug_u;
    fp_debug = fopen("Data_Files_Output/debug_stress.dat", "w");
    for(int point = 0; point < global.subdomain.N_point; point++){
        if(fabs(global.subdomain.point_XYZ[3*point+2] - 1.0) < 1.0e-5 && fabs(global.subdomain.point_XYZ[3*point]  - 0.0) < 1.0e-5){
            fprintf(fp_debug, "%5d      ", point);
            fprintf(fp_debug, "%+15.14e %+15.14e %+15.14e  ", global.subdomain.point_XYZ[option.dim * point], global.subdomain.point_XYZ[option.dim * point + 1], global.subdomain.point_XYZ[option.dim * point + 2]);
            fprintf(fp_debug, "%+15.14e %+15.14e %+15.14e %+15.14e %+15.14e %+15.14e\n", stress_all[6 * point], stress_all[6 * point + 1], stress_all[6 * point + 2], stress_all[6 * point + 3],  stress_all[6 * point + 4],  stress_all[6 * point + 5]);
        }
    }
    fclose(fp_debug);
    #endif 
    //FILE *fp_debug_u;
    //fp_debug_u = fopen("Data_Files_Output/debug_displacement.dat", "w");
    //for(int point = 0; point < global.subdomain.N_point; point++){
        //if(fabs(global.subdomain.point_XYZ[3*point] - 0.0) < 1.0e-5 && fabs(global.subdomain.point_XYZ[3*point + 1]  -0.0) < 1.0e-5){
            //fprintf(fp_debug_u, "%5d      ", point);
            //fprintf(fp_debug_u, "%+15.14e %+15.14e %+15.14e  ", global.subdomain.point_XYZ[option.dim * point], global.subdomain.point_XYZ[option.dim * point + 1], global.subdomain.point_XYZ[option.dim * point + 2]);
            //fprintf(fp_debug_u, "%+15.14e %+15.14e %+15.14e\n", du[3 * point], du[3 * point + 1], du[3 * point + 2]);
        //}
    //}
    //fclose(fp_debug_u);
     
    #if 0
    for(int point = 0; point < global.subdomain.N_point; point++){
        int N_support = global.subdomain.support_offset[point + 1] - global.subdomain.support_offset[point];
        for(int i = 0; i < N_support; i++)
            support[i] = global.subdomain.support[global.subdomain.support_offset[point] + i];
        G = matrix(9, 3*(N_support + 1));
        calc_G(option.dim, point, global.subdomain.point_XYZ, global.subdomain.support_offset, global.subdomain.support, G);
        //for(int i = 0; i < 9; i++){
            //for(int j = 0; j < 3*(N_support+1); j++){
                //printf("%+5.4e  ", G[i][j]);
            //}
            //printf("\n");
        //}
        //exit(-1);
        for(int i = 0; i < 9; i++){
            double deformation_gradient_i = 0.;
            for(int j = 0; j < N_support; j++){
                for(int k = 0; k < option.dim; k++){
                    deformation_gradient_i += G[i][3 * (j + 1) + k] * du[3 * support[j] + k];
                }
            }
            for(int j = 0;  j < option.dim; j++)
                deformation_gradient_i += G[i][j] * du[option.dim * point + j];
            deformation_gradient[i] = deformation_gradient_i;
        }
        
        free_matrix(G);
        printf("%5d  ", point);
        for(int i = 0; i < 9; i++)
            printf("%+15.14e ", deformation_gradient[i]);
        printf("\n");
    }
    #endif

    #if 0
    for(int point = 0; point < global.subdomain.N_point; point++){
        printf("%5d  ", point);
        for(int i = 0; i < 6; i++){
            double stress_i = 0;
            for(int j = 0; j < 6; j++){
                stress_i += d_matrix[i][j] * strain_all[6 * point + j];
            }
            stress[i] = stress_i;
        }
        for(int i = 0; i < 6; i++){
            printf("%+8.7e   ", stress[i]);
        }
        printf("\n");
    }
    #endif
    #if 0
    double internal_force[12000000];
    for(int i = 0; i < option.dim*global.subdomain.N_point; i++){
        double internal_i = 0; 
        for(int j = 0; j < option.dim*global.subdomain.N_point; j++){
            internal_i += global.subdomain.Global_K[3*global.subdomain.N_point * i + j] * du[j];
        }
        internal_force[i] = internal_i;
    }
    for(int i = 0; i < global.subdomain.N_point; i++){
        printf("%5d  ", i);
        for(int j = 0; j < option.dim; j++){
            printf("%+15.14e ", internal_force[3*i+j]-global.subdomain.global_residual_force[3*i+j]);
        }
        printf("\n");
    }
    #endif

    free(stress_all);
    free(strain_all);

}
