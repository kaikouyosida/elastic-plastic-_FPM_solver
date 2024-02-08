#include<stdio.h>
#include<stdlib.h>
#include"type.h"
#include"fpm.h"
#include"field.h"
#include"coefficient_matrix.h"
#include"internal_force.h"
#include"external_force.h"
#include"ImposeDirichretCondition.h"
#include"LU_decomposition.h"

extern Global global;
extern Option option;

void analize_by_NewtonRapdon(){
    double time;
    double time_increment;
    FILE *fp_debug;
    char FILE_name[128];
    double residual_norm;
    double *du;

    global.count = 0.;
    init_field();
    
    time_increment = option.Delta_time;
    time = time_increment;
    for(int time_step = 0; time_step < option.N_timestep; time_step++){
        for(int iteration_step = 0; iteration_step < 1000; iteration_step++){
            
            update_field_and_internal_forces();
            update_external_force(time_step);
            residual_norm = calc_global_force_residual_norm();
            if(residual_norm <= option.NR_tol) break;
            generate_coefficient_matrix();

            
            ImposeDirichretResidual(iteration_step + 1);
            ImposeDirichletTangentialMatrix();

            #if 0
            snprintf(FILE_name, 128,"Data_Files_Output/debag%d.dat", iteration_step);
            fp_debug = fopen(FILE_name,"w");
            for(int i = 0; i < global.subdomain.N_point; i++){
                for(int j = 0; j < 3; j++){
                    fprintf(fp_debug, "%+4.3e  ", global.subdomain.global_residual_force[i*3+j]);
                }
                fprintf(fp_debug, "\n");
            }
            fclose(fp_debug);
            #endif

            //求解用の変数ベクトルを用意
            if((du = (double *)calloc(option.dim * global.subdomain.N_point, sizeof(double))) == NULL){
                printf("Error:du's memory is not enough\n");
            }
            solver_LU_decomposition(global.subdomain.Global_K, du, global.subdomain.global_residual_force, option.dim * global.subdomain.N_point);

            for(int i = 0; i < global.subdomain.N_point; i++)
                for(int j = 0; j < option.dim; j++)
                    global.subdomain.displacement_increment[i][j] += du[option.dim * i + j];
            free(du);

            update_nodal_displacement_increment();

            printf("error norm : %+15.14e\n", residual_norm);

            #if 0
            snprintf(FILE_name, 128,"Data_Files_Output/debag%d.dat", iteration_step);
            fp_debug = fopen(FILE_name,"w");
            for(int i = 0; i < global.subdomain.N_point; i++){
                for(int j = 0; j < 3; j++){
                    fprintf(fp_debug, "%+15.14e  ", global.subdomain.global_residual_force[i*3+j]);
                }
                fprintf(fp_debug, "\n");
            }
            fclose(fp_debug);
            #endif
            
        }
    }
    
    break_field();
}