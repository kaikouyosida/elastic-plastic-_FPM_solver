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
            
            ImposeDirichretResidual(time_step + 1);
            ImposeDirichletTangentialMatrix();
            
            //求解用の変数ベクトルを用意
            if((du = (double *)calloc(3 * global.subdomain.N_point, sizeof(double))) == NULL){
                printf("Error:du's memory is not enough\n");
            }
            solver_LU_decomposition(global.subdomain.Global_K, du, global.subdomain.global_residual_force,option.dim * global.subdomain.N_point);

            for(int i = 0; i < global.subdomain.N_point; i++)
                for(int j = 0; j < option.dim; j++)
                    global.subdomain.displacement_increment[i][j] += du[option.dim * i + j];
            free(du);

            update_nodal_displacement_increment();
            printf("error norm : %+15.14e\n", residual_norm);

            snprintf(FILE_name, 128,"debag%d.dat", iteration_step);
            fp_debug = fopen(FILE_name,"w");
            for(int i = 0; i < global.subdomain.N_point; i++){
                for(int j = 0; j < 3; j++){
                    fprintf(fp_debug, "%+15.14e  ", global.subdomain.displacement_increment[i][j]);
                }
                fprintf(fp_debug, "\n");
            }

            

        }
    }
    
    break_field();
}