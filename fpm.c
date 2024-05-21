#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"type.h"
#include"fpm.h"
#include"field.h"
#include"coefficient_matrix.h"
#include"internal_force.h"
#include"external_force.h"
#include"ImposeDirichretCondition.h"
#include"LU_decomposition.h"
#include"Output.h"
#include"Output_data.h"

extern Global global;
extern Option option;

void analize_by_NewtonRapdon(){
    FILE *fp_debug;
    char FILE_name[128];
    double residual_norm;
    double *du;

    global.count = 0.;
    init_field();           //変数を初期化
    for(int time_step = 0; time_step < option.N_timestep; time_step++){
        option.time = option.Delta_time * (time_step + 1);
        option.time_old = option.Delta_time * time_step;
        for(int iteration_step = 0; iteration_step < 1000; iteration_step++){   //反復計算が１０００回を超えたら強制終了
            update_field_and_internal_forces();
            update_external_force(time_step);
            residual_norm = calc_global_force_residual_norm(iteration_step);
            if(residual_norm <= option.NR_tol && iteration_step != 0){
                printf("Step %d: %d time: residual norm %+15.14e\n", time_step, iteration_step, residual_norm);
                break;
            }
            generate_coefficient_matrix();

            #if 0
            if(iteration_step == 1){
                fp_debug = fopen("coefficient_for_debug.dat", "w");
                if(fp_debug == NULL)
                    printf("FILE is not enough\n");
                for(int i = 0; i < 3*global.subdomain.N_point; i++){
                    for(int j = 0; j < 3*global.subdomain.N_point; j++){
                        fprintf(fp_debug, "%+8.7e   ", global.subdomain.Global_K[3*global.subdomain.N_point*i+j]);
                    }
                    fprintf(fp_debug, "\n");
                }
            }
            if(iteration_step > 1){
                fp_debug = fopen("coefficient_for_debug.dat", "r");
                if(fp_debug == NULL)
                    printf("File is not open\n");
                for(int i = 0; i < 3*global.subdomain.N_point; i++){
                    for(int j = 0; j < 3*global.subdomain.N_point; j++){
                        fscanf(fp_debug, "%lf   ", &global.subdomain.Global_K[3*global.subdomain.N_point*i+j]);
                    }
                    fgetc(fp_debug);
                }
                fclose(fp_debug);
            }
            #endif
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
            
            #if 1
            snprintf(FILE_name, 128,"Data_Files_Output/debag%d.dat", iteration_step);
            fp_debug = fopen(FILE_name,"w");
            fprintf(fp_debug, "point        /displacement           x           y           z\n");
            for(int i = 0; i < global.subdomain.N_point; i++){
                fprintf(fp_debug, "%5d  ", i);
                for(int j = 0; j < option.dim; j++){
                    fprintf(fp_debug, "%+15.14e  ", du[option.dim * i + j]);
                }
                fprintf(fp_debug, "\n");
            }
            fclose(fp_debug);
            #endif
            free(du);
            #if 0
            double u_norm = 0.;
            for(int i = 0; i < global.subdomain.N_point; i++){
                for(int j = 0; j < option.dim; j++){
                    u_norm += du[option.dim * i + j] * du[option.dim * i + j];
                }
            }
            u_norm = sqrt(u_norm);
            printf("%+15.14e\n", u_norm);
            #endif

            update_nodal_displacement_increment();
            //for(int i = 0 ; i < global.subdomain.N_node; i++){
                //for(int j = 0; j < 3; j++){
                    //printf("%+15.14e   ", global.subdomain.nodal_displacement_increments[i][j] + global.subdomain.node_XYZ[3*i+j]);
                //}
                //printf("\n");
            //}
            //for(int i = 0 ; i < global.subdomain.N_point; i++){
                //for(int j = 0; j < 3; j++){
                    //printf("%+15.14e   ", global.subdomain.displacement_increment[i][j]);
                //}
                //printf("\n");
            //}
            

            printf("%5d   %+15.14e\n", iteration_step+1, residual_norm);
            if(iteration_step == 1000){
                printf("Iteration is not converged\n");
                exit(-1);
            }
            
        }
        Output_data(time_step);
        increment_field();
        paraview_node_data(time_step);
    }
    
    break_field();          //変数のメモリを開放
}

void Linear_analization(){
    FILE *fp_debug;
    double *f_ext;
    double *du;

    option.time_old = 0.;
    option.time = 1.0;

    init_field();
    update_external_force(0);
    global.buf = calc_global_force_residual_norm(0);
    generate_coefficient_linear();
    ImposeDirichletTangentialMatrix();
    
    //求解用の変数ベクトルを用意
     if((du = (double *)calloc(option.dim * global.subdomain.N_point, sizeof(double))) == NULL){
        printf("Error: du's Memory is not enough\n");
        exit(-1);
    }
    printf("Now solving!!\n");
    solver_LU_decomposition(global.subdomain.Global_K, du, global.subdomain.global_residual_force, option.dim * global.subdomain.N_point);
    #if 0
        fp_debug = fopen("debug.dat", "w");
        for(int i = 0; i < global.subdomain.N_point; i++){
            for(int j = 0; j < 3; j++){
                fprintf(fp_debug, "%+15.14e  ", global.subdomain.global_residual_force[i*3+j]);
            }
            fprintf(fp_debug, "\n");
        }
        fclose(fp_debug);
        exit(-1);
    #endif
    printf("status1\n");
    Output_Linear_strain_data(du);
        

    free(du);
    break_field();
}