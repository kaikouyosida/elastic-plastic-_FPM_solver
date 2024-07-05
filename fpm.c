#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"type.h"
#include"fpm.h"
#include"scalar.h"
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

int flag = 0;
int count = 0;

void analize_by_NewtonRaphson(){
    FILE *fp_debug;         //デバッグ用のファイルポインタ
    char FILE_name[128];    //デバッグ用の文字配列
    double error_old = 100.0;       //１反復前の残差ノルム
    double residual_norm;   //残差ノルム（収束判定）
    double du_norm = 0;   //変位修正量ノルム（収束判定）
    
    double *du;             //求解用の変位増分ベクトル
    double *K_u;            //求解用の接線剛性マトリクス
    double *r;              //求解用の残差ベクトル
    double *current_point_xyz;

    global.count = 0.;
    //変数のメモリ確保＋初期値を格納
    init_field();

    for(int time_step = 0; time_step < option.N_timestep; time_step++){
        int stop_count = 0;

        option.time = option.Delta_time * (time_step + 1);
        option.time_old = option.Delta_time * time_step;
        
        for(int iteration_step = 0; iteration_step < 1000; iteration_step++){   //反復計算が１０００回を超えたら強制終了

            //変形勾配テンソル、応力、歪みを更新＋内力ベクトルの更新
            update_field_and_internal_forces();

            //外力ベクトルの更新
            update_external_force(time_step);

            //係数マトリクスの更新
            generate_coefficient_matrix();

            #if 0
            if(iteration_step == 0){
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
            if(iteration_step > 0){
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
            //残差ベクトルの更新＋収束判定パラメータの更新
            residual_norm = calc_global_force_residual_norm(iteration_step);
            if(residual_norm <= option.NR_tol){// && iteration_step != 0){
                printf("Step %d: %d time: residual norm %+15.14e\n", time_step, iteration_step, residual_norm);
                break;
            }
            printf("%5d   %+15.14e\n", iteration_step+1, residual_norm);
            #if 0
            //残差が著しく減少しない場合
            if(0.5 < residual_norm / error_old){
                stop_count++;
            }else{
                error_old = residual_norm;
            }
            if(stop_count == option.stop_count){
                printf("Iteration stop\n");
                break;
            }
            #endif
            //係数マトリクスにディリクレ境界条件を付与
            ImposeDirichletTangentialMatrix();
            

            //求解用の変数ベクトルと係数マトリクスを用意.LU分解で連立一次方程式を求解.
            int solver_DoF = global.subdomain.N_point * option.dim - global.bc.N_D_DoF;
            
            if((du = (double *)calloc(solver_DoF, sizeof(double))) == NULL){
                printf("Error:du's memory is not enough\n");
                exit(-1);
            }
            if((r = (double *)calloc(solver_DoF, sizeof(double))) == NULL){
                printf("Error:r's memory is not enough\n");
                exit(-1);
            }
            if((K_u = (double *)calloc(solver_DoF * solver_DoF, sizeof(double))) == NULL){
                printf("Error:K_u's memory is not enough\n");
                exit(-1);
            }
            assemble_matrix_and_vector_for_Dirichlet(K_u, r);
            //printf("Now solving!!\n");
            solver_LU_decomposition(K_u, du, r, solver_DoF);
            
            //ポイントの変位修正ベクトルの値をもとに変位増分を更新.
            if((current_point_xyz = (double *)calloc(option.dim * global.subdomain.N_point, sizeof(double))) == NULL){
                printf("current_point_xyz's is not neough\n");
                exit(-1);
            }
            for(int i = 0; i < global.subdomain.N_point; i++){
                for(int j = 0; j < option.dim; j++){
                    current_point_xyz[option.dim * i + j] = global.subdomain.point_XYZ[option.dim * i + j]
                                                        + global.subdomain.displacement[i][j]
                                                        + global.subdomain.displacement_increment[i][j];
                }
            }

            // ポイント変位の増分を更新
            update_point_displaecment_increment(current_point_xyz, du);
            //ノード変位の増分を更新
            update_nodal_displacement_increment(current_point_xyz);
            
            free(current_point_xyz);
            free(K_u);
            free(r);
            free(du);
            
            if(iteration_step == 1000){
                printf("Iteration is not converged\n");
                exit(-1);
            }
            
        }

        if((time_step + 1) % option.time_output == 0)
            Output_data(time_step);
        increment_field();
        if((time_step + 1) % option.time_output == 0)
            paraview_node_data(time_step);
    }
    
    break_field();          //変数のメモリを開放
}

void Linear_analization(){
    FILE *fp_debug;
    double *f_ext;
    double *deformation;
    double *du;
    double *K_u;            //求解用の接線剛性マトリクス
    double *r;              //求解用の残差ベクトル

    option.time_old = 0.;
    option.time = 1.0;
    //変数のメモリ確保＋初期値を格納
    init_field();

    //外力ベクトルの更新
    update_external_force(0);

    
    generate_coefficient_linear();
    global.buf = calc_global_force_residual_norm(0);
    ImposeDirichletTangentialMatrix();
    
    //求解用の変数ベクトルを用意
     if((du = (double *)calloc(option.dim * global.subdomain.N_point, sizeof(double))) == NULL){
        printf("Error: du's Memory is not enough\n");
        exit(-1);
    }
    printf("Now solving!!\n");
    //求解用の変数ベクトルと係数マトリクスを用意.LU分解で連立一次方程式を求解.
    int solver_DoF = global.subdomain.N_point * option.dim - global.bc.N_D_DoF;

    if((deformation = (double *)calloc(global.subdomain.N_point * option.dim, sizeof(double))) == NULL){
        printf("Error:deforamtion's memory is not enough\n");
        exit(-1);
    }
    if((du = (double *)calloc(solver_DoF, sizeof(double))) == NULL){
        printf("Error:du's memory is not enough\n");
        exit(-1);
    }
    if((r = (double *)calloc(solver_DoF, sizeof(double))) == NULL){
        printf("Error:r's memory is not enough\n");
        exit(-1);
    }
    if((K_u = (double *)calloc(solver_DoF * solver_DoF, sizeof(double))) == NULL){
        printf("Error:K_u's memory is not enough\n");
        exit(-1);
    }
    assemble_matrix_and_vector_for_Dirichlet(K_u, r);
    //printf("Now solving!!\n");
    solver_LU_decomposition(K_u, du, r, solver_DoF);
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
    for(int i = 0; i < global.subdomain.N_point; i++){
        for(int j = 0; j < option.dim; j++){
            for(int k = 0; k < global.bc.N_D_DoF; k++)
                if(option.dim * i + j == global.bc.fixed_dof[k]) flag++;
            if(flag == 0){
                deformation[option.dim * i + j] += du[count];
                count++;
            }
            flag = 0;   
        }
    }
    printf("status1\n");
    Output_Linear_strain_data(deformation);
        
    free(K_u);
    free(r);
    free(du);
    free(deformation);
    break_field();
}