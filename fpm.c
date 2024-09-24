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
    FILE *fp_residual;
    char FILE_name[128];    //デバッグ用の文字配列
    double residual_norm;   //残差ノルム（収束判定）
    double du_norm = 0;   //変位修正量ノルム（収束判定）
    double error_norm;      //ひずみエネルギノルムの相対誤差
    double *du;             //求解用の変位増分ベクトル
    double *K_u;            //求解用の接線剛性マトリクス
    double *r;              //求解用の残差ベクトル
    double *current_point_xyz;      //現在配置のポイント

    global.count = 0.;
    //変数のメモリ確保＋初期値を格納
    init_field();

    for(int time_step = 0; time_step < option.N_timestep; time_step++){

        option.time = option.Delta_time * (time_step + 1);
        option.time_old = option.Delta_time * time_step;
        #if 1
        snprintf(FILE_name, 128, "debug_for_residual/Residual_parameter%d.dat", time_step);
        fp_residual = fopen(FILE_name, "w");
        if(fp_residual == NULL){
            printf("residual file is not open\n");
            exit(-1);
        }
        fprintf(fp_residual, "iteration         /     error norm\n");
        #endif

        for(int iteration_step = 0; iteration_step < 1000; iteration_step++){   //反復計算が１０００回を超えたら強制終了
            
            //変形勾配テンソル、応力、歪みを更新＋内力ベクトルの更新
            update_field_and_internal_forces();

            //外力ベクトルの更新
            update_external_force(time_step);

            //係数マトリクスの更新
            generate_coefficient_matrix();

            //残差ベクトルの更新＋収束判定パラメータの更新
            residual_norm = calc_global_force_residual_norm(iteration_step);
            fprintf(fp_residual, "%5d   %+15.14e\n", iteration_step, residual_norm);
            printf("Error:%+15.14e\n", residual_norm);

            //収束判定（residual_normが閾値を超えたら反復計算を終了)
            #if 1
            if(residual_norm < option.NR_tol){
                printf("Step %d: %d time: residual norm %+15.14e\n", time_step, iteration_step, residual_norm);
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

            //連立一次方程式を求解
            solver_LU_decomposition(K_u, du, r, solver_DoF);


            //ポイントの変位修正ベクトルの値をもとに変位増分を更新.
            if((current_point_xyz = (double *)calloc(option.dim * global.subdomain.N_point, sizeof(double))) == NULL){
                printf("current_point_xyz's memory is not enough\n");
                exit(-1);
            }
            for(int i = 0; i < global.subdomain.N_point; i++){
                for(int j = 0; j < option.dim; j++){
                    global.subdomain.previous_displacement_increment[i][j] = global.subdomain.displacement_increment[i][j];
                }
            }

            // ポイント変位の増分を更新
            update_point_displaecment_increment(du);
            
            double deformation_norm = incremental_deformation_norm();

            //printf("%+15.14e\n", deformation_norm);
            #if 0
            snprintf(FILE_name, 128, "displacement_increment/displacement_increment%d.dat", iteration_step);
            fp_debug = fopen(FILE_name, "w");
            if(fp_debug == NULL){
                printf("Error:memory is not enough\n");
                exit(-1);
            }
            for(int i = 0; i < global.subdomain.N_point; i++){
                for(int j = 0; j < option.dim; j++){
                    fprintf(fp_debug, "%+8.7e   ", global.subdomain.displacement_increment[i][j]);
                }
                fprintf(fp_debug, "\n");
            }
            fclose(fp_debug);
            #endif

            for(int i = 0; i < global.subdomain.N_point; i++){
                for(int j = 0; j < option.dim; j++){
                    current_point_xyz[option.dim * i + j] = global.subdomain.point_XYZ[option.dim * i + j]
                                                            + global.subdomain.displacement[i][j]
                                                            + global.subdomain.displacement_increment[i][j];
                }
            }
            

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
        fclose(fp_residual);

        if((time_step + 1) % option.time_output == 0)
            Output_data(time_step);
        
        increment_field();
        if((time_step + 1) % option.time_output == 0){

            //出力用の節点変位の計算
            update_nodal_coordinate();

            //paraviewデータの出力
            paraview_node_data(time_step);
            
        }
            
    }
    
    break_field();          //変数のメモリを開放
}

void infinitesimal_analization(){
    const int max_iteration_count = 1000;
    double residual_norm; 
    double *du;
    double *r;
    double *K_u;

    FILE *fp_debug;
    char File_name[128];
    FILE *fp_residual;

    init_field();
    
    for(int time_step = 0; time_step < option.N_timestep; time_step++){
       #if 1
        snprintf(File_name, 128, "debug_for_residual/Residual_parameter%d.dat", time_step);
        fp_residual = fopen(File_name, "w");
        if(fp_residual == NULL){
            printf("residual file is not open\n");
            exit(-1);
        }
        fprintf(fp_residual, "iteration         /     error norm\n");
        #endif
        
        for(int iteration = 0; iteration < max_iteration_count; iteration++){
            //内力ベクトルの更新
            update_field_and_internal_infinitesimal();
            printf("Internal force updated!\n");

            //外力ベクトルの更新
            update_external_force(time_step);
            printf("External force updated!\n");

            //接線剛性マトリクスの計算
            generate_coefficient_matrix();
            printf("Coefficient matrix force updated!\n");

            //収束判定
            residual_norm = calc_global_force_residual_norm(time_step);
            fprintf(fp_residual, "%5d   %+15.14e\n", iteration, residual_norm);
            printf("Residual force updated!\n");
            printf("%5d  %+15.14e\n", iteration, residual_norm);
            

            if(residual_norm <= option.NR_tol){
                printf("time: %d  iteration: %d => norm: %+15.14e\n", time_step, iteration, residual_norm);
                break;
            }
            
            //係数マトリクスにディリクレ境界条件を付与
            ImposeDirichletTangentialMatrix();
            printf("Linear system updated!\n");

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
            printf("Matrix assemble updated!\n");

            //連立一次方程式を求解
            solver_LU_decomposition(K_u, du, r, solver_DoF);
            printf("Solveed!\n");

            // ポイント変位の増分を更新
            update_point_displaecment_increment(du);
            printf("Displacement increment updated!\n");

            free(K_u);
            free(r);
            free(du);
            
        }

        fclose(fp_residual);
        increment_field();

        if((time_step + 1) % option.time_output == 0)
            Output_data(time_step);

        if((time_step + 1) % option.time_output == 0){

            //出力用の節点変位の計算
            update_nodal_coordinate();

            //paraviewデータの出力
            paraview_node_data(time_step);
        }
        #if 0
        if(time_step == 0){
            fp_debug = fopen("Yield_stress_vs_plastic_strain_elemtn_1.dat", "w");
            if(fp_debug == NULL){
                printf("Error:File is not open\n");
                exit(-1);
            }
            fprintf(fp_debug, "plastic_strain   /   yiels stress\n");
        }
        fprintf(fp_debug, "%+15.14e       %+15.14e\n", global.subdomain.equivalent_plastic_strains[1], global.subdomain.yield_stresses[1]);
        if(time_step == option.N_timestep-1)
            fclose(fp_debug);
        #endif
    }
   
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