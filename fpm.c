#pragma warning(disable: 4100) // 引数が未使用の場合
#pragma warning(disable: 4189) // ローカル変数が未使用の場合
#pragma warning(disable: 4996) //fopenの警告番号

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
#include"MKL_solver.h"

extern Global global;
extern Option option;

int flag = 0;
int count = 0;



void analize_by_NewtonRaphson(){
    FILE *fp_debug;         //デバッグ用のファイルポインタ
    FILE *fp_residual;
    char FILE_name[128];    //デバッグ用の文字配列
    double residual_norm;   //残差ノルム（収束判定）
    double u_norm;           //変位修正量ノルム（収束判定）
    double error_norm;      //ひずみエネルギノルムの相対誤差
    double *du;             //求解用の変位増分ベクトル
    double *K_u;            //求解用の接線剛性マトリクス
    double *r;              //求解用の残差ベクトル
    double *current_point_xyz;      //現在配置のポイント
    double *Initial_point_XYZ;      //各ステップの現配置のポイント座標

    global.count = 0.;
    //変数のメモリ確保＋初期値を格納
    init_field();

    for(int i = 0; i < global.subdomain.N_point; i++){
        if(fabs(sqrt(global.subdomain.point_XYZ[3*i]* global.subdomain.point_XYZ[3*i] + global.subdomain.point_XYZ[3*i+1] * global.subdomain.point_XYZ[3*i+1]) - 1.0) < 1.0e-2
         && global.subdomain.point_XYZ[3*i+1] < 3.0 - 1.0e-5 && global.subdomain.point_XYZ[3*i+1] > 1.0e-5){
            global.subdomain.point_XYZ[3*i] += 0.01;
            global.subdomain.point_XYZ[3*i+1] += 0.01;
            global.subdomain.point_XYZ[3*i+2] -= 0.01;
        }else if(global.subdomain.point_XYZ[3*i+1] < 0.5 && global.subdomain.point_XYZ[3*i+1] >1.0e-5){
            global.subdomain.point_XYZ[3*i+2] -= 0.01;
        }
    }
    
    for(int time_step = 0; time_step < option.N_timestep; time_step++){
        
        option.time = option.Delta_time * (time_step + 1);
        option.time_old = option.Delta_time * time_step;
        u_norm = 100.0;

        #if 1
        snprintf(FILE_name, 128, "debug_for_residual/Residual_parameter%d.dat", time_step);
        fp_residual = fopen(FILE_name, "w");
        if(fp_residual == NULL){
            printf("residual file is not open\n");
            exit(-1);
        }
        fprintf(fp_residual, "iteration         /     error norm        /       u_norm\n");
        #endif

        for(int iteration_step = 0; iteration_step < 1000; iteration_step++){   //反復計算が１０００回を超えたら強制終了

            //変形勾配テンソル、応力、歪みを更新＋内力ベクトルの更新
            update_field_and_internal_forces();
            printf("debug1\n");
            //外力ベクトルの更新
            update_external_force(time_step);
            printf("debug2\n");
            //係数マトリクスの更新
            generate_coefficient_matrix();
            printf("debug3\n");
            //残差ベクトルの更新＋収束判定パラメータの更新
            residual_norm = calc_global_force_residual_norm(iteration_step);

            fprintf(fp_residual, "%5d  %+15.14e %+15.14e\n", iteration_step, residual_norm, u_norm);
            printf("debug4\n");
            printf("Error:%+15.14e %+15.14e\n", residual_norm, u_norm);
            printf("debug5\n");
            //収束判定（residual_normが閾値を超えたら反復計算を終了)
            #if 0
            if(residual_norm < option.NR_tol || u_norm < option.NR_tol){
                if(residual_norm < option.NR_tol){
                    printf("Step %d: %d time: residual norm %+15.14e\n", time_step, iteration_step, residual_norm);
                    break;
                }else if(u_norm < option.NR_tol){
                    printf("Step %d: %d time: u_norm %+15.14e\n", time_step, iteration_step, u_norm);
                    break;
                }
            }
            #endif

            if(u_norm < option.NR_tol){
                printf("Step %d: %d time: u_norm %+15.14e\n", time_step, iteration_step, u_norm);
                break;
            }
            
            //係数マトリクスにディリクレ境界条件を付与
            ImposeDirichletTangentialMatrix();
            printf("debug6\n");
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
            printf("debug7\n");
            #if 0
            int rank = calculate_rank(K_u, solver_DoF);
            printf("rank = %7d\n", rank);
            #endif
    
            //連立一次方程式を求解
            //solver_LU_decomposition(K_u, du, r, solver_DoF);
            #if 1
            //LU分解で疎行列の連立一次方程式を計算
            int NNZ = 0;
            double *a;
            int *ia, *ja;

            for(int i = 0; i < solver_DoF; i++)
                for(int j = 0; j < solver_DoF; j++)
                    if(K_u[solver_DoF *  i + j] != 0)
                        NNZ++;

            if((a = (double *)calloc(NNZ, sizeof(double))) == NULL){
                printf("Error:a's memory is not enough\n");
                exit(-1);
            }    
            if((ia = (int *)calloc(solver_DoF + 1, sizeof(int))) == NULL){
                printf("Error:ja's memory is not enough\n");
                exit(-1);
            }
            if((ja = (int *)calloc(NNZ, sizeof(int))) == NULL){
                printf("Error:ia's memory is not enough\n");
                exit(-1);
            }
            
            count = 0;
            ia[0] = 0;
            for(int i = 0; i < solver_DoF; i++){
                for(int j = 0; j < solver_DoF; j++){
                    if(K_u[solver_DoF * i + j] != 0){
                        a[count] = K_u[solver_DoF * i + j];
                        ja[count] = j;
                        count++;
                    }
                }
                ia[i+1] = count;
            }

            Paradiso(solver_DoF, NNZ, a, ia, ja, r, du);
            #endif
            printf("debug8\n");
            //ポイントの変位修正ベクトルの値をもとに変位増分を更新.
            if((current_point_xyz = (double *)calloc(option.dim * global.subdomain.N_point, sizeof(double))) == NULL){
                printf("current_point_xyz's memory is not enough\n");
                exit(-1);
            }
            if((Initial_point_XYZ = (double *)calloc(option.dim * global.subdomain.N_point, sizeof(double))) == NULL){
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
            printf("debug9\n");
            u_norm = 0;
            for(int i = 0; i < solver_DoF; i++){
                u_norm += du[i] * du[i];
            }

            for(int i = 0; i < global.subdomain.N_point; i++){
                for(int j = 0; j < option.dim; j++){
                    current_point_xyz[option.dim * i + j] = global.subdomain.point_XYZ[option.dim * i + j]
                                                            + global.subdomain.displacement[i][j]
                                                            + global.subdomain.displacement_increment[i][j];
                }
            }

            for(int i = 0; i < global.subdomain.N_point; i++)
                for(int j = 0; j < option.dim; j++)
                    Initial_point_XYZ[option.dim * i + j] = global.subdomain.point_XYZ[option.dim * i + j] + global.subdomain.displacement[i][j];

            //ノード変位の増分を更新
            update_nodal_displacement_by_inital_NT(Initial_point_XYZ);
            //update_nodal_displacement_increment(current_point_xyz, Initial_point_XYZ);
            printf("debug10\n");
            free(Initial_point_XYZ);
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
        
        increment_field();

        whether_points_is_in_the_subdomain();

        if((time_step + 1) % option.time_output == 0)
            Output_data(time_step);

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
    double u_norm;

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
        u_norm = 100.0;
        option.time = option.Delta_time * (time_step + 1);
        option.time_old = option.Delta_time * time_step;
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
            printf("Error:%+15.14e %+15.14e\n", residual_norm, u_norm);
            

            if(u_norm <= option.NR_tol){
                printf("Step %d: %d time: u_norm %+15.14e\n", time_step, iteration, u_norm);
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
            //solver_LU_decomposition(K_u, du, r, solver_DoF);
            #if 1
            //LU分解で疎行列の連立一次方程式を計算
            int NNZ = 0;
            double *a;
            int *ia, *ja;

            for(int i = 0; i < solver_DoF; i++)
                for(int j = 0; j < solver_DoF; j++)
                    if(K_u[solver_DoF *  i + j] != 0)
                        NNZ++;

            if((a = (double *)calloc(NNZ, sizeof(double))) == NULL){
                printf("Error:a's memory is not enough\n");
                exit(-1);
            }    
            if((ia = (int *)calloc(solver_DoF + 1, sizeof(int))) == NULL){
                printf("Error:ja's memory is not enough\n");
                exit(-1);
            }
            if((ja = (int *)calloc(NNZ, sizeof(int))) == NULL){
                printf("Error:ia's memory is not enough\n");
                exit(-1);
            }
            
            count = 0;
            ia[0] = 0;
            for(int i = 0; i < solver_DoF; i++){
                for(int j = 0; j < solver_DoF; j++){
                    if(K_u[solver_DoF * i + j] != 0){
                        a[count] = K_u[solver_DoF * i + j];
                        ja[count] = j;
                        count++;
                    }
                }
                ia[i+1] = count;
            }

            Paradiso(solver_DoF, NNZ, a, ia, ja, r, du);
            #endif
            printf("Solveed!\n");

            // ポイント変位の増分を更新
            update_point_displaecment_increment(du);
            printf("Displacement increment updated!\n");
            
            u_norm = 0;
            for(int i = 0; i < solver_DoF; i++){
                u_norm += du[i] * du[i];
            }

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