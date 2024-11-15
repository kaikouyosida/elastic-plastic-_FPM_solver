#pragma warning(disable: 4100) // 引数が未使用の場合
#pragma warning(disable: 4189) // ローカル変数が未使用の場合
#pragma warning(disable: 4996) //fopenの警告番号

#include<stdio.h>
#include<stdlib.h>
#include"model.h"
#include"type.h"
#include"fpm.h"
#include"Output_data.h"
#include"matrix.h"
#include"stress.h"

Global global;
Option option;
SS_CURVE ss_curve;

#define NON_LINEAR_SOLVER 1
#define INFINITESIMAL 2
#define LINEAR_SOLVER 3

int main(){
    
    Init_model();
    option.solver_type = 1;

        
    // FILE *fp_extract;
    // double stress[10000][6];
    // double equivalent_stress = 0;
    
    // fp_extract = fopen("Data_Files_Output/Output_cauchy_stress_time19.dat", "r");
    // if(fp_extract == NULL){
    //     printf("File is not open\n");
    //     exit(-1);
    // }
    // fscanf(fp_extract, "%*[^\n]\n");
    
    // for(int i = 0; i < global.subdomain.N_point; i++){
    //     fscanf(fp_extract, "%*d");
    //     for(int j = 0; j < 6; j++){
    //         fscanf(fp_extract, "%lf", &stress[i][j]);
    //     }
    //     fscanf(fp_extract, "\n");
    // }
    // fclose(fp_extract);
    // fp_extract = fopen("D.dat","w");
    // for(int i = 0; i < global.subdomain.N_node; i++){
    //     double equivalent_stress = 0;
    //     int N_ar_node = global.subdomain.ar_node_offset[i+1] - global.subdomain.ar_node_offset[i];
    //     for(int j = 0; j < N_ar_node; j++){
    //         equivalent_stress += calc_equivalent_stress(stress[global.subdomain.ar_node[global.subdomain.ar_node_offset[i]+j]]);
    //     }
    //     equivalent_stress /= N_ar_node;
    //     fprintf(fp_extract, "%+15.14e\n", equivalent_stress);
    // }
    // fclose(fp_extract);
    // exit(0);

    
    switch(option.solver_type){
        case NON_LINEAR_SOLVER:
            analize_by_NewtonRaphson();
        break;
        case INFINITESIMAL:
            infinitesimal_analization();
        break;
        case LINEAR_SOLVER:
            Linear_analization();
        break;
    }


    
    //nit_model()で確保したメモリの開放
    break_model_memory();

    printf("complete!\n");
    return 0;
}